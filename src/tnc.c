/* tnc : truncated newton bound constrained minimization
   using gradient information, in C */

/*
 * Copyright (c) 2002-2005, Jean-Sebastien Roy (js@jeannot.org)
 *
 * Permission is hereby granted, free of charge, to any person obtaining a
 * copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject to
 * the following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
 * CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
 * TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
 * SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

/*
 * This software is a C implementation of TNBC, a truncated newton minimization
 * package originally developed by Stephen G. Nash in Fortran.
 *
 * The original source code can be found at :
 * http://iris.gmu.edu/~snash/nash/software/software.html
 *
 * Copyright for the original TNBC fortran routines:
 *
 *   TRUNCATED-NEWTON METHOD:  SUBROUTINES
 *     WRITTEN BY:  STEPHEN G. NASH
 *           SCHOOL OF INFORMATION TECHNOLOGY & ENGINEERING
 *           GEORGE MASON UNIVERSITY
 *           FAIRFAX, VA 22030
 */

/*
 * Conversion into C by Elisabeth Nguyen & Jean-Sebastien Roy
 * Modifications by Jean-Sebastien Roy, 2001-2002
 *
 * This code has been modified by L. Garrido and M. Kalmoun to 
 * adapt to the multigrid algorithm. In addiction, the line 
 * search algorithm has been changed by one proposed in the 
 * L-BFGS algorith by Nash since the current line search gave
 * sometimes weird problems in the multigrid approach.
 */

/* This code has been adapted by L. Garrido and M. Kalmoun for the
   purpose of multigrid optical flow. The main change is the removal
   of the line-search algorithm that was included by default. The
   line-search method included in the L-BFGS method has been used,
   see README file. */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "tnc.h"
#include "parameters.h"
#include "proto/tnc_utils.h"

#define MIN(a,b)    ((a)<(b)?(a):(b))

/**
 * Return code strings
 */

static char *tnc_rc_string[14] =
{
  "Memory allocation failed",
  "Invalid parameters (n<0)",
  "Infeasible (low bound > up bound)",
  "Local minima reach (|pg| ~= 0)",
  "Converged (|f_n-f_(n-1)| ~= 0)",
  "Converged (|x_n-x_(n-1)| ~= 0)",
  "Maximum number of function evaluations reached",
  "Linear search failed",
  "All lower bounds are equal to the upper bounds",
  "Unable to progress",
  "User requested end of minimization",
  "Presmoothing condition satisfied",
  "Correction is not a descent"
};

/**
 *
 * This routine solves the optimization problem
 *
 *   minimize   f(x)
 *     x
 *   subject to   low <= x <= up
 *
 * where x is a vector of n real variables. The method used is
 * a truncated-newton algorithm (see "newton-type minimization via
 * the lanczos algorithm" by s.g. nash (technical report 378, math.
 * the lanczos method" by s.g. nash (siam j. numer. anal. 21 (1984),
 * pp. 770-778).  this algorithm finds a local minimum of f(x). It does
 * not assume that the function f is convex (and so cannot guarantee a
 * global solution), but does assume that the function is bounded below.
 * it can solve problems having any number of variables, but it is
 * especially useful when the number of variables (n) is large.
 *
 */

int tnc(int n, double x[], double *f, double g[], 
    tnc_function *function, void *state, double low[], double up[], double scale[], double offset[],
    int messages, int maxiter, int maxCGit, int maxnfeval, double eta, double stepmx,
    double accuracy, double fmin, double ftol, double xtol, double pgtol,
    double rescale, int *nfeval)
{
  int frc __attribute__ ((unused));
  int rc, i, nc, nfeval_local,
      free_low = TNC_FALSE, free_up = TNC_FALSE,
      free_g = TNC_FALSE;
  double *xscale = NULL, fscale, epsmch, rteps, *xoffset = NULL;

  if(nfeval==NULL)
  {
    /* Ignore nfeval */
    nfeval = &nfeval_local;
  }
  *nfeval = 0;

  /* Check for errors in the input parameters */
  if (n == 0)
  {
    rc = TNC_CONSTANT;
    goto cleanup;
  }

  if (n < 0)
  {
    rc = TNC_EINVAL;
    goto cleanup;
  }

  /* Check bounds arrays */
  if (low == NULL)
  {
    low = malloc(n*sizeof(*low));
    if (low == NULL)
    {
      rc = TNC_ENOMEM;
      goto cleanup;
    }
    free_low = TNC_TRUE;
    for (i = 0 ; i < n ; i++) low[i] = -HUGE_VAL;
  }
  if (up == NULL)
  {
    up = malloc(n*sizeof(*up));
    if (up == NULL)
    {
      rc = TNC_ENOMEM;
      goto cleanup;
    }
    free_up = TNC_TRUE;
    for (i = 0 ; i < n ; i++) up[i] = HUGE_VAL;
  }

  /* Coherency check */
  for (i = 0 ; i < n ; i++)
  {
    if (low[i] > up [i])
    {
      rc = TNC_INFEASIBLE;
      goto cleanup;
    }
  }

  /* Coerce x into bounds */
  tnc_coercex(n, x, low, up);

  if (maxnfeval < 1)
  {
    rc = TNC_MAXFUN;
    goto cleanup;
  }

  /* Allocate g if necessary */
  if(g == NULL)
  {
    g = malloc(n*sizeof(*g));
    if (g == NULL)
    {
      rc = TNC_ENOMEM;
      goto cleanup;
    }
    free_g = TNC_TRUE;
  }

  /* Initial function evaluation */

  cac_tnc_set_uniform_sampling(1, state); 
  frc = function(x, f, g, state);
  (*nfeval) ++;
  if (frc)
  {
    rc = TNC_USERABORT;
    goto cleanup;
  }

  /* Constant problem ? */
  for (nc = 0, i = 0 ; i < n ; i++)
    if ((low[i] == up[i]) || (scale != NULL && scale[i] == 0.0))
      nc ++;

  if (nc == n)
  {
    rc = TNC_CONSTANT;
    goto cleanup;
  }

  /* Scaling parameters */
  xscale = malloc(sizeof(*xscale)*n);
  if (xscale == NULL)
  {
    rc = TNC_ENOMEM;
    goto cleanup;
  }
  xoffset = malloc(sizeof(*xoffset)*n);
  if (xoffset == NULL)
  {
    rc = TNC_ENOMEM;
    goto cleanup;
  }
  fscale = 1.0;

  for (i = 0 ; i < n ; i++)
  {
    if (scale != NULL)
    {
      xscale[i] = fabs(scale[i]);
      if (xscale[i] == 0.0)
	xoffset[i] = low[i] = up[i] = x[i];
    }
    else if (low[i] != -HUGE_VAL && up[i] != HUGE_VAL)
    {
      xscale[i] = 1.0;
      xoffset[i] = 0.0;
    }
    else
    {
      xoffset[i] = 0.0;
      xscale[i]  = 1.0;
    }
    if (offset != NULL)
      xoffset[i] = offset[i];
  }

  /* Default values for parameters */
  epsmch = tnc_mchpr1();
  rteps = sqrt(epsmch);

  if (stepmx < rteps * 10.0) stepmx = 1e20; //1.0e1;
  if (eta < 0.0 || eta >= 1.0) eta = 0.25;
  if (rescale < 0) rescale = 1.3;
  if (maxCGit < 0)
  {
    maxCGit = n / 2;
    if (maxCGit < 1) maxCGit = 1;
    else if (maxCGit > 50) maxCGit = 50;
  }
  if (maxCGit > n) maxCGit = n;
  if (accuracy <= epsmch) accuracy = rteps;
  if (ftol < 0.0) ftol = accuracy;
  if (pgtol < 0.0) pgtol = 1e-2 * sqrt(accuracy);
  if (xtol < 0.0) xtol = rteps;

  /* Optimisation */
  rc = tnc_minimize(n, x, f, g, function, state,
      xscale, xoffset, &fscale, low, up, messages,
      maxCGit, maxiter, maxnfeval, nfeval, eta, stepmx, accuracy, fmin, ftol, xtol, pgtol,
      rescale);

cleanup:
  if (messages & TNC_MSG_EXIT)
    fprintf(stderr, "tnc: %s\n", tnc_rc_string[rc - TNC_MINRC]);

  if (xscale) free(xscale);
  if (free_low) free(low);
  if (free_up) free(up);
  if (free_g) free(g);
  if (xoffset) free(xoffset);

  return rc;
}

/**
 *
 * Coerce x into bounds 
 *
 */

void tnc_coercex(int n, double x[], double low[], double up[])
{
  int i;

  for (i = 0 ; i < n ; i++)
  {
    if (x[i]<low[i]) x[i] = low[i];
    else if (x[i]>up[i]) x[i] = up[i];
  }
}

/**
 *
 * Unscale x 
 *
 */

void tnc_unscalex(int n, double x[], double xscale[], double xoffset[])
{
  int i;
  for (i = 0 ; i < n ; i++)
    x[i] = x[i]*xscale[i]+xoffset[i];
}

/**
 *
 * Scale x 
 *
 */

void tnc_scalex(int n, double x[], double xscale[], double xoffset[])
{
  int i;
  for (i = 0 ; i < n ; i++)
    if (xscale[i]>0.0)
      x[i] = (x[i]-xoffset[i])/xscale[i];
}

/** 
 * 
 * Scale g 
 *
 */

void tnc_scaleg(int n, double g[], double xscale[], double fscale)
{
  int i;
  for (i = 0 ; i < n ; i++)
    g[i] *= xscale[i]*fscale;
}

/**
 *
 * Caculate the pivot vector 
 *
 */
void tnc_setConstraints(int n, double x[], int pivot[], double xscale[],
    double xoffset[], double low[], double up[])
{
  int i;
  double epsmch, xval;

  epsmch = tnc_mchpr1();

  for (i = 0; i < n; i++)
  {
    /* tolerances should be better ajusted */
    if (xscale[i] == 0.0)
    {
      pivot[i] = 2;
    }
    else
    {
      xval = x[i] * xscale[i] + xoffset[i]; 

      if (low[i] != - HUGE_VAL && (xval - low[i] <= epsmch * 10.0 * (fabs(low[i]) + 1.0)))
	pivot[i] = -1;
      else
      {
	if (up[i] != HUGE_VAL && (up[i] - xval <= epsmch * 10.0 * (fabs(up[i]) + 1.0)))
	  pivot[i] = 1;
	else
	  pivot[i] = 0;
      }
    }
  }
}


/**
 * This routine is a bounds-constrained truncated-newton method.
 * the truncated-newton method is preconditioned by a limited-memory
 * quasi-newton method (this preconditioning strategy is developed
 * in this routine) with a further diagonal scaling
 * (see routine diagonalscaling).
 */

tnc_rc tnc_minimize(int n, double x[],
    double *f, double gfull[], tnc_function *function, void *state,
    double xscale[], double xoffset[], double *fscale,
    double low[], double up[], tnc_message messages,
    int maxCGit, int maxiter, int maxnfeval, int *nfeval, double eta, double stepmx,
    double accuracy, double fmin, double ftol, double xtol, double pgtol,
    double rescale)
{
  double fLastReset, difnew, epsmch, epsred, oldgtp,
	 difold, oldf, xnorm, 
	 gnorm, ustpmax, fLastConstraint, spe, yrsr, yksk, gu,
	 *temp = NULL, *sk = NULL, *yk = NULL, *diagb = NULL, *sr = NULL,
	 *yr = NULL, *oldg = NULL, *pk = NULL, *g = NULL, *oldx = NULL;
  double alpha = 0.0; /* Default unused value */
  int i, ninternal, nint_iter, icycle, niter = 0, oldnfeval, *pivot = NULL, frc;
  tnc_logical lreset, newcon, upd1, remcon;
  tnc_rc rc = TNC_ENOMEM; /* Default error */

  eta += 0.0;  /* Not used, but avoid complaining of compiler */
  rescale += 0.0; /* Not used, but avoid complaining of compiler */

  /* Allocate temporary vectors */
  oldg = malloc(sizeof(*oldg)*n);
  if (oldg == NULL) goto cleanup;
  g = malloc(sizeof(*g)*n);
  if (g == NULL) goto cleanup;
  temp = malloc(sizeof(*temp)*n);
  if (temp == NULL) goto cleanup;
  diagb = malloc(sizeof(*diagb)*n);
  if (diagb == NULL) goto cleanup;
  pk = malloc(sizeof(*pk)*n);
  if (pk == NULL) goto cleanup;

  sk = malloc(sizeof(*sk)*n);
  if (sk == NULL) goto cleanup;
  yk = malloc(sizeof(*yk)*n);
  if (yk == NULL) goto cleanup;
  sr = malloc(sizeof(*sr)*n);
  if (sr == NULL) goto cleanup;
  yr = malloc(sizeof(*yr)*n);
  if (yr == NULL) goto cleanup;

  pivot = malloc(sizeof(*pivot)*n);
  if (pivot == NULL) goto cleanup;

  oldx = malloc(sizeof(*oldx)*n);
  if (oldx == NULL) goto cleanup;

  /* Initialize variables */
  epsmch = tnc_mchpr1();

  difnew = 0.0;
  epsred = 0.05;
  upd1 = TNC_TRUE;
  icycle = n - 1;
  newcon = TNC_TRUE;

  /* Uneeded initialisations */
  lreset = TNC_FALSE;
  yrsr = 0.0;
  yksk = 0.0;

  /* Initial scaling */
  tnc_scalex(n, x, xscale, xoffset);
  (*f) *= *fscale;

  /* initial pivot calculation */
  tnc_setConstraints(n, x, pivot, xscale, xoffset, low, up);

  tnc_dcopy1(n, gfull, g);
  tnc_scaleg(n, g, xscale, *fscale);

  /* Test the lagrange multipliers to see if they are non-negative. */
  for (i = 0; i < n; i++)
    if (-pivot[i] * g[i] < 0.0)
      pivot[i] = 0;

  tnc_project(n, g, pivot);

  /* Set initial values to other parameters */
  gnorm = tnc_dnrm21(n, g);

  fLastConstraint = *f; /* Value at last constraint */
  fLastReset = *f; /* Value at last reset */

  tnc_dcopy1(n, gfull, temp);
  tnc_dneg1(n, temp);
  cac_tnc_normalize_gradient(n, temp);

  cac_progress(x, temp, *f, gnorm, 0.0, niter, state);

  /* Set the diagonal of the approximate hessian to unity. */
  for (i = 0; i < n; i++) diagb[i] = 1.0;

  /* We begin iterating ... */
  ninternal = 0;
  nint_iter = 0;

  /* Start of main iterative loop */
  while(TNC_TRUE)
  {
    /* Local minimum test */
    if (tnc_dnrm21(n, g) <= pgtol * (*fscale))
    {
      /* |PG| == 0.0 => local minimum */
      tnc_dcopy1(n, gfull, g);
      tnc_project(n, g, pivot);
      if (messages & TNC_MSG_INFO) 
	fprintf(stderr, "tnc: |pg| = %g -> local minimum\n", tnc_dnrm21(n, g) / (*fscale));
      rc = TNC_LOCALMINIMUM;
      break;
    }

    /* Terminate if more than maxiter iterations have been made */
    if (niter >= maxiter) 
    {
      rc = TNC_MAXFUN;
      break;
    }

    tnc_dcopy1(n, x, temp);
    tnc_project(n, temp, pivot);
    xnorm = tnc_dnrm21(n, temp);
    oldnfeval = *nfeval;

    /* Compute the new search direction by solving the Newon equation.
     * If only a line search has to be performed, get the serach direction
     * and store it in pk */

#if VERBOSE > 1
    cac_tnc_print_msg("\nComputing search direction...\n", state);
#endif
    frc = tnc_direction(pk, diagb, x, g, n, maxCGit, &nint_iter, maxnfeval, nfeval,
	upd1, yksk, yrsr, sk, yk, sr, yr,
	lreset, function, state, xscale, xoffset, *fscale,
	pivot, accuracy, gnorm, xnorm, low, up);

    ninternal += nint_iter;

    if (frc == -1)
    {
      rc = TNC_ENOMEM;
      break;
    }

    if (frc)
    {
      rc = TNC_USERABORT;
      break;
    }

    if (!newcon)
    {
      if (!lreset)
      {
	/* Compute the accumulated step and its corresponding gradient
	   difference. */
	tnc_dxpy1(n, sk, sr);
	tnc_dxpy1(n, yk, yr);
	icycle++;
      }
      else
      {
	/* Initialize the sum of all the changes */
	tnc_dcopy1(n, sk, sr);
	tnc_dcopy1(n, yk, yr);
	fLastReset = *f;
	icycle = 1;
      }
    }

    tnc_dcopy1(n, g, oldg);
    oldf = *f;
    oldgtp = tnc_ddot1(n, pk, g);

    /* Check descent if we are doing a descent */

    if (cac_tnc_get_project_gradient(state))
      cac_tnc_project_gradient_line_vertex_to_center(x, pk, state);

    cac_tnc_normalize_gradient(n, pk);

    tnc_dcopy1(n, gfull, temp);
    tnc_scaleg(n, temp, xscale, *fscale);
    gu = tnc_ddot1(n, temp, pk);

    if (gu > 0.0)
    {
      /* Take the negative of the gradient to descent */

      tnc_dcopy1(n, gfull, pk);
      tnc_scaleg(n, pk, xscale, *fscale);
      tnc_dneg1(n, pk);
      gu = tnc_ddot1(n, temp, pk);
    }

    /* Maximum unconstrained step length */
    ustpmax = stepmx; 

    /* Maximum constrained step length */
    //spe = tnc_stepMax(ustpmax, n, x, pk, pivot, low, up, xscale, xoffset);

    spe = cac_tnc_get_maxstep_linsrch(n, x, pk, state);

#if VERBOSE > 1
    char msg[200];
    sprintf(msg, "Max step = %E\n", spe);
    cac_tnc_print_msg(msg, state); 
#endif

    if (spe > accuracy)
    {
      tnc_ls_rc lsrc;

      /* Set the initial step length */

      alpha = MIN(DEFAULTMOVE, spe);

#if VERBOSE > 1
      sprintf(msg, "Initial alpha = %E\n", alpha);
      cac_tnc_print_msg(msg, state);
#endif

      tnc_dcopy1(n, x, oldx);

#if VERBOSE > 1
      cac_tnc_print_msg("Performing line search... ", state);
#endif

      /* Perform the linear search */
      lsrc = tnc_linearSearch(n, low, up, xscale, xoffset, *fscale, pivot,
	  x, f, gfull, pk, &alpha, spe, gu, xtol, maxnfeval, nfeval, function, state, 0);

#if VERBOSE > 1
      sprintf(msg, "alpha = %E\n\n", alpha);
      cac_tnc_print_msg(msg, state);
#endif

      if (lsrc == LS_ENOMEM)
      {
	rc = TNC_ENOMEM;
	break;
      }

      if (lsrc == LS_USERABORT)
      {
	rc = TNC_USERABORT;
	break;
      }

      if (lsrc == LS_FAIL)
      {
	rc = TNC_LSFAIL;
	break;
      }

      if (lsrc == LS_MAXFUN)
      {
	rc = TNC_MAXFUN;
	break;
      }

      if (alpha >= ustpmax)
      {
	fprintf(stderr, "\nERROR(tnc_minimize): alpha (%e) >= ustpmax (%e). Problem may be badly scaled.\n", alpha, ustpmax);
	exit(1);
      }

      /* If we went up to the maximum constrained step,
       * a new constraint was encountered */

      if (alpha - spe >= -epsmch * 10.0) 
      {
	newcon = TNC_TRUE;
      }
      else
      {
	/* Break if the linear search has failed to find a lower point */
	if (lsrc != LS_OK)
	{
	  rc = TNC_LSFAIL;
	  break;
	}
	newcon = TNC_FALSE;
      }
    }
    else
    {
      /* Maximum constrained step == 0.0 */ 

      newcon = TNC_TRUE;
    }

    if (newcon)
    {
      newcon = tnc_addConstraint(n, x, pk, pivot, low, up, xscale, xoffset);

      if (!newcon)
      {
	if(*nfeval == oldnfeval)
	{
	  rc = TNC_NOPROGRESS;
	  break;
	}
      }

      fLastConstraint = *f;
    }

    niter++;

    /* Set up parameters used in convergence and resetting tests */
    difold = difnew;
    difnew = oldf - *f;

    /* If this is the first iteration of a new cycle, compute the
       percentage reduction factor for the resetting test */
    if (icycle == 1)
    {
      if (difnew > difold * 2.0) epsred += epsred;
      if (difnew < difold * 0.5) epsred *= 0.5;
    }

    tnc_dcopy1(n, gfull, g);
    tnc_scaleg(n, g, xscale, *fscale);

    tnc_dcopy1(n, g, temp);
    tnc_project(n, temp, pivot);
    gnorm = tnc_dnrm21(n, temp);

    /* Reset pivot */
    remcon = tnc_removeConstraint(oldgtp, gnorm, pgtol * (*fscale), *f,
	fLastConstraint, g, pivot, n);

    /* If a constraint is removed */
    if (remcon)
    {
      /* Recalculate gnorm and reset fLastConstraint */
      tnc_dcopy1(n, g, temp);
      tnc_project(n, temp, pivot);
      gnorm = tnc_dnrm21(n, temp);
      fLastConstraint = *f;
    }

    if (!remcon && !newcon)
    {
      if (cac_tnc_get_project_gradient(state)) 
      {	
	if (alpha < 1.0)
	{
	  cac_tnc_print_msg("tnc: project gradient converged. setting project_gradient to zero!\n\n", state);
	  cac_tnc_set_project_gradient(0, state);
	  // maxCGit = 10; TODO: This has not been tested thoroughly. The 
	  // idea is to use conjugate gradient for the second step in order
	  // to in increase efficiency.
	}
      }
      else
      {
	if (alpha < THRSALPHA)
	{
	  cac_tnc_print_msg("tnc: alpha < THRSALPHA -> converged\n", state);
	  rc = TNC_XCONVERGED;
	  break;
	}
	/* No constraint removed & no new constraint : tests for convergence */
	if (fabs(difnew) <= ftol * (*fscale)) 
	{
	  cac_tnc_print_msg("tnc: |fn-fn-1] -> convergence\n", state);
	  rc = TNC_FCONVERGED;
	  break;
	}
	if (alpha * tnc_dnrm21(n, pk) <= xtol)
	{
	  cac_tnc_print_msg("tnc: |xn-xn-1] -> convergence\n", state);
	  rc = TNC_XCONVERGED;
	  break;
	}
      }
    }

    tnc_project(n, g, pivot);

    cac_progress(x, pk, *f, gnorm, alpha, niter, state);

    /* Compute the change in the iterates and the corresponding change in the
       gradients */
    if (!newcon)
    {
      for (i = 0; i < n; i++)
      {
	yk[i] = g[i] - oldg[i];
	sk[i] = alpha * pk[i];
      }

      /* Set up parameters used in updating the preconditioning strategy */
      yksk = tnc_ddot1(n, yk, sk);

      if (icycle == (n - 1) || difnew < epsred * (fLastReset - *f))
	lreset = TNC_TRUE;
      else
      {
	yrsr = tnc_ddot1(n, yr, sr);
	if (yrsr <= 0.0) lreset = TNC_TRUE;
	else lreset = TNC_FALSE;
      }
      upd1 = TNC_FALSE;
    }
  }
  /* Unscaling */
  tnc_unscalex(n, x, xscale, xoffset);
  tnc_coercex(n, x, low, up);
  (*f) /= *fscale;

  cac_progress(x, pk, *f, gnorm, alpha, niter, state);

cleanup:
  if (oldg) free(oldg);
  if (g) free(g);
  if (temp) free(temp);
  if (diagb) free(diagb);
  if (pk) free(pk);

  if (sk) free(sk);
  if (yk) free(yk);
  if (sr) free(sr);
  if (yr) free(yr);

  if (pivot) free(pivot);

  if (oldx) free(oldx);

  return rc;
}


/**
 *
 * Set x[i] = 0.0 if direction i is currently constrained
 *
 */

void tnc_project(int n, double x[], int pivot[])
{
  int i;
  for (i = 0; i < n; i++)
    if (pivot[i] != 0)
      x[i] = 0.0;
}

/**
 *
 * Set x[i] = 0.0 if direction i is constant
 *
 */

void tnc_projectConstants(int n, double x[], double xscale[])
{
  int i;
  for (i = 0; i < n; i++)
    if (xscale[i] == 0.0)
      x[i] = 0.0;
}

/**
 *
 * Compute the maximum allowable step length
 *
 */

double tnc_stepMax(double step, int n, double x[], double dir[],
    int pivot[], double low[], double up[], double xscale[], double xoffset[])
{
  int i;
  double t;

  /* Constrained maximum step */
  for (i = 0; i < n; i++)
  {
    if ((pivot[i] == 0) && (dir[i] != 0.0))
    {
      if (dir[i] < 0.0)
      {
	t = (low[i]-xoffset[i])/xscale[i] - x[i];
	if (t > step * dir[i]) step = t / dir[i];
      }
      else
      {
	t = (up[i]-xoffset[i])/xscale[i] - x[i];
	if (t < step * dir[i]) step = t / dir[i];
      }
    }
  }

  return step;
}

/**
 *
 * Update the constraint vector pivot if a new constraint is encountered
 *
 */

tnc_logical tnc_addConstraint(int n, double x[], double p[], int pivot[],
    double low[], double up[], double xscale[], double xoffset[])
{
  int i, newcon = TNC_FALSE;
  double tol, epsmch;

  epsmch = tnc_mchpr1();

  for (i = 0; i < n; i++)
  {
    if ((pivot[i] == 0) && (p[i] != 0.0))
    {
      if (p[i] < 0.0 && low[i] != - HUGE_VAL)
      {
	tol = epsmch * 10.0 * (fabs(low[i]) + 1.0);
	if (x[i]*xscale[i]+xoffset[i] - low[i] <= tol)
	{
	  pivot[i] = -1;
	  x[i] = (low[i]-xoffset[i])/xscale[i];
	  newcon = TNC_TRUE;
	}
      }
      else if (up[i] != HUGE_VAL)
      {
	tol = epsmch * 10.0 * (fabs(up[i]) + 1.0);
	if (up[i] - (x[i]*xscale[i]+xoffset[i]) <= tol)
	{
	  pivot[i] = 1;
	  x[i] = (up[i]-xoffset[i])/xscale[i];
	  newcon = TNC_TRUE;
	}
      }
    }
  }
  return newcon;
}

/**
 *
 * Check if a constraint is no more active
 *
 */

tnc_logical tnc_removeConstraint(double gtpnew, double gnorm, double pgtolfs,
    double f, double fLastConstraint, double g[], int pivot[], int n)
{
  double cmax, t;
  int imax, i;

  if (((fLastConstraint - f) <= (gtpnew * -0.5)) && (gnorm > pgtolfs))
    return TNC_FALSE;

  imax = -1;
  cmax = 0.0;

  for (i = 0; i < n; i++)
  {
    if (pivot[i] == 2)
      continue;
    t = -pivot[i] * g[i];
    if (t < cmax)
    {
      cmax = t;
      imax = i;
    }
  }

  if (imax != -1)
  {
    pivot[imax] = 0;
    return TNC_TRUE;
  }
  else
    return TNC_FALSE;

  /*
   * For details, see gill, murray, and wright (1981, p. 308) and
   * fletcher (1981, p. 116). The multiplier tests (here, testing
   * the sign of the components of the gradient) may still need to
   * modified to incorporate tolerances for zero.
   */
}

/**
 *
 * This routine performs a preconditioned conjugate-gradient
 * iteration in order to solve the newton equations for a search
 * direction for a truncated-newton algorithm.
 * When the value of the quadratic model is sufficiently reduced,
 * the iteration is terminated. Adapted by L. Garrido and M. Kalmoun Adapted by L. Garrido and M. Kalmoun.
 *
 */

int tnc_direction(double *zsol, double *diagb,
    double *x, double g[], int n,
    int maxCGit, int *niter, int maxnfeval, int *nfeval,
    tnc_logical upd1, double yksk, double yrsr,
    double *sk, double *yk, double *sr, double *yr,
    tnc_logical lreset, tnc_function *function, void *state,
    double xscale[], double xoffset[], double fscale,
    int *pivot, double accuracy,
    double gnorm, double xnorm, double low[], double up[])
{
  double alpha, beta, qold, qnew, rhsnrm, tol, vgv, rz, rzold, qtest, pr, gtp, gzsol, gzsoln;
  int i, k, frc;
  /* Temporary vectors */
  double *r = NULL, *zk = NULL, *v = NULL, *emat = NULL, *gv = NULL, *zsoln = NULL;

  *niter = 0;
  maxnfeval += 0; /* Avoid complaining of compiler */

  /* No CG it. => dir = -grad */
  if (maxCGit == 0)
  {
    tnc_dcopy1(n, g, zsol);
    tnc_dneg1(n, zsol);
    tnc_project(n, zsol, pivot);
    return 0;
  }

  /* General initialization */
  rhsnrm = gnorm;
  tol = 1e-12;
  qold = 0.0;
  rzold = 0.0; /* Uneeded */

  frc = -1; /* ENOMEM here */
  zsoln = malloc(sizeof(*zsol)*n); /* Residual */
  if (zsoln == NULL) goto cleanup;
  r = malloc(sizeof(*r)*n); /* Residual */
  if (r == NULL) goto cleanup;
  v = malloc(sizeof(*v)*n);
  if (v == NULL) goto cleanup;
  zk = malloc(sizeof(*zk)*n);
  if (zk == NULL) goto cleanup;
  emat = malloc(sizeof(*emat)*n); /* Diagonal preconditoning matrix */
  if (emat == NULL) goto cleanup;
  gv = malloc(sizeof(*gv)*n); /* hessian times v */
  if (gv == NULL) goto cleanup;

  /* Initialization for preconditioned conjugate-gradient algorithm */
  frc = tnc_initPreconditioner(diagb, emat, n, lreset, yksk, yrsr, sk, yk, sr, yr,
      upd1);
  if (frc) goto cleanup;

  for (i = 0; i < n; i++)
  {
    r[i]     = -g[i];
    v[i]     = 0.0;
    zsol[i]  = 0.0; /* Computed search direction */
    zsoln[i] = 0.0;
  }

  gzsol = 0.0;  

  /* Main iteration */
  for (k = 0; k < maxCGit; k++)
  {
    /* CG iteration to solve system of equations */
    tnc_project(n, r, pivot);
    frc = tnc_msolve(r, zk, n, sk, yk, diagb, sr, yr, upd1, yksk, yrsr, lreset);
    if (frc) goto cleanup;
    tnc_project(n, zk, pivot);
    rz = tnc_ddot1(n, r, zk);

    if (fabs(rz / rhsnrm) < tol) // Changed, L. Garrido
    {
      /* Truncate algorithm in case of an emergency
	 or too many function evaluations */
      if (k == 0)
      {
	tnc_dcopy1(n, g, zsol);
	tnc_dneg1(n, zsol);
	tnc_project(n, zsol, pivot);
      }
      break;
    }
    if (k == 0) beta = 0.0;
    else beta = rz / rzold;

    for (i = 0; i < n; i++)
      v[i] = zk[i] + beta * v[i];

    tnc_project(n, v, pivot);
    frc = tnc_hessianTimesVector(v, gv, n, x, g, function, state,
	xscale, xoffset, fscale, accuracy, xnorm, low, up);
    ++(*nfeval);
    if (frc) goto cleanup;
    tnc_project(n, gv, pivot);

    vgv = tnc_ddot1(n, v, gv);

    if (fabs(vgv / rhsnrm) < tol) // We use negative curvature afterwards.
    {
      /* Truncate algorithm in case of an emergency */
      if (k == 0)
      {
	tnc_dcopy1(n, zk, zsol);
	tnc_project(n, zsol, pivot);
      }
      break;
    }

    /* Compute linear step length */
    alpha = rz / vgv;

    /* Compute current solution and related vectors */
    tnc_daxpy1(n, alpha, v, zsoln);

    /* Check if we are in negative curvature */
    gzsoln = tnc_ddot1(n, g, zsoln);

    /* No descent direction. Replaces negative curvature. */
    if (gzsoln >= gzsol + tol) 
    {
      if (k == 0)
      {
	tnc_dcopy1(n, g, zsol);
	tnc_dneg1(n, zsol);
	tnc_project(n, zsol, pivot);
      }

      break;
    }

    gzsol = gzsoln;
    tnc_dcopy1(n, zsoln, zsol);

    /* Ok, let's continue */
    tnc_diagonalScaling(n, emat, v, gv, r);
    tnc_daxpy1(n, -alpha, gv, r);

    /* Test for convergence */
    gtp = tnc_ddot1(n, zsol, g);
    pr = tnc_ddot1(n, r, zsol);
    qnew = (gtp + pr) * 0.5;
    qtest = (k + 1) * (1.0 - qold / qnew);
    if (qtest <= 0.5) 
      break;

    /* Perform cautionary test */
    if (gtp > 0.0)
    {
      /* Truncate algorithm in case of an emergency */
      tnc_daxpy1(n, -alpha, v, zsol);
      break;
    }

    qold = qnew;
    rzold = rz;
  }

  *niter = k;

  /* Terminate algorithm */
  /* Store (or restore) diagonal preconditioning */
  tnc_dcopy1(n, emat, diagb);

cleanup:
  if (zsoln) free(zsoln);
  if (r) free(r);
  if (v) free(v);
  if (zk) free(zk);
  if (emat) free(emat);
  if (gv) free(gv);
  return frc;
}

/**
 *
 * Update the preconditioning matrix based on a diagonal version
 * of the bfgs quasi-newton update.
 *
 */

void tnc_diagonalScaling(int n, double e[], double v[], double gv[],
    double r[])
{
  int i;
  double vr, vgv;

  vr = 1.0/tnc_ddot1(n, v, r);
  vgv = 1.0/tnc_ddot1(n, v, gv);
  for (i = 0; i < n; i++)
  {
    e[i] += - r[i]*r[i]*vr + gv[i]*gv[i]*vgv;
    if (e[i] <= 1e-6) e[i] = 1.0;
  }
}

/**
 * 
 * Returns the length of the initial step to be taken along the
 * vector p in the next linear search.
 *
 */

double tnc_initialStep(double fnew, double fmin, double gtp, double smax)
{
  double d, alpha;

  d = fabs(fnew - fmin);
  alpha = 1.0;
  if (d * 2.0 <= -(gtp) && d >= tnc_mchpr1()) alpha = d * -2.0 / gtp;
  if (alpha >= smax) alpha = smax;

  return alpha;
}

/**
 *
 * Hessian vector product through finite differences
 *
 */

int tnc_hessianTimesVector(double v[], double gv[], int n,
    double x[], double g[], tnc_function *function, void *state,
    double xscale[], double xoffset[], double fscale,
    double accuracy, double xnorm, double low[], double up[])
{
  double dinv, f, delta, *xv;
  int i, frc;

  xv = malloc(sizeof(*xv)*n);
  if (xv == NULL) return -1;

  delta = accuracy * (xnorm + 1.0);
  for (i = 0; i < n; i++)
    xv[i] = x[i] + delta * v[i];

  tnc_unscalex(n, xv, xscale, xoffset);
  tnc_coercex(n, xv, low, up);
  // For the current implementation this function is NOT called
  cac_tnc_set_uniform_sampling(0, state); 
  frc = function(xv, &f, gv, state);
  free(xv);
  if (frc) return 1;
  tnc_scaleg(n, gv, xscale, fscale);

  dinv = 1.0 / delta;
  for (i = 0; i < n; i++)
    gv[i] = (gv[i] - g[i]) * dinv;

  tnc_projectConstants(n, gv, xscale);

  return 0;
}

/**
 *
 * This routine acts as a preconditioning step for the
 * linear conjugate-gradient routine. It is also the
 * method of computing the search direction from the
 * gradient for the non-linear conjugate-gradient code.
 * It represents a two-step self-scaled bfgs formula.
 *
 */

int tnc_msolve(double g[], double y[], int n,
    double sk[], double yk[], double diagb[], double sr[],
    double yr[], tnc_logical upd1, double yksk, double yrsr,
    tnc_logical lreset)
{
  double ghyk, ghyr, yksr, ykhyk, ykhyr, yrhyr, rdiagb, gsr, gsk;
  int i, frc;
  double *hg = NULL, *hyk = NULL, *hyr = NULL;

  if (upd1)
  {
    for (i = 0; i < n; i++) y[i] = g[i] / diagb[i];
    return 0;
  }

  frc = -1;
  gsk = tnc_ddot1(n, g, sk);
  hg = malloc(sizeof(*hg)*n);
  if (hg == NULL) goto cleanup;
  hyr = malloc(sizeof(*hyr)*n);
  if (hyr == NULL) goto cleanup;
  hyk = malloc(sizeof(*hyk)*n);
  if (hyk == NULL) goto cleanup;
  frc = 0;

  /* Compute gh and hy where h is the inverse of the diagonals */
  if (lreset)
  {
    for (i = 0; i < n; i++)
    {
      rdiagb = 1.0 / diagb[i];
      hg[i] = g[i] * rdiagb;
      hyk[i] = yk[i] * rdiagb;
    }
    ykhyk = tnc_ddot1(n, yk, hyk);
    ghyk = tnc_ddot1(n, g, hyk);
    tnc_ssbfgs(n, 1.0, sk, hg, hyk, yksk, ykhyk, gsk, ghyk, y);
  }
  else
  {
    for (i = 0; i < n; i++)
    {
      rdiagb = 1.0 / diagb[i];
      hg[i] = g[i] * rdiagb;
      hyk[i] = yk[i] * rdiagb;
      hyr[i] = yr[i] * rdiagb;
    }
    gsr = tnc_ddot1(n, g, sr);
    ghyr = tnc_ddot1(n, g, hyr);
    yrhyr = tnc_ddot1(n, yr, hyr);
    tnc_ssbfgs(n, 1.0, sr, hg, hyr, yrsr, yrhyr, gsr, ghyr, hg);
    yksr = tnc_ddot1(n, yk, sr);
    ykhyr = tnc_ddot1(n, yk, hyr);
    tnc_ssbfgs(n, 1.0, sr, hyk, hyr, yrsr, yrhyr, yksr, ykhyr, hyk);
    ykhyk = tnc_ddot1(n, hyk, yk);
    ghyk = tnc_ddot1(n, hyk, g);
    tnc_ssbfgs(n, 1.0, sk, hg, hyk, yksk, ykhyk, gsk, ghyk, y);
  }

cleanup:
  if (hg) free(hg);
  if (hyk) free(hyk);
  if (hyr) free(hyr);

  return frc;
}

/**
 *
 * Self-scaled BFGS
 *
 */

void tnc_ssbfgs(int n, double gamma, double sj[], double hjv[],
    double hjyj[], double yjsj,
    double yjhyj, double vsj, double vhyj, double hjp1v[])
{
  double beta, delta;
  int i;

  if (yjsj == 0.0)
  {
    delta = 0.0;
    beta = 0.0;
  }
  else
  {
    delta = (gamma * yjhyj / yjsj + 1.0) * vsj / yjsj - gamma * vhyj / yjsj;
    beta = -gamma * vsj / yjsj;
  }

  for (i = 0; i < n; i++)
    hjp1v[i] = gamma * hjv[i] + delta * sj[i] + beta * hjyj[i];
}

/**
 *
 * Initialize the preconditioner
 *
 */

int tnc_initPreconditioner(double diagb[], double emat[], int n,
    tnc_logical lreset, double yksk, double yrsr,
    double sk[], double yk[], double sr[], double yr[],
    tnc_logical upd1)
{
  double srds, yrsk, td, sds;
  int i;
  double *bsk;

  if (upd1)
  {
    tnc_dcopy1(n, diagb, emat);
    return 0;
  }

  bsk = malloc(sizeof(*bsk)*n);
  if (bsk == NULL) return -1;

  if (lreset)
  {
    for (i = 0; i < n; i++) bsk[i] = diagb[i] * sk[i];
    sds = tnc_ddot1(n, sk, bsk);
    if (yksk == 0.0) yksk = 1.0;
    if (sds == 0.0) sds = 1.0;
    for (i = 0; i < n; i++)
    {
      td = diagb[i];
      emat[i] = td - td * td * sk[i] * sk[i] / sds + yk[i] * yk[i] / yksk;
    }
  }
  else
  {
    for (i = 0; i < n; i++) bsk[i] = diagb[i] * sr[i];
    sds = tnc_ddot1(n, sr, bsk);
    srds = tnc_ddot1(n, sk, bsk);
    yrsk = tnc_ddot1(n, yr, sk);
    if (yrsr == 0.0) yrsr = 1.0;
    if (sds == 0.0) sds = 1.0;
    for (i = 0; i < n; i++)
    {
      td = diagb[i];
      bsk[i] = td * sk[i] - bsk[i] * srds / sds + yr[i] * yrsk / yrsr;
      emat[i] = td - td * td * sr[i] * sr[i] / sds + yr[i] * yr[i] / yrsr;
    }
    sds = tnc_ddot1(n, sk, bsk);
    if (yksk == 0.0) yksk = 1.0;
    if (sds == 0.0) sds = 1.0;
    for (i = 0; i < n; i++)
      emat[i] = emat[i] - bsk[i] * bsk[i] / sds + yk[i] * yk[i] / yksk;
  }

  free(bsk);
  return 0;
}

/**
 * 
 * Line search algorithm. Taken from L-BFGS code, see http://www.chokkan.org/software/liblbfgs/.
 * We have taken the line search implementation of the More-Thuente method and adapted it to the
 * interface of the TNC method.
 *
 */

#define min2(a, b)      ((a) <= (b) ? (a) : (b))
#define max2(a, b)      ((a) >= (b) ? (a) : (b))
#define max3(a, b, c)   max2(max2((a), (b)), (c));

tnc_ls_rc tnc_linearSearch(int n, double *low, double *up, double *xscale, double *xoffset, double fscale, int *pivot,
    double *x, double *f, double *gfull, double *s, double *stp, double max_step, double dginit, double xtol, int maxnfeval, 
    int *count, tnc_function *function, void *state, int only_line_search)
{
  int frc, rc, niter;
  int brackt, stage1, uinfo = 0;
  double dg;
  double stx, fx, dgx;
  double sty, fy, dgy;
  double fxm, dgxm, fym, dgym, fm, dgm;
  double finit, ftest1, dgtest;
  double width, prev_width;
  double min_step, stmin, stmax;

  double ftol, gtol;
  double *g, *wa, *tempx;

  pivot[0] = 0; /* Not used in this function. Used to avoid complaining of compiler */

  g     = (double *) malloc(sizeof(double) * n);
  wa    = (double *) malloc(sizeof(double) * n);
  tempx = (double *) malloc(sizeof(double) * n);

  if ((!g) || (!wa) || (!tempx))
  {
    rc = LS_ENOMEM;
    goto cleanup;
  }

  /**
   * The minimum step of the line search routine.
   *  The default value is \c 1e-20. This value need not be modified unless
   *  the exponents are too large for the machine being used, or unless the
   *  problem is extremely badly scaled (in which case the exponents should
   *  be increased).
   */

  min_step = 1e-20;

  /**
   * The maximum step of the line search.
   *  The default value is \c 1e+20. This value need not be modified unless
   *  the exponents are too large for the machine being used, or unless the
   *  problem is extremely badly scaled (in which case the exponents should
   *  be increased).
   */

  //max_step is automatically set by calling function
  //max_step = 1e+20;

  /**
   * A parameter to control the accuracy of the line search routine.
   *  The default value is \c 1e-4. This parameter should be greater
   *  than zero and smaller than \c 0.5.
   */

  ftol = 1e-04;

  /**
   * A parameter to control the accuracy of the line search routine.
   *  The default value is \c 0.9. If the function and gradient
   *  evaluations are inexpensive with respect to the cost of the
   *  iteration (which is sometimes the case when solving very large
   *  problems) it may be advantageous to set this parameter to a small
   *  value. A typical small value is \c 0.1. This parameter shuold be
   *  greater than the \ref ftol parameter (\c 1e-4) and smaller than
   *  \c 1.0.
   */

  gtol = 0.9;

  /* Just to avoid complaining of the compiler */

  rc = LS_FAIL;

  /* Check the input parameters for errors. */
  if (*stp <= 0.) {
    rc = LS_FAIL;
    goto cleanup;
  }

  //not needed, is computed outside
  //dginit = tnc_ddot1(n, g, s);

  /* Make sure that s points to a descent direction. */
  if (0 < dginit) {
    rc = LS_FAIL;
    goto cleanup;
  }

  /* Initialize local variables. */
  brackt = 0;
  stage1 = 1;
  finit = *f;
  dgtest = ftol * dginit;
  width = max_step - min_step;
  prev_width = 2.0 * width;

  /* Copy the value of x to the work area. */
  tnc_dcopy1(n, x, wa);

  /*
     The variables stx, fx, dgx contain the values of the step,
     function, and directional derivative at the best step.
     The variables sty, fy, dgy contain the value of the step,
     function, and derivative at the other endpoint of
     the interval of uncertainty.
     The variables stp, f, dg contain the values of the step,
     function, and derivative at the current step.
   */
  stx = sty = 0.;
  fx = fy = finit;
  dgx = dgy = dginit;

  for (niter = 0;; niter++) {
    /*
       Set the minimum and maximum steps to correspond to the
       present interval of uncertainty.
     */
    if (brackt) {
      stmin = min2(stx, sty);
      stmax = max2(stx, sty);
    } else {
      stmin = stx;
      stmax = *stp + 4.0 * (*stp - stx);
    }

    /* Clip the step in the range of [stpmin, stpmax]. */
    if (*stp < min_step) *stp = min_step;
    if (max_step < *stp) *stp = max_step;

    /*
       If an unusual termination is to occur then let
       stp be the lowest point obtained so far.
     */
    if ((brackt && ((*stp <= stmin || stmax <= *stp) || maxnfeval <= *count + 1 || uinfo != 0)) || 
	(brackt && (stmax - stmin <= xtol * stmax))) {
      *stp = stx;
    }

    /*
       Compute the current value of x:
       x <- x + (*stp) * s.
     */

    if (isnan(*stp))
    {
      fprintf(stderr, "\nERROR: stp == nan in line search\n");
      exit(1);
    }

    tnc_dcopy1(n, wa, x);
    tnc_daxpy1(n, *stp, s, x);

    /* Scale and constraint current x value */

    tnc_dcopy1(n, x, tempx);
    tnc_unscalex(n, tempx, xscale, xoffset);
    tnc_coercex(n, tempx, low, up);

    /* Once the values have been scaled, evaluate function updating lists */

    cac_tnc_set_uniform_sampling(1, state);
    cac_tnc_setup_omega1_omega2(tempx, state);
    frc = function(tempx, f, gfull, state);
    if (frc)
    {
      rc = LS_USERABORT;
      goto cleanup;
    }
    ++(*count);

    /* Scale function and gradient */

    *f *= fscale;

    tnc_dcopy1(n, gfull, g);
    tnc_scaleg(n, g, xscale, fscale);

    dg = tnc_ddot1(n, g, s);
    ftest1 = finit + *stp * dgtest;

    /* Test for errors and convergence. */
    if (*stp < THRSALPHA) {
      rc = LS_OK;
      goto cleanup;
    }
    if (brackt && ((*stp <= stmin || stmax <= *stp) || uinfo != 0)) {
      /* Rounding errors prevent further progress. */
      rc = LS_OK;
      goto cleanup;
    }
    if (*stp == max_step && *f <= ftest1 && dg <= dgtest) {
      /* The step is the maximum value. */
      rc = LS_MAXSTP; 
      goto cleanup;
    }
    if (*stp == min_step && (ftest1 < *f || dgtest <= dg)) {
      /* The step is the minimum value. We'll return a LS_OK since
       * the there may be some inaccuracies in the direction of descent. */
      rc = LS_OK;
      goto cleanup;
    }
    if (brackt && (stmax - stmin) <= xtol * stmax) {
      /* Relative width of the interval of uncertainty is at most xtol. */
      rc = LS_FAIL;
      goto cleanup;
    }
    if (maxnfeval <= *count) {
      /* Maximum number of iteration. */
      rc = LS_MAXFUN;
      goto cleanup;
    }
    if (*f <= ftest1 && fabs(dg) <= gtol * (-dginit)) {
      /* The sufficient decrease condition and the directional derivative condition hold. */
      rc = LS_OK;
      goto cleanup; // count contains the number of function evaluations
    }

    /*
       In the first stage we seek a step for which the modified
       function has a nonpositive value and nonnegative derivative.
     */
    if (stage1 && *f <= ftest1 && min2(ftol, gtol) * dginit <= dg) {
      stage1 = 0;
    }

    /*
       A modified function is used to predict the step only if
       we have not obtained a step for which the modified
       function has a nonpositive function value and nonnegative
       derivative, and if a lower function value has been
       obtained but the decrease is not sufficient.
     */
    if (stage1 && ftest1 < *f && *f <= fx) {
      /* Define the modified function and derivative values. */
      fm = *f - *stp * dgtest;
      fxm = fx - stx * dgtest;
      fym = fy - sty * dgtest;
      dgm = dg - dgtest;
      dgxm = dgx - dgtest;
      dgym = dgy - dgtest;

      /*
	 Call update_trial_interval() to update the interval of
	 uncertainty and to compute the new step.
       */
      uinfo = tnc_update_trial_interval(
	  &stx, &fxm, &dgxm,
	  &sty, &fym, &dgym,
	  stp, &fm, &dgm,
	  stmin, stmax, &brackt
	  );

      /* Reset the function and gradient values for f. */
      fx = fxm + stx * dgtest;
      fy = fym + sty * dgtest;
      dgx = dgxm + dgtest;
      dgy = dgym + dgtest;
    } else {
      /*
	 Call update_trial_interval() to update the interval of
	 uncertainty and to compute the new step.
       */
      uinfo = tnc_update_trial_interval(
	  &stx, &fx, &dgx,
	  &sty, &fy, &dgy,
	  stp, f, &dg,
	  stmin, stmax, &brackt
	  );
    }

    /*
       Force a sufficient decrease in the interval of uncertainty.
     */
    if (brackt) {
      if (0.66 * prev_width <= fabs(sty - stx)) {
	*stp = stx + 0.5 * (sty - stx);
      }
      prev_width = width;
      width = fabs(sty - stx);
    }
  }

cleanup:

  if (g)     free((char *) g);
  if (wa)    free((char *) wa);
  if (tempx) free((char *) tempx);

  return rc;
}

/**
 * Define the local variables for computing minimizers.
 */
#define USES_MINIMIZER \
  double a, d, gamma, theta, p, q, r, s;

/**
 * Find a minimizer of an interpolated cubic function.
 *  @param  cm      The minimizer of the interpolated cubic.
 *  @param  u       The value of one point, u.
 *  @param  fu      The value of f(u).
 *  @param  du      The value of f'(u).
 *  @param  v       The value of another point, v.
 *  @param  fv      The value of f(v).
 *  @param  du      The value of f'(v).
 */
#define CUBIC_MINIMIZER(cm, u, fu, du, v, fv, dv) \
  d = (v) - (u); \
theta = ((fu) - (fv)) * 3 / d + (du) + (dv); \
p = fabs(theta); \
q = fabs(du); \
r = fabs(dv); \
s = max3(p, q, r); \
/* gamma = s*sqrt((theta/s)**2 - (du/s) * (dv/s)) */ \
a = theta / s; \
gamma = s * sqrt(a * a - ((du) / s) * ((dv) / s)); \
if ((v) < (u)) gamma = -gamma; \
p = gamma - (du) + theta; \
q = gamma - (du) + gamma + (dv); \
r = p / q; \
(cm) = (u) + r * d;

/**
 * Find a minimizer of an interpolated cubic function.
 *  @param  cm      The minimizer of the interpolated cubic.
 *  @param  u       The value of one point, u.
 *  @param  fu      The value of f(u).
 *  @param  du      The value of f'(u).
 *  @param  v       The value of another point, v.
 *  @param  fv      The value of f(v).
 *  @param  du      The value of f'(v).
 *  @param  xmin    The maximum value.
 *  @param  xmin    The minimum value.
 */
#define CUBIC_MINIMIZER2(cm, u, fu, du, v, fv, dv, xmin, xmax) \
  d = (v) - (u); \
theta = ((fu) - (fv)) * 3 / d + (du) + (dv); \
p = fabs(theta); \
q = fabs(du); \
r = fabs(dv); \
s = max3(p, q, r); \
/* gamma = s*sqrt((theta/s)**2 - (du/s) * (dv/s)) */ \
a = theta / s; \
gamma = s * sqrt(max2(0, a * a - ((du) / s) * ((dv) / s))); \
if ((u) < (v)) gamma = -gamma; \
p = gamma - (dv) + theta; \
q = gamma - (dv) + gamma + (du); \
r = p / q; \
if (r < 0. && gamma != 0.) { \
  (cm) = (v) - r * d; \
} else if (a < 0) { \
  (cm) = (xmax); \
} else { \
  (cm) = (xmin); \
}

/**
 * Find a minimizer of an interpolated quadratic function.
 *  @param  qm      The minimizer of the interpolated quadratic.
 *  @param  u       The value of one point, u.
 *  @param  fu      The value of f(u).
 *  @param  du      The value of f'(u).
 *  @param  v       The value of another point, v.
 *  @param  fv      The value of f(v).
 */
#define QUARD_MINIMIZER(qm, u, fu, du, v, fv) \
  a = (v) - (u); \
(qm) = (u) + (du) / (((fu) - (fv)) / a + (du)) / 2 * a;

/**
 * Find a minimizer of an interpolated quadratic function.
 *  @param  qm      The minimizer of the interpolated quadratic.
 *  @param  u       The value of one point, u.
 *  @param  du      The value of f'(u).
 *  @param  v       The value of another point, v.
 *  @param  dv      The value of f'(v).
 */
#define QUARD_MINIMIZER2(qm, u, du, v, dv) \
  a = (u) - (v); \
(qm) = (v) + (dv) / ((dv) - (du)) * a;

/**
 * Update a safeguarded trial value and interval for line search.
 *
 *  The parameter x represents the step with the least function value.
 *  The parameter t represents the current step. This function assumes
 *  that the derivative at the point of x in the direction of the step.
 *  If the bracket is set to true, the minimizer has been bracketed in
 *  an interval of uncertainty with endpoints between x and y.
 *
 *  @param  x       The pointer to the value of one endpoint.
 *  @param  fx      The pointer to the value of f(x).
 *  @param  dx      The pointer to the value of f'(x).
 *  @param  y       The pointer to the value of another endpoint.
 *  @param  fy      The pointer to the value of f(y).
 *  @param  dy      The pointer to the value of f'(y).
 *  @param  t       The pointer to the value of the trial value, t.
 *  @param  ft      The pointer to the value of f(t).
 *  @param  dt      The pointer to the value of f'(t).
 *  @param  tmin    The minimum value for the trial value, t.
 *  @param  tmax    The maximum value for the trial value, t.
 *  @param  brackt  The pointer to the predicate if the trial value is
 *                  bracketed.
 *  @retval int     Status value. Zero indicates a normal termination.
 *  
 *  @see
 *      Jorge J. More and David J. Thuente. Line search algorithm with
 *      guaranteed sufficient decrease. ACM Transactions on Mathematical
 *      Software (TOMS), Vol 20, No 3, pp. 286-307, 1994.
 */

#define fsigndiff(x, y) (*(x) * (*(y) / fabs(*(y))) < 0.)

int tnc_update_trial_interval(
    double *x,
    double *fx,
    double *dx,
    double *y,
    double *fy,
    double *dy,
    double *t,
    double *ft,
    double *dt,
    const double tmin,
    const double tmax,
    int *brackt
    )
{
  int bound;
  int dsign = fsigndiff(dt, dx);
  double mc; /* minimizer of an interpolated cubic. */
  double mq; /* minimizer of an interpolated quadratic. */
  double newt;   /* new trial value. */
  USES_MINIMIZER;     /* for CUBIC_MINIMIZER and QUARD_MINIMIZER. */

  /* Check the input parameters for errors. */
  if (*brackt) {
    if (*t <= min2(*x, *y) || max2(*x, *y) <= *t) {
      /* The trival value t is out of the interval. */
      return LS_FAIL;
    }
    if (0. <= *dx * (*t - *x)) {
      /* The function must decrease from x. */
      return LS_FAIL;
    }
    if (tmax < tmin) {
      /* Incorrect tmin and tmax specified. */
      return LS_FAIL;
    }
  }

  /*
     Trial value selection.
   */
  if (*fx < *ft) {
    /*
       Case 1: a higher function value.
       The minimum is brackt. If the cubic minimizer is closer
       to x than the quadratic one, the cubic one is taken, else
       the average of the minimizers is taken.
     */
    *brackt = 1;
    bound = 1;
    CUBIC_MINIMIZER(mc, *x, *fx, *dx, *t, *ft, *dt);
    QUARD_MINIMIZER(mq, *x, *fx, *dx, *t, *ft);
    if (fabs(mc - *x) < fabs(mq - *x)) {
      newt = mc;
    } else {
      newt = mc + 0.5 * (mq - mc);
    }
  } else if (dsign) {
    /*
       Case 2: a lower function value and derivatives of
       opposite sign. The minimum is brackt. If the cubic
       minimizer is closer to x than the quadratic (secant) one,
       the cubic one is taken, else the quadratic one is taken.
     */
    *brackt = 1;
    bound = 0;
    CUBIC_MINIMIZER(mc, *x, *fx, *dx, *t, *ft, *dt);
    QUARD_MINIMIZER2(mq, *x, *dx, *t, *dt);
    if (fabs(mc - *t) > fabs(mq - *t)) {
      newt = mc;
    } else {
      newt = mq;
    }
  } else if (fabs(*dt) < fabs(*dx)) {
    /*
       Case 3: a lower function value, derivatives of the
       same sign, and the magnitude of the derivative decreases.
       The cubic minimizer is only used if the cubic tends to
       infinity in the direction of the minimizer or if the minimum
       of the cubic is beyond t. Otherwise the cubic minimizer is
       defined to be either tmin or tmax. The quadratic (secant)
       minimizer is also computed and if the minimum is brackt
       then the the minimizer closest to x is taken, else the one
       farthest away is taken.
     */
    bound = 1;
    CUBIC_MINIMIZER2(mc, *x, *fx, *dx, *t, *ft, *dt, tmin, tmax);
    QUARD_MINIMIZER2(mq, *x, *dx, *t, *dt);
    if (*brackt) {
      if (fabs(*t - mc) < fabs(*t - mq)) {
	newt = mc;
      } else {
	newt = mq;
      }
    } else {
      if (fabs(*t - mc) > fabs(*t - mq)) {
	newt = (mc!=tmin)?mc:mq; // Modification by L. Garrido
      } else {
	newt = (mq!=tmin)?mq:mc; // Modification by L. Garrido
      }

      /* NOTA: Se ha modificado el codigo anterior, puesto que parece que en
       * algun caso particular no escoge bien cual de los dos valores debe
       * escoger. Tal y como estaba programado antes, escogia newt = mc, donde
       * mc = tmin y eso hacia que el algoritmo funcionase incorrectamente en
       * las iteraciones siguientes. */
    }
  } else {
    /*
       Case 4: a lower function value, derivatives of the
       same sign, and the magnitude of the derivative does
       not decrease. If the minimum is not brackt, the step
       is either tmin or tmax, else the cubic minimizer is taken.
     */
    bound = 0;
    if (*brackt) {
      CUBIC_MINIMIZER(newt, *t, *ft, *dt, *y, *fy, *dy);
    } else if (*x < *t) {
      newt = tmax;
    } else {
      newt = tmin;
    }
  }

  /*
     Update the interval of uncertainty. This update does not
     depend on the new step or the case analysis above.

     - Case a: if f(x) < f(t),
     x <- x, y <- t.
     - Case b: if f(t) <= f(x) && f'(t)*f'(x) > 0,
     x <- t, y <- y.
     - Case c: if f(t) <= f(x) && f'(t)*f'(x) < 0, 
     x <- t, y <- x.
   */
  if (*fx < *ft) {
    /* Case a */
    *y = *t;
    *fy = *ft;
    *dy = *dt;
  } else {
    /* Case c */
    if (dsign) {
      *y = *x;
      *fy = *fx;
      *dy = *dx;
    }
    /* Cases b and c */
    *x = *t;
    *fx = *ft;
    *dx = *dt;
  }

  /* Clip the new trial value in [tmin, tmax]. */
  if (tmax < newt) newt = tmax;
  if (newt < tmin) newt = tmin;

  /*
     Redefine the new trial value if it is close to the upper bound
     of the interval.
   */
  if (*brackt && bound) {
    mq = *x + 0.66 * (*y - *x);
    if (*x < *y) {
      if (mq < newt) newt = mq;
    } else {
      if (newt < mq) newt = mq;
    }
  }

  /* Return the new trial value. */
  *t = newt;
  return 0;
}

/**
 *
 * Return epsmch, where epsmch is the smallest possible
 * power of 2 such that 1.0 + epsmch > 1.0
 *
 */
double tnc_mchpr1(void)
{
  static double epsmch = 0.0;

  if (epsmch == 0.0)
  {
    double eps = 1.0;
    while((1.0 + (eps*0.5)) > 1.0)
      eps *= 0.5;
    epsmch = eps;
  }

  return epsmch;
}

/**
 *  
 *  Blas like routines 
 *
 */

/* dy+=dx */
void tnc_dxpy1(int n, double dx[], double dy[])
{
  int i;
  for (i = 0; i < n; i++)
    dy[i] += dx[i];
}

/* dy+=da*dx */
void tnc_daxpy1(int n, double da, double dx[], double dy[])
{
  int i;
  for (i = 0; i < n; i++)
    dy[i] += da*dx[i];
}

/* Copy dx -> dy */
/* Could use memcpy */
void tnc_dcopy1(int n, double dx[], double dy[])
{
  int i;
  for (i = 0; i < n; i++)
    dy[i] = dx[i];
}

/* Negate */
void tnc_dneg1(int n, double v[])
{
  int i;
  for (i = 0; i < n; i++)
    v[i] = -v[i];
}

/* Dot product */
double tnc_ddot1(int n, double dx[], double dy[])
{
  int i;
  double dtemp = 0.0;
  for (i = 0; i < n; i++)
    dtemp += dy[i]*dx[i];
  return dtemp;
}

/* Euclidian norm */
double tnc_dnrm21(int n, double dx[])
{
  int i;
  double dssq = 1.0, dscale = 0.0;

  for (i = 0; i < n; i++)
  {
    if (dx[i] != 0.0)
    {
      double dabsxi = fabs(dx[i]);
      if (dscale<dabsxi)
      {
	/* Normalization to prevent overflow */
	double ratio = dscale/dabsxi;
	dssq = 1.0 + dssq*ratio*ratio;
	dscale = dabsxi;
      }
      else
      {
	double ratio = dabsxi/dscale;
	dssq += ratio*ratio;
      }
    }
  }

  return dscale*sqrt(dssq/n);
}
