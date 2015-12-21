// This program is free software: you can use, modify and/or redistribute it
// under the terms of the simplified BSD License. You should have received a
// copy of this license along this program. If not, see
// <http://www.opensource.org/licenses/bsd-license.html>.
//
// Copyright (C) 2015, Lluis Garrido <lluis.garrido@ub.edu>
// All rights reserved.


#include "cac.h"
#include "tnc.h"

/**
 *
 *  Applies Newton algorithm using TNC library
 *
 */

void cac_newton_algorithm(
    int n,
    double *x,
    tnc_function *proc_evaluate,
    struct resources_common *res_common)
{
  int i, messages, rc, nfeval, ok;
  double f, *low, *up, *g;

  int maxiter   = MAXITER;
  int maxnfeval = 2000000;

  double fmin = 0.0;

  double eta      = -1.0;
  double stepmx   = -1.0;
  double accuracy = -1.0;
  double rescale  = -1.0;

  double ftol  = 1e-05; 
  double xtol  = 1e-05;
  double pgtol = 1e-05;

  // Gradient descent
  int maxCGit  = 0;

  g = cac_xmalloc(sizeof(double) * n);

  low = cac_xmalloc(sizeof(double) * n);
  up  = cac_xmalloc(sizeof(double) * n);

  for(i = 0; i < n; i++)
  {
    low[i] = - HUGE_VAL;
    up[i]  = HUGE_VAL;
  }

  messages = TNC_MSG_ALL;
  cac_print_msg(res_common->fp_log, "------------------------------------------------\n");
  cac_print_msg(res_common->fp_log, "Maximum number of iterations: %d\n", MAXITER);
  cac_print_msg(res_common->fp_log, "Maximum vertex displacement: %d\n", MAXMOVE);
  cac_print_msg(res_common->fp_log, "Default vertex displacement: %d\n", DEFAULTMOVE);
  cac_print_msg(res_common->fp_log, "Discrete sampling for constraints: %d\n", SAMPLECONSTRAINTS);
  cac_print_msg(res_common->fp_log, "Band size: %d\n", BANDSIZE);
  cac_print_msg(res_common->fp_log, "Allow crosssings: %d\n", ALLOWCROSSINGS);
  cac_print_msg(res_common->fp_log, "Alpha convergence criterion: %f\n", THRSALPHA);
  cac_print_msg(res_common->fp_log, "Tolerance ftol, xtol: %e, %e\n", ftol, xtol);
  cac_print_msg(res_common->fp_log, "------------------------------------------------\n");

  rc = tnc(n, x, &f, g, proc_evaluate,  (void *) res_common, low, up, NULL, NULL, messages,
      maxiter, maxCGit, maxnfeval, eta, stepmx, accuracy, fmin, ftol, xtol, pgtol,
      rescale, &nfeval);

  cac_print_msg(res_common->fp_log, "Iterations finished\n");

  ok = 0;

  switch (rc) {
    case TNC_FCONVERGED:
    case TNC_XCONVERGED:
    case TNC_LOCALMINIMUM:
      ok = 1;
      break;

    default:
      ok = 0;
  }

  if (!ok) {
    printf("TNC code: %d\n", rc);
    printf("ERROR: tnc not converged\n");
  }

  free(g);
  free(low);
  free(up);
}


