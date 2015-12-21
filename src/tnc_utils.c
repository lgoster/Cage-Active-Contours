// This program is free software: you can use, modify and/or redistribute it
// under the terms of the simplified BSD License. You should have received a
// copy of this license along this program. If not, see
// <http://www.opensource.org/licenses/bsd-license.html>.
//
// Copyright (C) 2015, Lluis Garrido <lluis.garrido@ub.edu>
// All rights reserved.

#include "cac.h"

/**
 * 
 *  These functions are interfaces from the TNC library to
 *  our functions.
 *
 */

/**
 *
 *  Print a message
 *
 */

void cac_tnc_print_msg(char *msg, void *instance)
{

  struct resources_common *res_common = ((struct resources_common *) instance);

  cac_print_msg(res_common->fp_log, "%s", msg);
}

/**
 *
 *  Setup omega1 (interior) and omega2 (exterior) regions
 *
 */

void cac_tnc_setup_omega1_omega2(double *xin, void *instance)
{
  struct resources_common *res_common = ((struct resources_common *) instance);

  cac_setup_omega1_omega2(xin, 0, res_common);
}

/**
 *
 *  Allows to setup uniform or non-uniform sampling. This function is
 *  indeed currently not used. It may be used in the future when the
 *  Newton direction (rather than the gradient descent direction) is
 *  used to minimize the energy. If you suceed, tell me!
 *
 */

void cac_tnc_set_uniform_sampling(int value, void *instance)
{
  struct resources_common *res_common = ((struct resources_common *) instance);

  res_common->uniform_sampling = value;
}

/**
 *
 *  Allows to enable or disable the projection of the gradient
 *
 */

void cac_tnc_set_project_gradient(int value, void *instance)
{
  struct resources_common *res_common = ((struct resources_common *) instance);

  res_common->project_gradient = value;
}

/**
 *
 *  Gets the value of project_gradient
 *
 */

int cac_tnc_get_project_gradient(void *instance)
{
  int rc;
  struct resources_common *res_common = ((struct resources_common *) instance);

  rc = res_common->project_gradient;

  return rc;
}

void cac_tnc_project_gradient_line_vertex_to_center(double *xin, double *g, void *instance)
{
  struct resources_common *res_common = ((struct resources_common *) instance);
  int i;

  double cx = res_common->cx;
  double cy = res_common->cy;
  double dx, dy, norm, prod;

  for(i = 0; i < res_common->cage->size; i++)
  {
    dx = cx - xin[2*i];
    dy = cy - xin[2*i+1]; 
    norm = sqrt(dx * dx + dy * dy);

    dx /= norm;
    dy /= norm;

    prod = g[2*i] * dx + g[2*i+1] * dy;

    g[2*i] = dx * prod;
    g[2*i+1] = dy * prod;
  }
}


/*
 *
 * Compute max step for the linear search
 *
 */

double cac_tnc_get_maxstep_linsrch(int n, double *x, double *pk, void *instance)
{
  int i, count, maxcount, n2;
  double normstep, cost, maxstep, norm, maxnorm;
  double *tmp_pk, *tmp_x;

  struct resources_common *res_common = ((struct resources_common *) instance);

  tmp_pk = (double *) cac_xmalloc(sizeof(double) * n);
  tmp_x  = (double *) cac_xmalloc(sizeof(double) * n);

  /* Normalize gradient so that maximum norm corresponds to 1.0 pixels */

  n2 = n / 2;

  maxnorm = 0.0;
  for(i = 0; i < n2; i++)
  {
    norm = sqrt(SQR(pk[2*i]) + SQR(pk[2*i+1]));
    maxnorm = MAX(maxnorm, norm);
  }

  for(i = 0; i < n2; i++)
  {
    tmp_pk[2*i] = pk[2*i] / maxnorm;
    tmp_pk[2*i+1] = pk[2*i+1] / maxnorm;
  }

  /* Compute maximum step */

  maxcount = (double) MAXMOVE / (double) SAMPLECONSTRAINTS;

  cost = 0.0;
  count = 0;

  while ((cost == 0.0) && (count <= maxcount))
  {
    normstep = (double) count * (double) SAMPLECONSTRAINTS;

    for(i = 0; i < n2; i++) 
    {
      tmp_x[2*i] = x[2*i] + normstep * tmp_pk[2*i];
      tmp_x[2*i+1] = x[2*i+1] + normstep * tmp_pk[2*i+1];
    }

    cost = cac_compute_dist_constraints(res_common, tmp_x);
    count++;
  }

  count--;

  normstep = (double) count * (double) SAMPLECONSTRAINTS;
  maxstep = normstep / maxnorm;

  free(tmp_pk);
  free(tmp_x);

  return maxstep;  
}

/**
 *
 * Normalitze gradient so that the maximum component of the
 * gradient has norm 1.0.
 *
 */

void cac_tnc_normalize_gradient(int n, double *g)
{
  int i;
  int n2 = n / 2;
  double value, maxvalue;

  double *g_norm = malloc(sizeof(double) * n);

  if (!g_norm)
  {
    printf("ERROR: could not allocate memory for g_norm\n");
    exit(1);
  }

  maxvalue = 0.0;
  for(i = 0; i < n2; i++)
  {
    value = sqrt(SQR(g[2*i]) + SQR(g[2*i+1]));
    maxvalue = MAX(value, maxvalue);
  }

  for(i = 0; i < n2; i++)
  {
    g_norm[2*i] = g[2*i] / maxvalue;
    g_norm[2*i+1] = g[2*i+1] / maxvalue;
  }

  value = 0;
  for(i = 0; i < n; i++)
    value += g[i] * g_norm[i];

  if (value > 0.0)
  {
    for(i = 0; i < n; i++)
      g[i] = g_norm[i];
  }
  else
    printf("ERROR: normalize_gradient, check_gnorm\n");

  free(g_norm);
}

/**
 *
 *  Print our progress to screen and save some data to disc
 *
 */


int cac_progress(double *x, double *pk, double fx, double gnorm, double alpha, int k, void *instance)
{
  int j;

  struct resources_common *res_common = ((struct resources_common *) instance);

  /* Log data */

  FILE *fp;

  char sz_string[1024], out_name[256];

  sprintf(sz_string, "Iteration %3d : fx = %E, gnorm = %E  alpha = %E\n", k, fx, gnorm, alpha);
  cac_print_msg(res_common->fp_log, sz_string);

  sprintf(out_name,"%s/contour_%03d.txt",LOGDIR, k);
  cac_write_points2d(res_common->contour, 1, out_name);

  sprintf(out_name,"%s/cage_%03d.txt",LOGDIR, k);

  fp = fopen(out_name, "w");
  for(int j = 0; j < res_common->cage->size; j++)
    fprintf(fp, "%E %E\n", x[2*j], x[2*j+1]);

  fprintf(fp, "%E %E\n", x[0], x[1]);
  fclose(fp);

  /* Copy current cage to resources */

  for(j = 0; j < res_common->cage->size; j++)
  {
    res_common->cage->values[2*j] = x[2*j];
    res_common->cage->values[2*j+1] = x[2*j+1];
  }

  return 0;
}


