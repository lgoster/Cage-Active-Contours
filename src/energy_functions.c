// This program is free software: you can use, modify and/or redistribute it
// under the terms of the simplified BSD License. You should have received a
// copy of this license along this program. If not, see
// <http://www.opensource.org/licenses/bsd-license.html>.
//
// Copyright (C) 2015, Lluis Garrido <lluis.garrido@ub.edu>
// All rights reserved.

#include "cac.h"

#include <sys/time.h>

/**
 *
 * Mean model. This is the energy function that uses the mean model. Updates
 * mean at each iteration, computes the gradient of the mean at each iteration
 * and does *not* divide by the number of pixels.
 *
 */

int cac_mean_evaluate(double *xin, double *f, double *deriv, void *instance)
{
  struct resources_common *res_common = (struct resources_common *) instance;

  int j, k, n;
  int size_int, size_ext;

  double B_energy = 0.0;

  n = res_common->n;

  double deriv_mean_int[n], deriv_mean_ext[n];

  /* set to zero the gradient */
  for(j = 0; j < n; j++)
  {
    deriv[j] = 0.0;
    deriv_mean_int[j] = 0.0;
    deriv_mean_ext[j] = 0.0;
  }

  if (res_common->uniform_sampling)
    cac_interpolate_points_uniform_sampling(res_common);
  else
  {
    cac_update_points_non_uniform_sampling(res_common, xin);
    cac_interpolate_points_non_uniform_sampling(res_common);
  }

  cac_compute_mean(res_common);

  /* ********** Compute energy *************** */

  size_int = res_common->pixels_int->size;
  size_ext = res_common->pixels_ext->size;

  // Energy interior
  double terme_int = 0;
  for(k = 0; k < size_int; k++)
  {
    if (res_common->insideim_int[k])
    {
      double I_p = res_common->interpolated_u_int[k];
      double dif = I_p - res_common->mean_int;

      terme_int = terme_int + 0.5*pow(dif,2.0);
    }
  }

  // Energy exterior
  double terme_ext = 0;
  for(k = 0; k < size_ext; k++)
  {
    if (res_common->insideim_ext[k])
    {
      double I_p = res_common->interpolated_u_ext[k];
      double dif = I_p - res_common->mean_ext;

      terme_ext = terme_ext + 0.5*pow(dif,2.0);
    }
  }

  // energy 
  B_energy = terme_int + terme_ext;

  /* ********** Compute gradient *************** */

  cac_compute_gradient_mean(res_common, deriv_mean_int, deriv_mean_ext);

  for(j = 0;j < res_common->cage->size; j++)
  {
    double terme_dx = 0;
    double terme_dy = 0;

    // Derivative of internal region with respect vj
    for( k = 0; k < size_int; k++)
    {
      if (res_common->insideim_int[k])
      {
	double I_p = res_common->interpolated_u_int[k];
	double dif = I_p - res_common->mean_int;

	terme_dx += dif * (res_common->interpolated_ux_int[k]*res_common->affcoord_int[k][j] - deriv_mean_int[2*j] );
	terme_dy += dif * (res_common->interpolated_uy_int[k]*res_common->affcoord_int[k][j] - deriv_mean_int[2*j+1] );
      }
    }

    deriv[2*j] = terme_dx;
    deriv[2*j+1] = terme_dy;

    // Derivative of external region with respect vj 
    terme_dx = 0.0;
    terme_dy = 0.0;

    for( k = 0; k < size_ext; k++)
    {
      if (res_common->insideim_ext[k])
      {
	double I_p = res_common->interpolated_u_ext[k];
	double dif = I_p - res_common->mean_ext;

	terme_dx += dif *(res_common->interpolated_ux_ext[k]*res_common->affcoord_ext[k][j] - deriv_mean_ext[2*j] );
	terme_dy += dif *(res_common->interpolated_uy_ext[k]*res_common->affcoord_ext[k][j] - deriv_mean_ext[2*j+1] );
      }
    }

    deriv[2*j] += terme_dx;
    deriv[2*j+1] += terme_dy;
  }

  *f = B_energy;

  return 0;
}

/**
 *
 * Gauss method. This is the energy function that uses the gaussian model. Mean
 * and variances are computed at each iteration inside and outside. We use the
 * log of the energy. We do not divide by the number of pixels. For the
 * gradient computation we assume that mean and variance are fixed.
 *
 */

int cac_gaussian_evaluate(double *xin, double *f, double *deriv, void *instance)
{
  struct resources_common *res_common = (struct resources_common *) instance;

  int j, k, n;

  double B_energy = 0.0;
  int size_int, size_ext; 

  n = res_common->n;

  /* set to zero the gradient */
  for(j = 0; j < n; j++)
    deriv[j] = 0.0;

  if (res_common->uniform_sampling)
  {
    cac_interpolate_points_uniform_sampling(res_common);
  }
  else
  {
    cac_update_points_non_uniform_sampling(res_common, xin);
    cac_interpolate_points_non_uniform_sampling(res_common);
  }

  cac_compute_mean(res_common);

  cac_compute_variance2(res_common, &(res_common->sigma2_int), &(res_common->sigma2_ext));

  /* ********** Compute energy *************** */

  size_int = res_common->pixels_int->size;
  size_ext = res_common->pixels_ext->size;

  // Energy computation for inner pixels
  double terme_int = 0;
  double E_1;

  for(k = 0; k < size_int; k++)
  {
    if (res_common->insideim_int[k])
    {
      double I_p = res_common->interpolated_u_int[k];
      terme_int += pow(I_p - res_common->mean_int,2.0);
    }
  }

  E_1 = 1.0/(2.0*res_common->sigma2_int) * terme_int;
  terme_int = size_int * log(sqrt(2.0 * M_PI * res_common->sigma2_int)) + E_1;

  // Energy computation for outer pixels
  double terme_ext = 0;
  double E_2;

  for(k = 0; k < size_ext; k++)
  {
    if (res_common->insideim_ext[k])
    {
      double I_p = res_common->interpolated_u_ext[k];
      terme_ext += pow(I_p - res_common->mean_ext,2.0);
    }
  }

  E_2 = 1.0/(2.0*res_common->sigma2_ext) * terme_ext;
  terme_ext = size_ext * log(sqrt(2.0 * M_PI * res_common->sigma2_ext)) + E_2;

  // Energy 
  B_energy = terme_int + terme_ext;

  /* **************** Compute gradient *************** */

  for(j = 0;j < res_common->cage->size; j++)
  {
    double terme_dx = 0.0;
    double terme_dy = 0.0;

    double deriv_int_vj_x = 0.0;
    double deriv_int_vj_y = 0.0;

    // Derivative of internal energy respect vj 
    for( k = 0; k< size_int; k++)
    {
      if (res_common->insideim_int[k])
      {
	double I_p = res_common->interpolated_u_int[k];
	double I_p_d_x = res_common->interpolated_ux_int[k];
	double I_p_d_y = res_common->interpolated_uy_int[k];

	deriv_int_vj_x +=  (I_p - res_common->mean_int) * I_p_d_x * res_common->affcoord_int[k][j];
	deriv_int_vj_y +=  (I_p - res_common->mean_int) * I_p_d_y * res_common->affcoord_int[k][j];
      }
    }
    deriv_int_vj_x = deriv_int_vj_x / res_common->sigma2_int;
    deriv_int_vj_y = deriv_int_vj_y / res_common->sigma2_int;

    // Derivative of internal term with respect vj
    terme_dx = deriv_int_vj_x;
    terme_dy = deriv_int_vj_y;

    // External region
    double deriv_ext_vj_x = 0.0;
    double deriv_ext_vj_y = 0.0;

    // Derivative of external energy with respect vj
    for( k = 0; k< size_ext; k++)
    {
      if (res_common->insideim_ext[k])
      {
	double I_p = res_common->interpolated_u_ext[k];
	double I_p_d_x = res_common->interpolated_ux_ext[k];
	double I_p_d_y = res_common->interpolated_uy_ext[k];
	deriv_ext_vj_x +=  (I_p - res_common->mean_ext) * I_p_d_x * res_common->affcoord_ext[k][j];
	deriv_ext_vj_y +=  (I_p - res_common->mean_ext) * I_p_d_y * res_common->affcoord_ext[k][j];
      }
    }
    deriv_ext_vj_x = deriv_ext_vj_x / res_common->sigma2_ext;
    deriv_ext_vj_y = deriv_ext_vj_y / res_common->sigma2_ext;

    // Derivative of external term with respect vj
    terme_dx += deriv_ext_vj_x;
    terme_dy += deriv_ext_vj_y;

    // Update derivatives
    deriv[2*j] += terme_dx;
    deriv[2*j+1] += terme_dy;
  }

  // Energy is returned through argument
  *f = B_energy;

  return 0;
}

/**
 *
 * Given a set of pixels of interior and exterior pixels (stored in
 * res_common), compute the mean associated to the interior and exterior model.
 *
 */

void cac_compute_mean(struct resources_common *res_common)
{
  int k;

  int size_int = res_common->pixels_int->size;
  int size_ext = res_common->pixels_ext->size;

  //Establish to 0 the mean values
  double mean_int = 0;
  double mean_ext = 0;

  /* mean interior */
  for( k = 0; k< size_int; k++)
    if (res_common->insideim_int[k])
      mean_int += res_common->interpolated_u_int[k];

  mean_int = mean_int / (double) res_common->size_eff_int;

  /* mean exterior */
  for( k = 0; k< size_ext; k++ )
    if (res_common->insideim_ext[k])
      mean_ext += res_common->interpolated_u_ext[k];

  mean_ext = mean_ext / (double) res_common->size_eff_ext;

  res_common->mean_int = mean_int;
  res_common->mean_ext = mean_ext;
}

/**
 *
 * Given a set of pixels of interior and exterior pixels (stored in
 * res_common), compute the gradient of the mean associated to the
 * interior and exterior model.
 *
 */

void cac_compute_gradient_mean(struct resources_common *res_common, double *deriv_int, double *deriv_ext)
{
  int j, k;
  int size_int = res_common->pixels_int->size;
  int size_ext = res_common->pixels_ext->size;

  for(j = 0; j < res_common->cage->size; j++)
  {
    // Derivative of internal mean with respect vj
    double terme_int_dx = 0.0;
    double terme_int_dy = 0.0;

    double terme_ext_dx = 0.0;
    double terme_ext_dy = 0.0;

    for( k = 0; k < size_int; k++)
    {
      if (res_common->insideim_int[k])
      {
	terme_int_dx += res_common->interpolated_ux_int[k] * res_common->affcoord_int[k][j];
	terme_int_dy += res_common->interpolated_uy_int[k] * res_common->affcoord_int[k][j];
      }
    }

    terme_int_dx = terme_int_dx / (double) res_common->size_eff_int;
    terme_int_dy = terme_int_dy / (double) res_common->size_eff_int;

    // Derivative of external ext with respect vj
    for( k = 0; k < size_ext; k++)
    {
      if (res_common->insideim_ext[k])
      {
	terme_ext_dx +=  res_common->interpolated_ux_ext[k] * res_common->affcoord_ext[k][j];
	terme_ext_dy +=  res_common->interpolated_uy_ext[k] * res_common->affcoord_ext[k][j];
      }
    }

    terme_ext_dx = terme_ext_dx / (double) res_common->size_eff_ext;
    terme_ext_dy = terme_ext_dy / (double) res_common->size_eff_ext;

    deriv_int[2*j]   = terme_int_dx;
    deriv_int[2*j+1] = terme_int_dy;

    deriv_ext[2*j]   = terme_ext_dx;
    deriv_ext[2*j+1] = terme_ext_dy;
  }
}

/**
 *
 * Given a set of pixels of interior and exterior pixels (stored in
 * res_common), compute the variance (sigma^2) associated to the
 * interior and exterior model.
 *
 */

void cac_compute_variance2(struct resources_common *res_common, float *s2int, float *s2ext)
{
  int k;
  int size_int = res_common->pixels_int->size;
  int size_ext = res_common->pixels_ext->size;

  // Calculate variance of the interior.
  double sigma2_int = 0.0;
  for( k = 0; k< size_int; k++ )
  {
    if (res_common->insideim_int[k])
    {
      double dif = res_common->interpolated_u_int[k] - res_common->mean_int;
      sigma2_int += pow(dif,2.0);
    }
  }
  sigma2_int = sigma2_int / (double) res_common->size_eff_int;

  // Calculate variance of the exterior.
  double sigma2_ext = 0.0;
  for( k = 0; k< size_ext; k++ )
  {
    if (res_common->insideim_ext[k])
    {
      double dif = res_common->interpolated_u_ext[k] - res_common->mean_ext;
      sigma2_ext += pow(dif,2.0);
    }
  }
  sigma2_ext = sigma2_ext / (double) res_common->size_eff_ext;

  *s2int = sigma2_int;
  *s2ext = sigma2_ext;
}

/**
 *
 * This function is used to avoid crossings of edges of the cage. It computes,
 * for each vertex, its distance to other edges and vertexs. The exact
 * algorithm is explained in the paper. 
 *
 */

double cac_compute_dist_constraints(struct resources_common *res_common, double *xin)
{
  double mindist, r;
  int j, k;

  double dist_cost = 0.0;

  mindist = 2.0 * BANDSIZE;

  for(j = 0; j < res_common->cage->size; j++)
  {
    if (!ALLOWCROSSINGS)
    {
      double x0 = xin[2 * j];
      double y0 = xin[2 * j + 1];

      for(k = 0; k < res_common->cage->size; k++)
      {
	if (k == j)
	  continue;

	double x1 = xin[2 * k];
	double y1 = xin[2 * k + 1];

	double dx = x0 - x1;
	double dy = y0 - y1;
	r = sqrt(dx * dx + dy * dy);

	if (r <= mindist)
	  dist_cost += pow(mindist - r, 2.0);

	int k_next = k+1;

	if (k_next == res_common->cage->size)
	  k_next = 0;

	if (k_next == j)
	  continue;

	double x2 = xin[2 * k_next];
	double y2 = xin[2 * k_next + 1];

	double dist      = sqrt(pow(x2 - x1, 2.0) + pow(y2 - y1, 2.0));
	double dist_edge = fabs((y2 - y1) * (x0 - x1) - (x2 - x1) * (y0 - y1))/dist;
	double proj_edge = ((x2 - x1) * (x0 - x1) + (y2 - y1) * (y0 - y1))/dist;

	if ((proj_edge >= 0) && (proj_edge <= dist) && (dist_edge <= mindist)) 
	  dist_cost += pow(mindist - dist_edge, 2.0);
      }
    }
  }

  return dist_cost;
}



