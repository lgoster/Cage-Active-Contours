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
 * Get the affine coordinates for a given pixel for each point of the cage
 *
 */

void cac_get_affcoord_pixel(float v0x, float v0y, float *cage, int cage_size, float *affcoord)
{
  int j, j_prev, j_next, flag_norm, flag_j;
  float dx1, dy1, dx2, dy2;
  float norm1, norm2;
  float sum_weight;
  float *dd, *a = NULL;
  dd = (float*) malloc(sizeof(float)*cage_size);
  a = (float*) malloc(sizeof(float)*cage_size);

  sum_weight = 0.0;

  float sig = 0;
  flag_norm = 0;
  flag_j    = 0;
  for(j = 0; j < cage_size; j++)
  {
    j_next = j + 1;
    if (j_next >= cage_size)
      j_next = 0;

    dx1   = cage[2*j] - v0x;
    dy1   = cage[2*j+1] - v0y;
    norm1 = sqrt(dx1 * dx1 + dy1 * dy1);

    dx2   = cage[2*j_next] - v0x;
    dy2   = cage[2*j_next+1] - v0y;
    norm2 = sqrt(dx2 * dx2 + dy2 * dy2);

    dd[j] = norm1;
    sig = 1.0;
    if((dx1*dy2-dy1*dx2)<0) sig = -1.0;

    if ((norm1 <= 1e-05) || (norm2 <= 1e-05))
    {
      if (norm1 <= 1e-05) {
	flag_j    = j;
	flag_norm = 1;
	break;
      }
    }
    else
    {
      float val = (dx1*dx2+dy1*dy2)/(norm1*norm2);
      if(val>1) val = 1;
      if(val<-1) val = -1;
      a[j] = sig*tan(0.5*acos(val));

      if (!isfinite(a[j]))
	cac_error("error in computing affine coordinates for a pixel");
    }
  }

  if (flag_norm)
  {
    for(j = 0; j < cage_size; j++)
      affcoord[j] = 0.0;

    affcoord[flag_j] = 1.0;
  }
  else
  {
    for(j = 0; j < cage_size; j++)
    {
      j_prev = j - 1;
      if (j_prev < 0)
	j_prev = cage_size - 1;

      affcoord[j] = (a[j]+a[j_prev])/dd[j];

      sum_weight = sum_weight+affcoord[j];
    }

    for(j = 0; j < cage_size; j++)
    {
      affcoord[j] = affcoord[j]/sum_weight;
    }
  }

  free(a);
  free(dd);
}

/**
 *
 * Compute affine coordinates ordered per point, i.e.
 *  a[i][j] -- affine coordinate for point i, vertex j 
 *
 */

float **cac_get_affcoord(struct vector *points, struct vector *cagepoints)
{
  float v0x, v0y;
  float **affcoord;

  int k;
  int points_size, cage_size;

  points_size = points->size;
  cage_size = cagepoints->size;

  affcoord = (float **) malloc(sizeof(float **) * points_size);
  affcoord[0] = (float *) malloc(sizeof(float) * cage_size * points_size);

  // For each point of the list... 
  for(k = 0; k < points_size; k++)
  {
    affcoord[k] = affcoord[0] + k * cage_size;

    v0x = points->values[2*k];
    v0y = points->values[2*k+1];

    cac_get_affcoord_pixel(v0x, v0y, cagepoints->values, cagepoints->size, affcoord[k]);
  }

  return affcoord;
}

/**
 *
 *  Setup affine coordinates, etc. for a given iteration. This function is called
 *  once if non uniform sasmpling is used or at each iteration if uniform sampling
 *  is used. 
 *
 */

void cac_setup_omega1_omega2(double *xin, int init, struct resources_common *res_common)
{
  struct vector *lst_int, *lst_ext, *cage;
  struct image *contour, *mask, *mask_inverted;

  int n, n2, k, j;
  double x, y;

  n  = res_common->n;
  n2 = n / 2;

  /* Copy xin to floating point format */
  cage = cac_vector_alloc(n2, 2);
  for(k = 0; k < n; k++) cage->values[k] = xin[k];

  /* Recover contour from their affine coordinates */
  for(k = 0;k < res_common->contour->size; k++)
  {
    x = 0.0;
    y = 0.0;

    for(j = 0;j < res_common->cage->size; j++)
    {
      x += res_common->affcoord_contour[k][j]* xin[2*j];
      y += res_common->affcoord_contour[k][j]* xin[2*j+1];
    }

    res_common->contour->values[2*k] = x;
    res_common->contour->values[2*k+1] = y;
  }

  /* Allocate a mask to paint the contour */
  contour = cac_image_alloc(res_common->u->nrow, res_common->u->ncol);
  cac_image_clear(contour, 0.0);

  /* Paint contour */
  {
    int k_next;
    float *values = res_common->contour->values;
    j = res_common->contour->size - 1; 
    for(k = 0; k < j; k++)
    {
      k_next = k+1;
      cac_image_line_draw(contour, round(values[2*k]), round(values[2*k+1]), round(values[2*k_next]), round(values[2*k_next+1]), 1.0);
    }

    k_next = 0;
    cac_image_line_draw(contour, round(values[2*k]), round(values[2*k+1]), round(values[2*k_next]), round(values[2*k_next+1]), 1.0);
  }

  /* Fill the exterior contour */
  mask = cac_image_alloc(res_common->u->nrow, res_common->u->ncol);
  cac_image_clear(mask, 0.0);

  cac_filloutside(contour, mask, res_common);

  /* Free previously allocated memory */

  if (!init)
  {
    cac_vector_delete(res_common->pixels_ext);

    free(res_common->insideim_ext); 
    free(res_common->interpolated_u_ext);
    free(res_common->interpolated_ux_ext);
    free(res_common->interpolated_uy_ext);
    free(res_common->interpolated_dist_ext);

    free(res_common->affcoord_ext[0]);
    free(res_common->affcoord_ext);
  }

  /* Extract exterior pixels */

  lst_ext = cac_imagetovector(mask);

  if ((init) || (res_common->uniform_sampling))
    res_common->pixels_ext = lst_ext;
  else
    res_common->pixels_ext = cac_vector_alloc(lst_ext->size, 2); 

  res_common->insideim_ext          = (int *)   cac_xmalloc(lst_ext->size*sizeof(int));
  res_common->interpolated_u_ext    = (float *) cac_xmalloc(lst_ext->size*sizeof(float));
  res_common->interpolated_ux_ext   = (float *) cac_xmalloc(lst_ext->size*sizeof(float));
  res_common->interpolated_uy_ext   = (float *) cac_xmalloc(lst_ext->size*sizeof(float));
  res_common->interpolated_dist_ext = (float *) cac_xmalloc(lst_ext->size*sizeof(float));

  res_common->affcoord_ext = cac_get_affcoord(lst_ext, cage);

  /* Create an inverted copy of the mask we have obtained previously */

  mask_inverted = cac_image_alloc(res_common->u->nrow, res_common->u->ncol);
  cac_mask_invert(mask, mask_inverted);

  /* Free previously allocated memory */

  if (!init)
  {
    cac_vector_delete(res_common->pixels_int);

    free(res_common->insideim_int); 
    free(res_common->interpolated_u_int);
    free(res_common->interpolated_ux_int);
    free(res_common->interpolated_uy_int);
    free(res_common->interpolated_dist_int);

    free(res_common->affcoord_int[0]);
    free(res_common->affcoord_int);
  }

  /* Extract interior pixels */

  lst_int = cac_imagetovector(mask_inverted);

  if (res_common->uniform_sampling)
    res_common->pixels_int = lst_int;
  else
    res_common->pixels_int = cac_vector_alloc(lst_int->size, 2);

  res_common->insideim_int          = (int *)   cac_xmalloc(lst_int->size*sizeof(int));
  res_common->interpolated_u_int    = (float *) cac_xmalloc(lst_int->size*sizeof(float));
  res_common->interpolated_ux_int   = (float *) cac_xmalloc(lst_int->size*sizeof(float));
  res_common->interpolated_uy_int   = (float *) cac_xmalloc(lst_int->size*sizeof(float));
  res_common->interpolated_dist_int = (float *) cac_xmalloc(lst_ext->size*sizeof(float));

  res_common->affcoord_int = cac_get_affcoord(lst_int, cage);

  /* Done. We may delete non  necessary memory */

  if (!((init) | (res_common->uniform_sampling)))
  {
    cac_vector_delete(lst_ext);
    cac_vector_delete(lst_int);
  }

  cac_vector_delete(cage);

  cac_image_delete(contour);
  cac_image_delete(mask);
  cac_image_delete(mask_inverted);
}

/**
 *
 *  Only used for the non-uniform sampling approach. Assume that a the affine
 *  coordinates for a given set of initial interior and exterior pixels have
 *  been computed for a given cage. When the cage is deformed the "deformed"
 *  pixels can be recovered using this function. 
 *
 */

void cac_update_points_non_uniform_sampling(struct resources_common *res_common, double *xin)
{
  int k, j;
  int size_int = res_common->pixels_int->size;
  int size_ext = res_common->pixels_ext->size;

  double x, y;

  /* interior region */ 

  for(k = 0; k < size_int; k++)
  {
    x = 0.0;
    y = 0.0;

    for(j = 0; j < res_common->cage->size; j++)
    {
      x += res_common->affcoord_int[k][j] * xin[2*j];
      y += res_common->affcoord_int[k][j] * xin[2*j+1];
    }

    res_common->pixels_int->values[2*k]   = x;
    res_common->pixels_int->values[2*k+1] = y;
  }

  /* exterior region */
  for(k = 0; k < size_ext; k++)
  {
    x = 0.0;
    y = 0.0;

    for(j = 0; j < res_common->cage->size; j++)
    {
      x += res_common->affcoord_ext[k][j] * xin[2*j];
      y += res_common->affcoord_ext[k][j] * xin[2*j+1];
    }

    res_common->pixels_ext->values[2*k]   = x;
    res_common->pixels_ext->values[2*k+1] = y;
  }

  /* contour */
  for(k = 0;k < res_common->contour->size; k++)
  {
    x = 0.0;
    y = 0.0;

    for(j = 0;j < res_common->cage->size; j++)
    {
      x += res_common->affcoord_contour[k][j]* xin[2*j];
      y += res_common->affcoord_contour[k][j]* xin[2*j+1];
    }

    res_common->contour->values[2*k] = x;
    res_common->contour->values[2*k+1] = y;
  }
}

/**
 *
 *  This function interpolates the gray-level values of the pixels given a set
 *  of "deformed" pixel positions recovered with the function
 *  cac_update_points_non_uniform_sampling. Only used for non-uniform sampling
 *  since it is only for non-uniform sampling that non-integer pixel positions
 *  are used
 *
 */

void cac_interpolate_points_non_uniform_sampling(struct resources_common *res_common)
{
  int k;

  int size_eff_int = 0;
  int size_eff_ext = 0;

  int size_int = res_common->pixels_int->size;
  int size_ext = res_common->pixels_ext->size;

  /* interior */
  for( k = 0; k< size_int; k++ )
  {
    double x = res_common->pixels_int->values[2*k];
    double y = res_common->pixels_int->values[2*k+1];

    if (cac_inside_image_support(res_common->u, x, y))
    {
      res_common->insideim_int[k] = 1;
      res_common->interpolated_u_int[k]    = cac_get_bilinear_interpolation_image(res_common->u,x,y);
      res_common->interpolated_ux_int[k]   = cac_get_bilinear_interpolation_image(res_common->ux,x,y);
      res_common->interpolated_uy_int[k]   = cac_get_bilinear_interpolation_image(res_common->uy,x,y);
      size_eff_int++;
    }
    else
    {
      res_common->insideim_int[k] = 0;
    }
  }

  /* exterior */
  for( k = 0; k< size_ext; k++ )
  {
    double x = res_common->pixels_ext->values[2*k];
    double y = res_common->pixels_ext->values[2*k+1];

    if (cac_inside_image_support(res_common->u, x, y))
    {
      res_common->insideim_ext[k] = 1;
      res_common->interpolated_u_ext[k]    = cac_get_bilinear_interpolation_image(res_common->u,x,y);
      res_common->interpolated_ux_ext[k]   = cac_get_bilinear_interpolation_image(res_common->ux,x,y);
      res_common->interpolated_uy_ext[k]   = cac_get_bilinear_interpolation_image(res_common->uy,x,y);
      size_eff_ext++;
    }
    else
    {
      res_common->insideim_ext[k] = 0;
    }
  }

  res_common->size_eff_int = size_eff_int;
  res_common->size_eff_ext = size_eff_ext;
}

/**
 *
 *  This function is used in case uniform sampling is used. It just
 *  creates the lists evaluating the image at integer pixel positions.
 *
 */

void cac_interpolate_points_uniform_sampling(struct resources_common *res_common)
{
  int k;

  int size_eff_int = 0;
  int size_eff_ext = 0;

  int size_int = res_common->pixels_int->size;
  int size_ext = res_common->pixels_ext->size;

  /* interior */
  for( k = 0; k< size_int; k++ )
  {
    double x = res_common->pixels_int->values[2*k];
    double y = res_common->pixels_int->values[2*k+1];

    res_common->insideim_int[k] = 1;
    res_common->interpolated_u_int[k]    = cac_get_value_image(res_common->u,x,y);
    res_common->interpolated_ux_int[k]   = cac_get_value_image(res_common->ux,x,y);
    res_common->interpolated_uy_int[k]   = cac_get_value_image(res_common->uy,x,y);
    size_eff_int++;
  }

  /* exterior */
  for( k = 0; k< size_ext; k++ )
  {
    double x = res_common->pixels_ext->values[2*k];
    double y = res_common->pixels_ext->values[2*k+1];

    res_common->insideim_ext[k] = 1;
    res_common->interpolated_u_ext[k]    = cac_get_value_image(res_common->u,x,y);
    res_common->interpolated_ux_ext[k]   = cac_get_value_image(res_common->ux,x,y);
    res_common->interpolated_uy_ext[k]   = cac_get_value_image(res_common->uy,x,y);
    size_eff_ext++;
  }

  res_common->size_eff_int = size_eff_int;
  res_common->size_eff_ext = size_eff_ext;
}



