// This program is free software: you can use, modify and/or redistribute it
// under the terms of the simplified BSD License. You should have received a
// copy of this license along this program. If not, see
// <http://www.opensource.org/licenses/bsd-license.html>.
//
// Copyright (C) 2013, Lluis Garrido <lluis.garrido@ub.edu>
// All rights reserved.

#include "cac.h"

#include <string.h>

/*
 *
 * CAC mean energy interface
 *
 */

void cac_mean(struct image *img, struct image *mask_in, struct vector *cage_init, struct vector *cage_curr, struct vector *cage_out)
{
  struct resources_common res_common;

  cac_common(ALGO_MEAN, img, mask_in, cage_init, cage_curr, cage_out, &res_common);
}

/*
 *
 * CAC gaussian energy interface
 *
 */

void cac_gaussian(struct image *img, struct image *mask_in, struct vector *cage_init, struct vector *cage_curr, struct vector *cage_out)
{
  struct resources_common res_common;

  cac_common(ALGO_GAUSSIAN, img, mask_in, cage_init, cage_curr, cage_out, &res_common);
}

/*
 *
 *  Common interface for all energies
 *
 */

void cac_common(
    int algorithm, 
    struct image *img, 
    struct image *mask_in, 
    struct vector *cage_init,
    struct vector *cage_curr,
    struct vector *cage_out, 
    struct resources_common *res_common)
{
  struct vector *contour_init;

  int j, n;
  double *x, *g;

  if (cage_init->dim != 2)
    cac_error("in function cac_common cage dimension should be two");

  /* Log */
  char fname[256];
  sprintf(fname, "%s/%s", LOGDIR, LOGFILE);
  printf("%s\n", fname);
  res_common->fp_log = fopen(fname, "w");
  if (!res_common->fp_log)
    cac_error("Could not open log file. Does the log directory exist?\n");

  sprintf(fname, "%s/image.img", LOGDIR);
  cac_write_raw_image(img, fname);

  cac_print_msg(res_common->fp_log, "Performing initializations...\n");

  /* Assign input data to resources */
  res_common->cage = cage_init;
  res_common->u = img;
  res_common->mask_in = mask_in;

  /* Compute center of cage. Used to project gradien in first phase */
  res_common->cx = 0.0;
  res_common->cy = 0.0;

  for(j = 0; j < cage_init->size; j++)
  {
    res_common->cx += cage_init->values[2*j];
    res_common->cy += cage_init->values[2*j+1];
  }

  res_common->cx /= (float) cage_init->size;
  res_common->cy /= (float) cage_init->size;

  /* Extract contour and Compute associated affine coordinates */

  contour_init = cac_vector_new();
  cac_contour_get_interior_contour(contour_init, mask_in, 8);  

  cac_contour_blur(contour_init);

  res_common->affcoord_contour = cac_get_affcoord(contour_init, cage_init);

  /* Allocate memory for deformed contour */

  res_common->contour = cac_vector_alloc(contour_init->size, 2);

  /* Compute gradient images */

  cac_get_gradient_images(res_common);

  /* Variable x will hold the current cage, g the gradient */

  n  = cage_init->size * 2;
  res_common->n = n;

  x = (double *) malloc(n * sizeof(double));
  g = (double *) malloc(n * sizeof(double));

  if (cage_curr)
  {
    for(j = 0; j < n; j++)
      x[j] = cage_curr->values[j];
  }
  else
  {
    for(j = 0; j < n; j++)
      x[j] = cage_init->values[j];
  }

  for(j = 0; j < n; j++)
    *g = 0;

  res_common->uniform_sampling = 1;
  res_common->project_gradient = 1;

  cac_setup_omega1_omega2(x, 1, res_common);

  /* Now proceed to optimization */

  switch (algorithm)
  {
    case ALGO_MEAN:
      cac_print_msg(res_common->fp_log, "Model: mean energy\n");
      cac_newton_algorithm(n, x, cac_mean_evaluate, (void *) res_common);
      break;
    case ALGO_GAUSSIAN:
      cac_print_msg(res_common->fp_log, "Model: gaussian energy\n");
      cac_newton_algorithm(n, x, cac_gaussian_evaluate, (void *) res_common);
      break;
    default:
      cac_error("No valid algorithm in cac_common");
  }

  /* Store output cage */

  for(j = 0; j < cage_init->size; j++)
  {
    cage_out->values[2*j] = x[2*j];
    cage_out->values[2*j+1] = x[2*j+1];
  }

  /* Free memory */

  fclose(res_common->fp_log);

  cac_vector_delete(contour_init);
  cac_vector_delete(res_common->contour);

  cac_vector_delete(res_common->pixels_int);
  cac_vector_delete(res_common->pixels_ext);

  cac_image_delete(res_common->uyy); 
  cac_image_delete(res_common->uxy); 
  cac_image_delete(res_common->uxx);
  cac_image_delete(res_common->uy); 
  cac_image_delete(res_common->ux);

  free(res_common->interpolated_dist_int);
  free(res_common->interpolated_u_int);
  free(res_common->interpolated_ux_int);
  free(res_common->interpolated_uy_int);
  free(res_common->insideim_int); 

  free(res_common->interpolated_dist_ext);
  free(res_common->interpolated_u_ext);
  free(res_common->interpolated_ux_ext);
  free(res_common->interpolated_uy_ext);
  free(res_common->insideim_ext); 

  free(res_common->affcoord_int[0]);
  free(res_common->affcoord_int);

  free(res_common->affcoord_ext[0]);
  free(res_common->affcoord_ext);

  free(res_common->affcoord_contour[0]);
  free(res_common->affcoord_contour);

  free(x);
  free(g);
}


