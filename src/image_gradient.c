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
 *  Computes the 1st and 2nd order gradients of the image. For the mean
 *  and gaussian model only the 1st order gradient is enough. Just use
 *  the 2nd order gradient if you need it for other types of energies.
 *
 */

void cac_get_gradient_images(struct resources_common *res)
{
  struct vector *grad0_x, *grad1_x, *grad0_y, *grad1_y;

  /* Get filters for fist derivative */
  grad0_x = cac_get_filter_gradient(FILTER_TYPE, D_FIRST, FILTER_SIZE, 0);
  grad0_y = cac_get_filter_gradient(FILTER_TYPE, D_FIRST, FILTER_SIZE, 0);
  grad1_x = cac_get_filter_gradient(FILTER_TYPE, D_FIRST, FILTER_SIZE, 1);
  grad1_y = cac_get_filter_gradient(FILTER_TYPE, D_FIRST, FILTER_SIZE, 1);

  /* Compute gradients */
  res->ux = cac_gradient_images_finite_differences_convolution(res->u, grad1_x, grad0_y);
  res->uy = cac_gradient_images_finite_differences_convolution(res->u, grad0_x, grad1_y);

  /* Delete */
  cac_vector_delete(grad0_x);
  cac_vector_delete(grad0_y);
  cac_vector_delete(grad1_x);
  cac_vector_delete(grad1_y);

  /* Get filters for second derivative */
  grad0_x = cac_get_filter_gradient(FILTER_TYPE, D_SECOND, FILTER_SIZE, 0);
  grad0_y = cac_get_filter_gradient(FILTER_TYPE, D_SECOND, FILTER_SIZE, 0);
  grad1_x = cac_get_filter_gradient(FILTER_TYPE, D_SECOND, FILTER_SIZE, 2);
  grad1_y = cac_get_filter_gradient(FILTER_TYPE, D_SECOND, FILTER_SIZE, 2);

  /* Compute gradients */
  res->uxx = cac_gradient_images_finite_differences_convolution(res->u, grad1_x, grad0_y);
  res->uyy = cac_gradient_images_finite_differences_convolution(res->u, grad0_x, grad1_y);

  /* Delete */
  cac_vector_delete(grad0_x);
  cac_vector_delete(grad0_y);
  cac_vector_delete(grad1_x);
  cac_vector_delete(grad1_y);

  /* Get filters for second mixed derivative */
  grad1_x = cac_get_filter_gradient(FILTER_TYPE, D_SECOND, FILTER_SIZE, 1);
  grad1_y = cac_get_filter_gradient(FILTER_TYPE, D_SECOND, FILTER_SIZE, 1);

  /* Compute gradients */
  res->uxy = cac_gradient_images_finite_differences_convolution(res->u, grad1_x, grad1_y);

  /* De-allocate memory */
  cac_vector_delete(grad1_x);
  cac_vector_delete(grad1_y);
}

/**
 *
 * Performs convolution using finite differences. Used to compute gradient.
 *
 */


struct image *cac_gradient_images_finite_differences_convolution(
    struct image *in, 
    struct vector *filter_x, 
    struct vector *filter_y)
{
  struct image *out;

  double term_value, value;

  int i, j, k_i, k_j, d_i, d_j;
  int mid_sample_i, mid_sample_j;
  int nrow, ncol;

  nrow = in->nrow;
  ncol = in->ncol;

  /* Allocate memory for output image */

  out = cac_image_alloc(nrow, ncol);

  /* Apply filter */

  mid_sample_i = (filter_y->size / 2);
  mid_sample_j = (filter_x->size / 2);

  /* Filter the image */

  for(i = 0; i < nrow ; i++)
    for(j = 0; j < ncol ; j++)
    {
      value = 0.0;

      for(k_i = 0; k_i < filter_y->size; k_i++)
      {
	for(k_j = 0; k_j < filter_x->size; k_j++)
	{
	  /* Value of the pixel to take */

	  d_i = i + k_i - mid_sample_i;  
	  d_j = j + k_j - mid_sample_j; 

	  /* Check if the pixel is outside the image */

	  if (d_i < 0)
	    d_i = 0;

	  if (d_i >= nrow)
	    d_i = nrow - 1;

	  if (d_j < 0)
	    d_j = 0;

	  if (d_j >= ncol)
	    d_j = ncol - 1;

	  /* Filter */

	  term_value = in->gray[d_i * ncol + d_j] * filter_y->values[k_i] * filter_x->values[k_j];
	  value += term_value;
	}/*End for k_j*/
      } /* End for k_i */

      out->gray[i * ncol + j] = value;

    }/*End for (i,j)*/

  return out;
}


/**
 * 
 *  Get filter coefficients for gradient computation
 *
 */

struct vector *cac_get_filter_gradient(
    int filter_type, 
    int derivative_order, 
    int length, 
    int index)
{
  struct vector *filter;

  int i;
  double *my_filter; 

  /* Simoncelli filter taps, paper "Differentiation of multi-dimensional
   * signals" */

  double simoncelli_first_5taps[2][5] = 
  {{0.037659, 0.249153, 0.426375, 0.249153, 0.037659},
    {-0.109604, -0.276691, 0.0, 0.276691, 0.109604}};

  double simoncelli_first_7taps[2][7] = 
  {{0.005412,0.069591, 0.244560, 0.360875, 0.244560, 0.069591, 0.005412},
    {-0.019479, -0.123915, -0.193555, 0.0, 0.193555, 0.123915, 0.019479}};

  double simoncelli_first_9taps[2][9] =
  {{0.000738, 0.015530, 0.090260, 0.234469, 0.318007, 0.234469, 0.090260, 0.015530, 0.000738},
    {-0.003032, -0.035241, -0.118879, -0.144383, 0.0, 0.144383, 0.118879, 0.035241, 0.003032}};

  double simoncelli_first_11taps[2][11] = 
  {{0.000097, 0.003042, 0.026178, 0.103249, 0.223725, 0.287419, 0.223725, 0.103249, 0.026178, 0.003042, 0.000097},
    {-0.000440, -0.008083, -0.044978, -0.108831, -0.112871, 0.0, 0.112871, 0.108831, 0.044978, 0.008083, 0.000440}};

  double simoncelli_second_5taps[3][5] = 
  {{0.030320, 0.249724, 0.439911, 0.249724, 0.030320},
    {-0.104550, -0.292315, 0.0, 0.292315,  0.104550},
    {0.232905, 0.002668, -0.471147, 0.002668,  0.232905}};

  double simoncelli_second_7taps[3][7] = 
  {{0.004711, 0.069321, 0.245410, 0.361117, 0.245410, 0.069321, 0.004711},
    {-0.018708, -0.125376, -0.193091, 0.0, 0.193091, 0.125376, 0.018708},
    {0.055336, 0.137778, -0.056554, -0.273118, -0.056554, 0.137778, 0.055336}};

  double simoncelli_second_9taps[3][9] = 
  {{0.000734, 0.015591, 0.090295, 0.234405, 0.317952, 0.234405, 0.090295, 0.015591, 0.000734},
    {-0.003080, -0.035311, -0.118783, -0.144393, 0.0,  0.144393, 0.118783, 0.035311, 0.003080},
    {0.010316, 0.061361, 0.085644, -0.061329, -0.191986, -0.061329, 0.085644, 0.061361, 0.010316}};

  double simoncelli_second_11taps[3][11] = 
  {{0.000097, 0.003042, 0.026178, 0.103249, 0.223725, 0.287419, 0.223725, 0.103249, 0.026178, 0.003042, 0.000097},
    {-0.000440, -0.008083, -0.044978, -0.108831, -0.112871, 0.0, 0.112871, 0.108831, 0.044978, 0.008083, 0.000440},
    {0.001688, 0.017848, 0.057382, 0.053662, -0.059065, -0.143029, -0.059065, 0.053662, 0.057382, 0.017848, 0.001688}};

  /* 
   * Gaussian with a sigma = 0.96. The Matlab code to compute the coeficients
   * for the pre-filter, first order filter and second order filter is as follows:
   *
   * sigma = 0.96; pre = 1/sqrt(2*pi*sigma^2)*exp(-x.^2/(2*sigma^2)); 
   * disp(pre); disp(x/sigma^2 .* pre); disp((x.^2/sigma^4 - 1/sigma^2) .* pre)
   *
   */

  double gauss_3taps[3][3] =
  {{0.24156, 0.41556, 0.24156},
    {-0.26211, 0.00000, 0.26211},
    {0.022297, -0.450917, 0.022297}};

  double gauss_5taps[3][5] = 
  {{0.047442, 0.241557, 0.415565, 0.241557, 0.047442},
    {-0.10295, -0.26211, 0.00000, 0.26211, 0.10295},
    {0.171949, 0.022297, -0.450917, 0.022297, 0.171949}};

  double gauss_7taps[3][7] = 
  {{0.0031482, 0.0474416, 0.2415566, 0.4155649, 0.2415566, 0.0474416, 0.0031482},
    {-0.01025, -0.10295, -0.26211, 0.00000, 0.26211, 0.10295, 0.01025},
    {0.029943, 0.171949, 0.022297, -0.450917, 0.022297, 0.171949, 0.029943}};

  double gauss_9taps[3][9] =
  {{7.0586e-05, 3.1482e-03, 4.7442e-02, 2.4156e-01, 4.1556e-01, 2.4156e-01, 4.7442e-02, 3.1482e-03, 7.0586e-05},
    {-0.00031, -0.01025, -0.10295, -0.26211, 0.00000, 0.26211, 0.10295, 0.01025, 0.00031},
    {0.0012531, 0.0299434, 0.1719490, 0.0222972, -0.4509167, 0.0222972, 0.1719490, 0.0299434, 0.0012531}};

  double gauss_11taps[3][11] = 
  {{5.3474e-07, 7.0586e-05, 3.1482e-03, 4.7442e-02, 2.4156e-01, 4.1556e-01, 2.4156e-01, 4.7442e-02, 3.1482e-03, 7.0586e-05, 5.3474e-07},
    {-0.00000, -0.00031, -0.01025, -0.10295, -0.26211, 0.00000, 0.26211, 0.10295, 0.01025, 0.00031, 0.00000},
    {1.5160e-05, 1.2531e-03, 2.9943e-02, 1.7195e-01, 2.2297e-02, -4.5092e-01, 2.2297e-02, 1.7195e-01, 2.9943e-02, 1.2531e-03, 1.5160e-05}};

  /* Code begins here */

  my_filter = NULL;

  if (filter_type == FILTER_SIMONCELLI)
  {
    if (derivative_order == D_FIRST)
    {
      if (length == 5)
	my_filter = &simoncelli_first_5taps[index][0];
      else if (length == 7)
	my_filter = &simoncelli_first_7taps[index][0];
      else if (length == 9)
	my_filter = &simoncelli_first_9taps[index][0];
      else if (length == 11)
	my_filter = &simoncelli_first_11taps[index][0];
      else {
	fprintf(stderr, "no valid filter\n");
	exit(1);
      }
    }
    else if (derivative_order == D_SECOND)
    {
      if (length == 5)
	my_filter = &simoncelli_second_5taps[index][0];
      else if (length == 7)
	my_filter = &simoncelli_second_7taps[index][0];
      else if (length == 9)
	my_filter = &simoncelli_second_9taps[index][0];
      else if (length == 11)
	my_filter = &simoncelli_second_11taps[index][0];
      else {
	fprintf(stderr, "no valid filter\n");
	exit(1);
      }
    }
    else {
      fprintf(stderr, "no valid filter\n");
      exit(1);
    }
  }
  else if (filter_type == FILTER_GAUSS)
  {
    if (length == 3)
      my_filter = &gauss_3taps[index][0];
    else if (length == 5)
      my_filter = &gauss_5taps[index][0];
    else if (length == 7)
      my_filter = &gauss_7taps[index][0];
    else if (length == 9)
      my_filter = &gauss_9taps[index][0];
    else if (length == 11)
      my_filter = &gauss_11taps[index][0];
    else {
      fprintf(stderr, "no valid filter\n");
      exit(1);
    }
  }
  else {
    fprintf(stderr, "no valid filter\n");
    exit(1);
  }

  filter = cac_vector_alloc(length, 1);

  for(i = 0; i < length; i++)
    filter->values[i] = my_filter[i];

  return filter;
}


