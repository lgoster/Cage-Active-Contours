// This program is free software: you can use, modify and/or redistribute it
// under the terms of the simplified BSD License. You should have received a
// copy of this license along this program. If not, see
// <http://www.opensource.org/licenses/bsd-license.html>.
//
// Copyright (C) 2015, Lluis Garrido <lluis.garrido@ub.edu>
// All rights reserved.

#include "cac.h"

/*******************************************************************************/

/**
 *
 * Malloc with error checking
 *
 */

void *cac_xmalloc(size_t size)
{
  void *p = malloc(size);
  if (!p) 
    cac_error("xmalloc: out of memory!\n");

  return p;
}

/**
 *
 * Write a message to stdout and a file
 *
 */

void cac_print_msg(FILE *fp, const char *pFormat, ... )
{
  va_list args;
  va_start(args, pFormat);
  vprintf(pFormat, args );
  va_end( args );

  va_start(args, pFormat);
  if (fp)
    vfprintf(fp, pFormat, args);
  va_end( args );
}

/**
 *
 * Print error message
 *
 */

void cac_error(char *str)
{
  fprintf(stderr, "%s\n", str);
  exit(EXIT_FAILURE);
}

/*******************************************************************************/

/**
 *
 * Check if coordinate (x,y) is inside image support
 *
 */

int cac_inside_image_support(struct image *u, double x, double y)
{
  int max_x, max_y;

  max_x = u->ncol - 1;
  max_y = u->nrow - 1;

  if ((x < 0) || (x > max_x))
    return 0;

  if ((y <  0) || (y > max_y))
    return 0;

  return 1;
}

/**
 *
 * Project coordinate (x,y) inside image support 
 *
 */

void cac_project_to_imsupport(struct image *u, double *x, double *y)
{
  int max_x, max_y;

  max_x = u->ncol - 1;
  max_y = u->nrow - 1;

  if (*x < 0)
    *x = 0;
  else if (*x > max_x)
    *x = max_x;

  if (*y < 0)
    *y = 0;
  else if (*y > max_y)
    *y = max_y;
}

/**
 *
 * Get pixel value
 *
 */

double cac_get_value_image(struct image *u, double x, double y)
{
  int i, j;
  double ret;

  j = x;
  i = y;

  ret = u->gray[i * u->ncol + j];

  return ret;
}

/**
 * 
 * Bilinear interpolation
 *
 */

double cac_get_bilinear_interpolation_image(struct image *in, double xin, double yin)
{
  int l, k, offset, ncol;
  float *x0, *x1, *x2, *x3;
  double x, y, a, b, b_1, a_1, phi;

  x = xin;
  y = yin;

  /* Note that the next two lines do not use the floor function.
     Using the floor function is very slow! */

  l = x;
  k = y;

  a = x-l;
  b = y-k;

  ncol   = in->ncol;
  offset = k*ncol+l;

  a_1 = 1.0 - a;
  b_1 = 1.0 - b;

  /* This code has been optimized */

  x0 = x1 = in->gray + offset;
  x1++;
  x2 = x3 = x0 + ncol;
  x3++;

  if ((!a) || (!b))
  {
    if ((!a) && (!b))
      phi = *x0;
    else if (!a)
      phi = b_1 * (*x0) + b * (*x2);
    else
      phi = a_1 * (*x0) + a * (*x1);
  }
  else
    phi =
      b_1 * (a_1 * (*x0) + a * (*x1)) +
      b   * (a_1 * (*x2) + a * (*x3));

  return(phi);
}

/*******************************************************************************/

/** 
 *
 * Get contour of a mask store in img variable. Assume image is binary (0 and
 * any other value).  The interior contour is obtained.
 *
 */

#define EXPAND_SEQUENCE_CONTOUR  10000

void cac_contour_get_interior_contour(
    struct vector *contour, 
    struct image *img, 
    int conn)
{
  float p, *f, *values;
  int j, k, y, x, delta, npoints;
  int sizx, sizy, sizimg;

  int cy[8] = { -1,-1, 0, 1, 1, 1, 0,-1};
  int cx[8] = {  0, 1, 1, 1, 0,-1,-1,-1};

  /* Check if connectivity is 4 or 8 ... */

  delta = 0; /* Avoid complaining of compiler */

  if (conn == 4)
    delta = 2;
  else if (conn == 8)
    delta = 1;
  else
    cac_error("ERROR(contour_get_sequence): connectivity must be 4 or 8");

  /* Find first point of contour */

  sizx   = img->ncol;
  sizy   = img->nrow;
  sizimg = sizx * sizy;

  f = img->gray;
  j = 0;

  while ((f[j] == 0.0) && (j < sizimg)) j++;

  if (j == sizimg)
    cac_error("ERROR(contour_get_sequence): image has no non zero pixel");

  /* Allocate memory for sequence of points of the contour */

  values = (float *) cac_xmalloc(sizeof(float) * EXPAND_SEQUENCE_CONTOUR);

  values[0] = (float) (j % sizx); /* x coordinate */
  values[1] = (float) (j / sizx); /* y coordinate */

  /* Seek second contour point */

  k = 2;
  npoints = 0;

  do
  {
    x = (int) (values[0] + cx[k]); /* x coordinate */
    y = (int) (values[1] + cy[k]); /* y coordinate */

    p = f[y * sizx + x];

    if (p == 0)
    {
      k += delta;
      if (k >= 8) k -= 8;
    }
  }
  while ((p == 0) && (k != 0));

  if (p != 0)  /* Contour may have one pixel */
  {
    values[2] = (float) x; /* x coordinate */
    values[3] = (float) y; /* y coordinate */

    npoints = 1;

    while (1)
    {
      k += 6 - delta;
      if (k >= 8) k -= 8;

      do
      {
	k += delta;
	if (k >= 8) k -= 8;

	x = (int) (values[2*npoints]   + cx[k]);
	y = (int) (values[2*npoints+1] + cy[k]);

	p = f[y * sizx + x];
      }
      while (p == 0); 

      if ((x == values[0]) && (y == values[1])) 
      {
	if (conn == 4)
	  if (!((k == 6) && (f[(y + 1) * sizx + x] != 0)))
	    break;

	if (conn == 8)
	  if (!(((k == 6) || (k == 7)) && (f[(y + 1) * sizx + (x - 1)] != 0)))
	    break;
      }

      if (++npoints % EXPAND_SEQUENCE_CONTOUR == 0)
	values = realloc(values, (npoints + EXPAND_SEQUENCE_CONTOUR) * sizeof(float) * 2);

      values[2*npoints]   = (float) x; /* x coordinate */
      values[2*npoints+1] = (float) y; /* y coordinate */
    }
  }

  /* Set parameters */

  contour->size   = ++npoints;
  contour->dim    = 2;
  contour->values = realloc(values, 2 * npoints * sizeof(float));
}

/**
 *
 * Generate a 1D gaussian filter
 *
 */

void cac_gaussian_filter(float sigma, float ratio,
    int max_samples, float **filter, int *nsamples)
{
  int i, M;
  float sigma2, value1, value2, *p1, *p2;

  if (!(max_samples & 1))
    max_samples--;

  sigma2    = 2 * sigma * sigma;

  /* Compute how many samples we need. 
   Condition: gaussian(value1) >= value0 / ratio */     

  value1 = sqrt(- sigma2 * log(1.0 / ratio));

  /* And the total needed number of samples is ... */

  M         = (int) ceil(value1);
  *nsamples = M * 2 + 1;

  if (*nsamples > max_samples)
  {
    M = (max_samples - 1) >> 1;
    *nsamples = M * 2 + 1;
  }

  /* Allocate memory */

  *filter = (float *) cac_xmalloc(sizeof(float) * *nsamples);
  *filter = *filter + M;

  /* Compute samples */

  p1 = *filter;
  p2 = p1;

  *p1 = 1.0;

  p1++; 
  p2--;

  value2 = 1.0;

  for(i = 1; i <= M; i++, p1++, p2--)
  {
    value1  = exp(- SQR((double) i) / sigma2);
    value2 +=  2.0 * value1;

    *p1 = value1;
    *p2 = value1;
  }

  /* Let's normalize */

  p1 = *filter;
  p2 = p1;

  *p1 /= value2;

  p1++;
  p2--;

  for(i = 1; i <= M; i++, p1++, p2--)
  {
    *p1 /= value2;
    *p2 /= value2;
  } 
}

/*
 *
 * Used internally by cac_contour_blur
 *
 */

void cac_circular_convolution(float *input, int npoints, 
    float *filter, int length_fil, float *output)
{
  int M;
  float seq, *p_fil_begin, *p_fil, *p_fil_end;
  float *p_out, *p_out_end, *p_begin, *p_in, *p_in_begin, *p_in_end;

  /* Implementation of

     y[m] = SUM(filter[k] * x[m - k])

   */
  /* Initialize pointers */

  M = (length_fil - 1) / 2;

  /* For m = 0; y[0] = SUM(filter[k] * x[-k]); and for k = - M;
     we are multiplying filter[-M] * x[M]. p_begin points to x[M],
     each time m increases p_begin increases also. */

  p_begin  = input + (M % npoints);

  /* Where does the filter begin ? */

  p_fil_begin = filter - M;   /* filter[-M] */
  p_fil_end   = filter + M;   /* filter[M]  */

  /* Filter */

  p_out     = output;
  p_out_end = output + npoints;

  p_in_begin = input;
  p_in_end   = input + npoints;

  for(; p_out < p_out_end; p_out++)
  {
    seq   = 0;

    p_in  = p_begin;
    p_fil = p_fil_begin; 

    for(; p_fil <= p_fil_end; p_fil++)
    {
      seq += *p_fil * *p_in;

      if (--p_in < p_in_begin)     /* Circular rotation */
      {
	p_in = p_in_end;
	p_in--;
      }
    }

    *p_out = seq;

    if (++p_begin == p_in_end)
      p_begin = p_in_begin;
  }
}

/*
 *
 * Blur the contour
 *
 */

void cac_contour_blur(struct vector *contour)
{
  int i, size, M;
  float *filter;

  float *x, *x_out, *y, *y_out;

  /* Allocate for input */

  x = (float *) cac_xmalloc(contour->size * sizeof(float));
  y = (float *) cac_xmalloc(contour->size * sizeof(float));

  for(i = 0; i < contour->size; i++)
  {
    x[i] = contour->values[2*i];
    y[i] = contour->values[2*i+1];
  }

  /* Allocate for output */

  x_out = (float *) cac_xmalloc(contour->size * sizeof(float));
  y_out = (float *) cac_xmalloc(contour->size * sizeof(float));

  /* Filter */

  cac_gaussian_filter(1.5, 1000, 9, &filter, &size);

  cac_circular_convolution(x, contour->size, filter, size, x_out);
  cac_circular_convolution(y, contour->size, filter, size, y_out);

  /* Copy to output */

  for(i = 0; i < contour->size; i++)
  {
    contour->values[2*i] = x_out[i];
    contour->values[2*i+1] = y_out[i];
  }

  /* Free memory */

  free(x);
  free(y);
  free(x_out);
  free(y_out);

  /* Free filter. Not nice */

  M = (size - 1) / 2;
  free(filter-M);

}

/*******************************************************************************/

/**
 *
 *  Count how many elements of the mask image are different 
 *  from zero. This function is used to know if there is
 *  at least one pixel that belongs to the exterior/interior.
 *
 */

int cac_imagetovector_count_non_zero(struct image *mask)
{
  int count;
  float *p, *end;

  count = 0;

  end = mask->gray + mask->nrow * mask->ncol;
  for(p = mask->gray; p < end; p++)
  {
    if (*p != 0.0) count++; 
  }

  if (count == 0)
    cac_error("count_non_zero: count == 0\n"); 

  return count;
}

/**
 *
 *  Given a set of pixels that are marked in the
 *  mask image, extract them and put them into a vector.
 *
 */

struct vector *cac_imagetovector(struct image *mask)
{
  float *p;
  int n;

  struct vector *list;

  n = cac_imagetovector_count_non_zero(mask);
  list = cac_vector_alloc(n, 2); 

  p = mask->gray;

  int i = 0;
  int row, col;
  for(row = 0; row < mask->nrow; ++row) 
  {
    for(col = 0; col < mask->ncol; ++col, p++)
    {
      if (*p != 0) {
	list->values[2*i] = col;
	list->values[2*i+1] = row;
	i++;
      }
    }
  }

  if (i != n) 
    cac_error("imagetovector: i != n"); 

  return list;
}

/**
 *
 *  This function is called by cac_filloutside. The function implements a FIFO
 *  queue to propagate a set of pixels (stored in queue q) in the image.
 *  Propagation is performed to neighboring pixels if the mask is not marked as
 *  visited.
 *
 */

void cac_filloutside_queue(struct image *contour, struct image *mask, struct queue *q)
{
  struct point point;

  int k;
  int xp, yp;
  int new_xp, new_yp, new_offset;

  int dx[4] = {1, -1, 0, 0}, dy[4] = {0, 0, 1, -1};

  while (!cac_is_queue_empty(q))
  {
    cac_get_elem_queue((char *) &point, q);

    xp = point.x;
    yp = point.y;

    for(k = 0; k < 4; k++)
    {
      new_xp = xp + dx[k];
      new_yp = yp + dy[k];

      if ((new_xp < 0) || (new_yp < 0) || (new_xp >= mask->ncol)  || (new_yp >= mask->nrow))
	continue;

      new_offset = new_yp * mask->ncol + new_xp;

      if ((contour->gray[new_offset] == 0) && (mask->gray[new_offset] == 0))
      {
	mask->gray[new_offset] = 1.0;

	point.x = new_xp;
	point.y = new_yp;

	cac_put_elem_queue((char *) &point, q);
      }
    }
  }
}

/**
 *
 *  Given an image in which the contour is painted the function returns in mask
 *  the exterior pixels in mask image. The contour pixels are not painted in
 *  the mask image (i.e. the exterior pixels are strictly exterior)
 *
 */

void cac_filloutside(struct image *contour, struct image *mask, struct resources_common *res_common)
{
  int i, offset;

  struct queue *q;
  struct point point;

  q = cac_new_queue(sizeof(struct point), 4 * (mask->nrow + mask->ncol));

  /* We run through the border pixels to seek pixels outside
     the contour */

  /* Upper row */

  for(i = 0; i < mask->ncol; i++)
  {
    offset = i;

    if ((contour->gray[offset] == 0) && (mask->gray[offset] == 0))
    {
      mask->gray[offset] = 1.0;

      point.x = i;
      point.y = 0;

      cac_put_elem_queue((void *) &point, q);
    }
  }

  /* Lower row */

  for(i = 0; i < mask->ncol; i++)
  {
    offset = (mask->nrow - 1) * mask->ncol + i;
   
    if ((contour->gray[offset] == 0) && (mask->gray[offset] == 0))
    {
      mask->gray[offset] = 1.0;

      point.x = i;
      point.y = mask->nrow - 1;

      cac_put_elem_queue((void *) &point, q);
    }
  }

  /* Left border */

  for(i = 0; i < mask->nrow; i++)
  {
    offset = i * mask->ncol;

    if ((contour->gray[offset] == 0) && (mask->gray[offset] == 0))
    {
      mask->gray[offset] = 1.0;

      point.x = 0;
      point.y = i;

      cac_put_elem_queue((void *) &point, q);
    }
  }

  /* Right border */

  for(i = 0; i < mask->nrow; i++)
  {
    offset = i * mask->ncol + (mask->ncol - 1);

    if ((contour->gray[offset] == 0) && (mask->gray[offset] == 0))
    {
      mask->gray[offset] = 1.0;

      point.x = mask->ncol - 1;
      point.y = i;

      cac_put_elem_queue((void *) &point, q);
    }
  }

  /* Check if we have inserted any pixels */
  if (cac_is_queue_empty(q))
    cac_error("ERROR: initial queue is empty. Cannot fill outside pixels. Check the shape of the contour.\n");

  /* Perform holefilling */
  cac_filloutside_queue(contour, mask, q);

  /* Delete queue */
  cac_delete_queue(q);
}

/**
 * 
 *  Invert mask image
 *
 */

void cac_mask_invert(struct image *in, struct image *out)
{
  int i, size;
  float *pin, *pout;

  size = in->nrow * in->ncol;

  pin  = in->gray;
  pout = out->gray;

  for(i = 0; i < size; i++, pin++, pout++)
  {
    if (*pin)
      *pout = 0;
    else
      *pout = 1.0;
  }
}

/*******************************************************************************/


