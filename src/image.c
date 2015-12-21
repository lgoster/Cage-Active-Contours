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
 *  Allocates a new image 
 *
 */

struct image *cac_image_new()
{
  struct image *image;

  image = (struct image *) cac_xmalloc(sizeof(struct image));

  image->nrow = image->ncol = 0;
  image->gray = NULL;

  return (image);
}

/**
 *
 *  Allocates the gray level array  
 *
 */

struct image *cac_image_alloc(
     int nrow,
     int ncol)
{
  int size;
  struct image *tmp;
  
  tmp = cac_image_new(); 

  size = nrow * ncol * sizeof(double);
  if (size <= 0)
  {
    fprintf(stderr, "Cannot allocate an image with zero or negative size\n");
    exit(1);
  }

  tmp->nrow = nrow;
  tmp->ncol = ncol;
  tmp->gray = (float *) cac_xmalloc(size);

  return tmp;
}

/**
 *
 *  Deallocates an image 
 *
 */

void cac_image_delete(
    struct image *image)

{
  if (image == NULL)
  {
    fprintf(stderr, "Cannot delete image: structure is NULL\n");
    exit(1);
  }

  if (image->gray != NULL) 
  {
    free(image->gray);
    image->gray = NULL;
  }
  else
    printf("image->gray is NULL! cannot delete\n");

  free(image);
}

/**
 *
 * Copy an image from in to out 
 *
 */

void cac_image_copy(
    struct image *in, 
    struct image *out)
{
  if ((!in) || (!out) || (!in->gray) || (!out->gray) 
      || (in->ncol != out->ncol) || (in->nrow != out->nrow)) 
  {
    fprintf(stderr, "Images cannot be copied: they are NULL or images are of different sizes\n");
    exit(1);
  }

  memcpy(out->gray, in->gray, sizeof(float) * in->ncol*in->nrow);
}


/**
 *
 * Initialize an image to a constant value 
 *
 */

void cac_image_clear(
  struct image *im,
  double value)
{
  float *p;
  int i, size;

  if ((!im) || (!im->gray))
    cac_error("Image is not allocated.\n");

  size = im->nrow * im->ncol;
  p    = im->gray;

  for(i = 0; i < size; i++, p++)
    *p = value;
}

/**
 *
 * Draw a line on an image with gray level c between (a0,b0) and (a1,b1) 
 *
 */

void cac_image_line_draw(struct image *image, int a0, int b0, int a1, int b1, float c)
{ 
  int bdx,bdy;
  int sx,sy,dx,dy,x,y,z;

  if ((!image) || (!image->gray)) 
    cac_error("NULL image struct or NULL gray plane\n");

  bdx = image->ncol;
  bdy = image->nrow;

  if (a0<0) a0=0; else if (a0>=bdx) a0=bdx-1;
  if (a1<0) 
  { 
    a1=0; 
  }
  else 
    if (a1>=bdx) 
    {
      a1=bdx-1; 
    }
  if (b0<0) b0=0; else if (b0>=bdy) b0=bdy-1;
  if (b1<0) 
  { 
    b1=0; 
  }
  else if (b1>=bdy) 
  { b1=bdy-1; 
  }

  if (a0<a1) { sx = 1; dx = a1-a0; } else { sx = -1; dx = a0-a1; }
  if (b0<b1) { sy = 1; dy = b1-b0; } else { sy = -1; dy = b0-b1; }
  x=0; y=0;

  if (dx>=dy) 
  {
    z = (-dx) / 2;
    while (abs(x) <= dx) 
    {
      image->gray[(y+b0)*bdx+x+a0] = c;
      x+=sx;
      z+=dy;
      if (z>0) { y+=sy; z-=dx; }
    } 
  }
  else 
  {
    z = (-dy) / 2;
    while (abs(y) <= dy) {
      image->gray[(y+b0)*bdx+x+a0] = c;
      y+=sy;
      z+=dx;
      if (z>0) { x+=sx; z-=dy; }
    }
  }
}

