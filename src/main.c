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
 *  Main function
 *
 */

int main(int argc, char *argv[])
{
  int model, i;

  char *fname_image, *fname_mask, *fname_cage1, *fname_cage2;

  struct image *image;
  struct image *mask;
  struct vector *cage1, *cage2, *cage_out;

  if ((argc != 5) && (argc != 6)) {
    printf("Us: %s model image mask cage_init [cage_curr]\n", argv[0]);
    exit(1);
  }

  i = 1;
  model        = atoi(argv[i]); i++;
  fname_image  = argv[i]; i++;
  fname_mask   = argv[i]; i++;
  fname_cage1  = argv[i]; i++;
  fname_cage2  = NULL;

  if (argc == 6) {
    fname_cage2 = argv[i];
  }

  image = cac_pgm_read_image(fname_image);
  mask  = cac_pgm_read_image(fname_mask);
  cage1 = cac_read_points2d(fname_cage1);
  cage2 = NULL;

  if (fname_cage2) {
    cage2 = cac_read_points2d(fname_cage2);
    if (cage1->size != cage2->size) 
      cac_error("Both cages should have the same size.\n");
  }

  cage_out = cac_vector_alloc(cage1->size, 2);

  switch (model)
  {
    case ALGO_MEAN:
      cac_mean(image, mask, cage1, cage2, cage_out);
      break;
    case ALGO_GAUSSIAN:
      cac_gaussian(image, mask, cage1, cage2, cage_out);
      break;
    default:
      cac_error("No valid model selected.\n");
  }

  cac_vector_delete(cage_out);
  cac_vector_delete(cage1);
  if (cage2)
    cac_vector_delete(cage2);

  cac_image_delete(mask);
  cac_image_delete(image);
}
