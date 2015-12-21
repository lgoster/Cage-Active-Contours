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
 * Creates a new empty vector structure 
 *
 */

struct vector *cac_vector_new()
{
  struct vector *vector;

  vector = (struct vector *) cac_xmalloc(sizeof(struct vector));

  vector->size = 0;
  vector->dim  = 0;
  vector->values = NULL;

  return (vector);
}

/** 
 *
 * Allocates a vector structure with N positions
 * and dimension dim
 *
 */

struct vector *cac_vector_alloc(
     int N,
     int dim)
{
  int size;
  struct vector *tmp;

  tmp = cac_vector_new();

  size = N * dim * sizeof(float);
  if (size <= 0)
  {
    fprintf(stderr, "Cannot allocate a vector with zero or negative size.\n");
    exit(1);
  }

  tmp->size = N;
  tmp->dim  = dim;
  tmp->values = (float *) cac_xmalloc(size);

  return tmp;
}

/**
 *
 * Deallocates a vector 
 *
 */

void cac_vector_delete(
    struct vector *vector)
{
  if (vector == NULL)
  {
    fprintf(stderr, "Cannot delete vector: structure is NULL\n");
    exit(1);
  }

  if (vector->values != NULL) 
  {
    free(vector->values);
    vector->values = NULL;
  }

  free(vector);
}

/** 
 *
 * Clear the array of a fsignal with value
 *
 */

void cac_vector_clear(
    struct vector *vector,
    float value)
{
  float *p;
  int i, j;

  if ((!vector) || (!vector->values)) 
  {
    fprintf(stderr, "Vector is not allocated.\n");
    exit(1);
  }

  p = vector->values;
  for (i = 0; i < vector->size; i++)
    for(j = 0; j < vector->dim; j++, p++)
      *p = value;
}






