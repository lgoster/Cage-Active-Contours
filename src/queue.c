// This program is free software: you can use, modify and/or redistribute it
// under the terms of the simplified BSD License. You should have received a
// copy of this license along this program. If not, see
// <http://www.opensource.org/licenses/bsd-license.html>.
//
// Copyright (C) 2015, Lluis Garrido <lluis.garrido@ub.edu>
// All rights reserved.


#include "cac.h"

/**
 *  Queue functions. These functions are used for propagation purposes, i.e.
 *  to obtain the interior and exterior of the region.
 *
 */

/**
 *
 *  Create queue
 *
 */

struct queue *cac_new_queue(size_elem, expand_size)
  int size_elem, expand_size;
{ 
  int size;
  struct queue *q;

  size = sizeof(char) * size_elem * expand_size;
  
  q = (struct queue *) cac_xmalloc(sizeof(struct queue)); 

  q->b = (char *) cac_xmalloc(size);
  q->e = q->b + size; 

  q->r = q->w = q->b;

  q->size_elem   = size_elem;
  q->nb_elem     = 0;
  q->expand_size = expand_size;
  q->max_elem    = expand_size;

  return q;
}

/**
 *
 *  Delete queue
 *
 */

void cac_delete_queue(q)
  struct queue *q;
{
  if (!q)
    cac_error("[cac_delete_queue]: queue is not allocated.\n");

  free((char *) q->b);
  free(q);
}

/**
 *
 *  Check if queue is empty
 *
 */

int cac_is_queue_empty(q)
  struct queue *q;
{
  int rt;

  if (!q)
    cac_error("[cac_is_queue_empty]: queue is not allocated.\n");

  rt = 0;

  if (!q->nb_elem)
    rt = 1;

  return rt;
}

/**
 *
 *  Insert element into the queue
 *
 */

void cac_put_elem_queue(elem, q)
  char *elem;
  struct queue *q;
{
  if (!q)
    cac_error("[cac_put_queue]: queue is not allocated.\n");

  if (cac_is_queue_full(q))
    cac_expand_queue(q);

  /* Copy data into the queue */

  memcpy(q->w, elem, sizeof(char) * q->size_elem);

  /* Increment write pointer and number of elements */

  q->nb_elem++;
  q->w += sizeof(char) * q->size_elem;

  /* Check if we are out of range in the queue */

  if (q->w == q->e)
    q->w = q->b;
}

/**
 * 
 *  Get element from queue 
 *
 */

void cac_get_elem_queue(elem, q)
  char *elem; 
  struct queue *q;
{
  if (!q)
    cac_error("[cac_get_queue]: queue is not allocated.\n");

  if (cac_is_queue_empty(q))
    cac_error("[cac_get_queue]: queue is empty, no elements to extact.\n");

  /* Copy from queue to data */

  memcpy(elem, q->r, sizeof(char) * q->size_elem);

  /* Decrement number of elements, increment read pointer */

  q->nb_elem--;
  q->r += sizeof(char) * q->size_elem;

  /* Check if we are out of range in the queue */

  if (q->r == q->e)
    q->r = q->b;
}

/**
 *
 *  Get number of elements of the queue
 *
 */

int cac_get_queue_nb_elements(q)
  struct queue *q;
{
  if (!q)
    cac_error("[cac_get_queue_nb_elements]: queue is not allocated.\n");

  return (q->nb_elem);
}

/**
 *
 *  Check if queue is full. Used by cac_put_elem_queue
 *
 */

int cac_is_queue_full(q)
  struct queue *q;
{
  int rt;

  if (!q)
    cac_error("[cac_is_queue_full]: queue is not allocated.\n");

  rt = 0;

  if (q->nb_elem == q->max_elem)
    rt = 1;

  return rt;

}

/**
 *
 *  Expands size of the queue. Used by cac_put_elem_queue
 *
 */

void cac_expand_queue(q)
  struct queue *q;
{
  char *elem;
  struct queue *new_q;

  /* Create new queue */

  new_q = cac_new_queue(q->size_elem, q->max_elem + q->expand_size);

  /* Copy contents of the current queue to the new one */

  elem = (char *) cac_xmalloc(q->size_elem);

  while (!cac_is_queue_empty(q))
  {
    cac_get_elem_queue(elem, q);
    cac_put_elem_queue(elem, new_q);
  }

  free(elem);

  /* Now reassign pointers and free memory */

  free(q->b);

  q->b = new_q->b;
  q->e = new_q->e;
  q->r = new_q->r;
  q->w = new_q->w;

  q->max_elem = new_q->max_elem;
}

