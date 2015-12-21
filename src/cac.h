#ifndef CAC_H
#define CAC_H

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include <values.h>
#include <ctype.h>

#include "parameters.h"

/* 
 * Do not touch below unless you know what you are doing
 */

#include "tnc.h"

/*
 * Basic structures needed by the algorithm 
 */

/* The image structure */

struct image {
  int nrow;       /* Number of rows (dy) */
  int ncol;       /* Number of columns (dx) */
  float *gray;    /* The Gray level plane (may be NULL) */
};

/* A vector  */

struct vector {
  int size;
  int dim;
  float *values;
};

/* A point, for use in the queue */
struct point {
  int x;
  int y;
};


/* Queue */
struct queue {  /* Structure for fifo queue */
  int size_elem, nb_elem, expand_size, max_elem;
  char *b,*e; /* Begin and end of the queue */
  char *r, *w; /* Read and write pointers of the queue */
};


/* Common properties of all algorithms */

struct resources_common 
{
  /* Some parameters  */

  FILE *fp_log;
  int uniform_sampling;
  int project_gradient;

  /* Dimension of problem. Should be size of cage times 2 */
  int n; 

  /* Cage */
  struct vector *cage;

  /* Input image */
  struct image *u;

  /* Mask image */
  struct image *mask_in;

  /* Interior region and exterior pixel coordinates */
  struct vector *pixels_int;
  struct vector *pixels_ext;

  /* Affine coordinates */
  float **affcoord_int;
  float **affcoord_ext;

  /* Pixel coordinates of contour */
  struct vector *contour;

  /* Affine coordinates for mask interior reconstruction */
  float **affcoord_contour;

  /* Gradient images */
  struct image *ux,*uy,*uxx,*uxy,*uyy;

  /* Label to know if pixel is inside or outside the image support*/
  int *insideim_int;
  int *insideim_ext;

  /* Effective size */
  int size_eff_int;
  int size_eff_ext;

  /* Saved precalculated interpolate original values to speed up */
  float *interpolated_u_int;
  float *interpolated_ux_int;
  float *interpolated_uy_int;

  float *interpolated_u_ext;
  float *interpolated_ux_ext;
  float *interpolated_uy_ext;

  float *interpolated_dist_int;
  float *interpolated_dist_ext;

  /* Mean values */
  float mean_int;
  float mean_ext;

  /* Variances */
  float sigma2_int;
  float sigma2_ext;

  /* center of the contour, used for gradient projection */
  float cx, cy;
};

/* pi */
#ifndef M_PI  
#define M_PI		3.14159265358979323846
#endif

/* Algorithm type */
#define ALGO_MEAN       1
#define ALGO_GAUSSIAN   2

/* Other functions */
#define SQR(a)   ((a)*(a))
#define MAX(a,b) ((a)>(b)?(a):(b))
#define MIN(a,b) ((a)<(b)?(a):(b))

/* For gradient computation */
#define D_NONE                     0
#define D_FIRST                    1
#define D_SECOND                   2
#define D_SECOND_MIXED             3

/* Filter types */
#define FILTER_SIMONCELLI          0
#define FILTER_GAUSS               1

/* Filter to apply */
#define FILTER_SIZE                7
#define FILTER_TYPE                0

/* Other includes */
#include "prototypes.h"

#endif

