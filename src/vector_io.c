#include "cac.h"

/**
 * 
 *  Reads a list fo 2D points from disc. Used to read the cage
 *  points from disc.
 *
 */

struct vector *cac_read_points2d(char *fname)
{
  FILE *fp;

  struct vector *vec;
  float *p_float;

  int i, num_lines;
  char line[100], *p;

  fp = fopen(fname, "r");
  if (!fp) {
    fprintf(stderr, "ERROR: could not open '%s' for reading.\n", fname);
    exit(1);
  }

  num_lines = 0;
  while (fgets(line, 100, fp)) num_lines++;

  rewind(fp);

  vec = cac_vector_alloc(num_lines, 2);
  cac_vector_clear(vec, 0.0);

  i = 0;
  p_float = vec->values;

  while (fgets(line, 100, fp))
  {
    p = strtok(line, " \t\n");
    if (!p) {
      fprintf(stderr, "ERROR: file '%s' does not contain enough data\n", fname);
      exit(1);
    }

    *p_float = atof(p);
    p_float++;

    p = strtok(NULL, " \t\n");
    if (!p) {
      fprintf(stderr, "ERROR: file '%s' does not contain enough data\n", fname);
      exit(1);
    }

    *p_float = atof(p);
    p_float++;

    i++;
  }

  if (i != num_lines) {
    fprintf(stderr, "ERROR: file '%s' does not contain enough lines\n", fname);
    exit(1);
  }

  fclose(fp);

  return vec;
}

/**
 *
 * Writes a list of 2D points to disc
 *
 */

void cac_write_points2d(struct vector *vector, int flag, char *fname)
{
  FILE *fp;

  int i;

  fp = fopen(fname, "w");
  if (!fp) 
    cac_error("cac_write_points2D: could not open file for writing.");

  for(i = 0; i < vector->size; i++) 
    fprintf(fp, "%f %f\n", vector->values[2*i], vector->values[2*i+1]);

  if (flag)
    fprintf(fp, "%f %f\n", vector->values[0], vector->values[1]);

  fclose(fp);
}

