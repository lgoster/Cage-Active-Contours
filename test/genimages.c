/**
 * 
 *  A small utility to generate gnuplot commands that allow
 *  to see the iterations performed by the CAC algorithm.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

void print_commands(int i, int ncol, int nrow)
{
  char *image, *cadena, *coordext, *coordint, *contour, *cage, *gradient;

  image    = malloc(1000);
  cadena   = malloc(1000);
  coordext = malloc(1000);
  coordint = malloc(1000);
  contour  = malloc(1000);
  cage     = malloc(1000);
  gradient = malloc(1000);


  sprintf(cadena, "set out \"gnuplot_%03d.png\"\n", i);
  printf(cadena);

  sprintf(image, "\"image.img\" binary array=(%d,%d) with image", ncol, nrow);
  sprintf(coordext, "\"coord_ext_%03d.txt\" with dots", i);
  sprintf(coordint, "\"coord_int_%03d.txt\" with dots", i);
  sprintf(contour,  "\"contour_%03d.txt\" with lines lw 2", i);
  sprintf(cage,  "\"cage_%03d.txt\" with lines lw 2", i);
  sprintf(gradient,  "\"gradient_%03d.txt\" with vectors", i);
  sprintf(cadena, "plot [0:%d][0:%d] %s, %s, %s\n", ncol-1, nrow-1, image, cage, contour);
  printf(cadena);

  free(image);
  free(cadena);
  free(coordext);
  free(coordint);
  free(contour);
  free(cage);
  free(gradient);
}

int main(int argc, char **argv)
{
  char *cadena;
  int i, nrow, ncol;

  if (argc != 3) {
    printf("Us: %s ncol nrow\n", argv[0]);
    exit(1);
  }

  ncol = atoi(argv[1]);
  nrow = atoi(argv[2]);

  cadena   = malloc(1000);

  printf("set palette gray\n");
  printf("set term png crop\n");
  printf("unset xtics\n");
  printf("unset ytics\n");
  printf("unset key\n");
  printf("unset colorbox\n");
  printf("set yrange [*:*] reverse\n");

  i = 0;
  while (1)
  {
    sprintf(cadena, "cage_%03d.txt", i);
    if (access(cadena, F_OK) == -1)
      break;

    print_commands(i, ncol, nrow);
    i++;
  }

  free(cadena);
}
