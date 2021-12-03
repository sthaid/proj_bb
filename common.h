#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <stdarg.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include <errno.h>

int binary_search(double element, double *array, int max);
void gnuplot(char *filename, char *xrange, char *yrange, int num_plot, ...);
