#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <stdarg.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include <errno.h>
#include <time.h>

void gnuplot(char *title, char *filename, char *xlabel, char *xrange, char *ylabel, char *yrange, 
             char **extra_cmds, ...);

int binary_search(double element, double *array, int max);
double max_array_val(int max_elem, ...);
