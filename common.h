#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <stdarg.h>
#include <stdint.h>
#include <unistd.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include <errno.h>
#include <time.h>

void *probdist_create(double func(double), double start, double end);
void probdist_destroy(void *hndl_arg);
double probdist_get_value(void *hndl_arg);
void probdist_test(void *hndl_arg);

void gnuplot(char *title, char *filename, char *xlabel, char *xrange, char *ylabel, char *yrange, 
             char **extra_cmds, int num_plot, ...);

int binary_search(double element, double *array, int max);
double max_array_val(int max_elem, ...);
uint64_t microsec_timer(void);
void random_vector(double magnitude, double * x, double * y, double * z);
double hypotenuse(double x, double y, double z);

static inline char *bool2str(bool x) {return x ? "true" : "false";}
