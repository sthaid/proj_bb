#include "common.h"

// -----------------  RANDOM VALUES FROM PROB DISTRIBUTION  -------------

typedef struct {
    double  start;
    double  end;
    double  delta;
    int     max;
    double *cum_prob;
    double  (*func)(double);
} handle_t;

void *probdist_create(double func(double), double start, double end)
{
    handle_t *hndl;
    double sum_p = 0;
    int max = 10000; //xxx
    double delta = (end-start)/max;

    hndl = malloc(sizeof(handle_t));
    hndl->start = start;
    hndl->end = end;
    hndl->delta = delta;
    hndl->max = max; 
    hndl->cum_prob = calloc(max, sizeof(double));
    hndl->func = func;

    for (int idx = 0; idx < max; idx++) {
        double p = func(start + idx * delta);
        sum_p += p * delta;
        hndl->cum_prob[idx] = sum_p;
    }

    //printf("max = %d\n", max);
    //printf("sum_p = %f\n", sum_p);
    assert(sum_p > 0.999 && sum_p <= 1);

    return hndl;
}

void probdist_destroy(void *hndl_arg)
{
    handle_t *hndl = hndl_arg;
    free(hndl->cum_prob);
    free(hndl);
}

double probdist_get_value(void *hndl_arg)
{
    handle_t *hndl = hndl_arg;
    double rand_range_0_1;
    int idx;

    do {
        rand_range_0_1 = (double)random() / ((double)RAND_MAX + 1);
    } while (rand_range_0_1 >= hndl->cum_prob[hndl->max-1]);

    idx = binary_search(rand_range_0_1, hndl->cum_prob, hndl->max);
    if (idx < 0 || idx >= hndl->max-1) {
        printf("ERROR get, bad idx %d\n", idx);
        exit(1);
    }

    return hndl->start + idx * hndl->delta;
}

void probdist_test(void *hndl_arg)
{
    handle_t *hndl = hndl_arg;
    int *histogram, histogram_cnt=0;
    
    FILE *fp;
    char title[100];

    #define MAX_TEST 10000000

    histogram = calloc(hndl->max, sizeof(int));

    for (int i = 0; i < MAX_TEST; i++) {
        double value = probdist_get_value(hndl_arg);
        int idx = (value - hndl->start) * (hndl->max / (hndl->end - hndl->start));
        if (idx < 0 || idx >= hndl->max) {
            printf("ERROR idx=%d max=%d\n", idx, hndl->max);
            continue;
        }
        histogram[idx]++;
        histogram_cnt++;
    }

    fp = fopen("test.dat", "w");
    double sum_p1=0, sum_p2=0;
    for (int i = 0; i < hndl->max; i++) {
        double value = hndl->start + i * hndl->delta;
        double p1 = (double)histogram[i] / histogram_cnt;
        double p2 = hndl->func(value) * hndl->delta;  // xxx func name
        fprintf(fp, "%e %0.9f %0.9f\n", value, p1, p2);
        sum_p1 += p1;
        sum_p2 += p2;
    }
    fclose(fp);

    if (sum_p1 < 0.999 || sum_p1 > 1.001 || sum_p2 < 0.999 || sum_p2 > 1.001) {
        printf("probdist_test ERROR: sum_p1=%0.20f sum_p2=%0.6f\n", sum_p1, sum_p2);
    }

    sprintf(title, "probdist_test");
    gnuplot(title, "test.dat", 
            "Value", "[*:*]", 
            "Probability", "[*:*]", 
            NULL,
            "", "1:2", "green", 
            "", "1:3", "blue",
            NULL, NULL, NULL);
}

// -----------------  GNUPLOT  ------------------------------------------

void gnuplot(char *title, char *filename, char *xlabel, char *xrange, char *ylabel, char *yrange, 
             char **extra_cmds, ...)
{
    va_list ap;
    char cmd[1000], *p=cmd;
    int rc, i;

    p += sprintf(p, "gnuplot-qt -p ");
    p += sprintf(p, "  -e \"set term qt size 1200,400\"");
    p += sprintf(p, "  -e \"set title '%s'\"", title);
    p += sprintf(p, "  -e \"set xlabel '%s'\"", xlabel);
    p += sprintf(p, "  -e \"set ylabel '%s'\"", ylabel);
    p += sprintf(p, "  -e \"set xrange %s\"", xrange);
    p += sprintf(p, "  -e \"set yrange %s\"", yrange);
    for (i = 0; extra_cmds && extra_cmds[i]; i++) {
        p += sprintf(p, "  -e \"%s\"", extra_cmds[i]);
    }
    p += sprintf(p, "  -e \"plot ");
    va_start(ap, extra_cmds);
    while (true) {
        char *title    = va_arg(ap, char*);
        if (title == NULL) break;
        char *elements = va_arg(ap, char*);
        char *color    = va_arg(ap, char*);
        p += sprintf(p, "'%s' using %s title '%s' with lines linewidth 2 linecolor rgb '%s', ",
                     filename, elements, title, color);
    }
    va_end(ap);
    p += sprintf(p, "\"");

    //printf("CMD %s\n", cmd);

    rc = system(cmd);
    if (rc < 0) {
        printf("ERROR failed to run gnuplot, %s\n", strerror(errno));
    }
}

// -----------------  MISC  ---------------------------------------------

int binary_search(double element, double *array, int max)
{
    int start=0, end=max-1, middle;

    while (start <= end) {
        middle = start + (end - start) / 2;
        if (middle == max-1) break;

        if (element < array[middle]) {
            end = middle-1;
        } else if (element >= array[middle+1]) {
            start = middle+1;
        } else {
            return middle;
        }
    }

    return -1;
}

double max_array_val(int max_elem, ...)
{
    double max = -INFINITY, *array;
    va_list ap;
    int i;

    va_start(ap, max_elem);
    while (true) {
        array = va_arg(ap, double*);
        if (array == NULL) break;
        for (i = 0; i < max_elem; i++) {
            if (array[i] > max) max = array[i];
        }
    }
    va_end(ap);

    return max;
}

uint64_t microsec_timer(void)
{
    struct timespec ts;

    clock_gettime(CLOCK_MONOTONIC,&ts);
    return  ((uint64_t)ts.tv_sec * 1000000) + ((uint64_t)ts.tv_nsec / 1000);
}

void random_vector(double magnitude, double * x, double * y, double * z)
{
    double x_try, y_try, z_try, hypot, f;

    // compute x/y/z_try within a spherical shell 
    while (true) {
        x_try = random() - (RAND_MAX/2.);
        y_try = random() - (RAND_MAX/2.);
        z_try = random() - (RAND_MAX/2.);
        hypot = hypotenuse(x_try,y_try,z_try);
        if (hypot >= (RAND_MAX/10.) && hypot <= (RAND_MAX/2.)) {
            break;
        }
    }

    // scale the random vector to the caller's specified magnitude
    f = magnitude / hypot;
    *x = x_try * f;
    *y = y_try * f;
    *z = z_try * f;
}

double hypotenuse(double x, double y, double z)
{
    return sqrt(x*x + y*y + z*z);
}

