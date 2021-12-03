#include "common.h"

void gnuplot(char *filename, char *xrange, char *yrange, int num_plot, ...)
{
    va_list ap;
    char cmd[1000], *p=cmd;
    int i, rc;

    p += sprintf(p, "gnuplot -p ");
    p += sprintf(p, "  -e \"set term x11 size 1200,400\"");
    p += sprintf(p, "  -e \"set xrange %s\"", xrange);
    p += sprintf(p, "  -e \"set yrange %s\"", yrange);
    p += sprintf(p, "  -e \"plot ");
    va_start(ap, num_plot);
    for (i = 0; i < num_plot; i++) {
        char *elements = va_arg(ap, char*);
        char *color    = va_arg(ap, char*);
        p += sprintf(p, "'%s' using %s with lines linewidth 2 linecolor rgb '%s'%s ",
                     filename, elements, color, i < num_plot-1 ? "," : "");
    }
    va_end(ap);
    p += sprintf(p, "\"");

    rc = system(cmd);
    if (rc < 0) {
        printf("ERROR failed to run gnuplot, %s\n", strerror(errno));
    }
}

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

