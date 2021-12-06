#include "common.h"

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

    printf("CMD %s\n", cmd);

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

// xxx used?
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
