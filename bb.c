#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>

#define   T    6000.
#define   K    1.380649E-23
#define   C    299792458.
#define   h    6.62607004e-34
#define   hbar 1.054571817e-34

#define KT  (K*T)

double calc_rj(double f);
double calc_planck(double f);
double calc_mine(double f);
int test(void); //xxx
int test1(void); //xxx

int main(int argc, char **argv)
{
    int max, i;
    double logf;
    double log_freq[1000], rj[1000], planck[1000], mine[1000];
    FILE *fp;

    test();

    // calculate the energy density vs frequency using:
    // - Rayleigh–Jeans law
    // - Planck's Law
    for (max = 0, logf = 12; logf <= 16; logf += .10) {  // xxx 14 16
        double f = pow(10, logf);

        log_freq[max] = logf;
        rj[max] = calc_rj(f);
        planck[max] = calc_planck(f);
        mine[max] = calc_mine(f);
        max++;
    }

    // print results to file plot.dat, for gnuplot
    fp = fopen("plot.dat", "w");
    for (i = 0; i < max; i++) {
        double ratio = mine[i] / rj[i];
        double ratio2 = mine[i] / planck[i];
        fprintf(fp, "%0.3f %10.3e %10.3e %10.3e  # %10f %10f\n",
                log_freq[i], 
                planck[i], 
                mine[i],
                rj[i],
                ratio, ratio2);
    }
    fclose(fp);

    // run gnuplot
#if 1
    system("gnuplot -p \
              -e \"set term x11 size 1000,600\" \
              -e \"set xrange [*:*]\" \
              -e \"set yrange [0:5e-15]\" \
              -e \"plot 'plot.dat' using 1:2 with lines linewidth 2, \
                        'plot.dat' using 1:3 with lines linewidth 2, \
                        'plot.dat' using 1:4 with lines linewidth 2\" \
                    ");
#else
#if 0
    system("gnuplot -p \
              -e \"set term x11 size 1000,600\" \
              -e \"plot 'plot.dat' using 1:2 with lines linewidth 2\" \
                    ");
    system("gnuplot -p \
              -e \"set term x11 size 1000,600\" \
              -e \"plot 'plot.dat' using 1:3 with lines linewidth 2\" \
                    ");
#else
    system("gnuplot -p \
              -e \"set term x11 size 1000,600\" \
              -e \"plot 'plot.dat' using 1:2 with lines linewidth 2, \
                        'plot.dat' using 1:3 with lines linewidth 2\" \
                    ");
#endif
#endif

    return 0;
}

double calc_rj(double f)
{
    double mode_density, energy_density;

    mode_density = (8 * M_PI * (f*f)) / (C*C*C);
    energy_density = mode_density * KT;

    return energy_density;
}

double calc_planck(double f)
{
    return ((8 * M_PI * h * (f*f*f)) / (C*C*C)) * (1 / (exp((h * f) / KT) - 1));
}

// ------------------------------

//#define MASS (4.002603 * 1.6603145E-27)  // 4He mass in kg
#define MASS (20 * 1.6603145E-27)
double array[1000000];
int max_array=0;
double maxwell_boltzmann(double velocity);
int binary_search(double element, double *array, int max);

double calc_mine(double f)
{
    double mode_density, energy_density;
    double sum_energy_quantized, avg_energy_quantized;
    int i;
    double hf = 1.22 * h * f;

    mode_density = (8 * M_PI * (f*f)) / (C*C*C);

#define MAX 5000000
#if 1
    sum_energy_quantized = 0;
    for (i = 0; i < MAX; i++) {
        double p = (double)random() / ((double)RAND_MAX + 1);
        double velocity = binary_search(p, array, max_array);
        double energy = .5 * MASS * (velocity * velocity);
        sum_energy_quantized += floor(energy / hf) * hf;
        //sum_energy_quantized += energy;
    }
    avg_energy_quantized = sum_energy_quantized / MAX;
#else
    sum_energy_quantized = 0;
    for (i = 0; i < MAX; i++) {
        double p = (double)random() / ((double)RAND_MAX + 1);
        double velocity = binary_search(p, array, max_array);
        double energy = .5 * MASS * (velocity * velocity);
        double tmp = floor(energy / hf) * hf;
        sum_energy_quantized += tmp*tmp;
    }
    avg_energy_quantized = sqrt(sum_energy_quantized / MAX);

#endif

    energy_density = mode_density * avg_energy_quantized;
    energy_density *= .66666;

    return energy_density;
}

//#undef KT
//#define KT  (K*293)


int test1(void)
{
    double array[] = {0, .1, .2, .3, .7, .8, .9, 1.0};
    char s[100];
    double value;
    int idx;

    while (true) {
        printf("? ");
        if (fgets(s, sizeof(s), stdin) == NULL) break;
        if (sscanf(s, "%lf", &value) != 1) {
            printf("INVALID\n");
            continue;
        }

        idx = binary_search(value, array, 8);

        if (idx == -1) {
            printf("NOT FOUND\n");
            continue;
        }
        if (idx < 0 || idx >= 8-1) {
            printf("ERROR bad idx %d\n", idx);
            continue;
        }
        printf("got idx=%d  value=%f  array=%f %f\n", idx, value, array[idx], array[idx+1]);
        if (value < array[idx] || value >= array[idx+1]) {
            printf("OOOPS\n");
            continue;
        }
    }
    return 0;
}

int test(void)
{
    double v;   // meter per second
    double p;
    double sum_p=0;


    array[max_array++] = 0;

    double maxv = sqrt(2 * (50*KT) / MASS);

    printf("KT = %0.3e\n", KT);
    printf(".5mv2 for v=%f = %0.3e\n", maxv, .5*MASS*3000*3000);

    for (v = 0; v < maxv; v++) {
        p = maxwell_boltzmann(v);
        //printf("%f %12.10f\n", v, p);
        sum_p += p;
        array[max_array++] = sum_p;
    }

    array[max_array++] = 1;
    printf("sum_p = %f\n", sum_p);
    printf("max array = %d\n", max_array);
    printf("max v     = %f\n", maxv);

#if 0
    printf("xxxxxxxxxxxxxxxxxxxxxxx\n");
    for (int i = 0; i < max_array; i++) {
        printf("%d - %f\n", i, array[i]);
    }
#endif

    // 1000 - 0.349904
    // 1001 - 0.350643
    int idx;
#if 0
    idx = binary_search(0.349904+.00001, array, max_array);
    printf("idx %d\n", idx);
    idx = binary_search(0.350643-.00001, array, max_array);
    printf("idx %d\n", idx);
#endif

    static int histogram[1000000];

    printf("starting\n");
    for (int i = 0; i < 10000000; i++) {
        p = (double)random() / ((double)RAND_MAX + 1);
        idx = binary_search(p, array, max_array);
        if (idx < 0 || idx >= max_array-1) {
            printf("ERROR bad idx %d\n", idx);
            exit(1);
        }
        //printf("p = %f, speed = %d\n", p, idx);

        histogram[idx]++;
    }
    printf("done\n");

    FILE *fp = fopen("test1.dat", "w");
    for (int i = 0; i < maxv; i++) {
        fprintf(fp, "%d %d\n", i, histogram[i]);
    }
    fclose(fp);

#if 0
    system("gnuplot -p \
              -e \"set term x11 size 1000,600\" \
              -e \"set xrange [*:*]\" \
              -e \"set yrange [*:*]\" \
              -e \"plot 'test1.dat' using 1:2 with lines linewidth 2\" \
                    ");
#endif

    return 0;
}

double maxwell_boltzmann(double velocity)
{
    double velocity_squared = velocity * velocity;
    double probability;

    // xxx try eliminating the mass, and just use energy
    probability = pow(MASS / (2*M_PI*KT), 1.5) * 
                 (4*M_PI) * velocity_squared *
                 exp(-MASS*(velocity_squared) / (2*KT));
    return probability;
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
