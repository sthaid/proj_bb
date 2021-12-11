// usage: ./bb [-t <temp_deg_k] [-m <mass_amu>] [-z]
//   -t <temp_deg_k> : black body temperature
//   -m <mass_amu>   : particle mass, this value does not affect the black body spectrum
//   -z              : display test plot of maxwell-boltzmann distribution

#include "common.h"

//
// defines
//

#define K   1.38064852E-23   // boltzmann constant
#define C   299792458.       // speed of light
#define h   6.62607004e-34   // planck constant

#define AMU_TO_KG(amu)  ((amu) * 1.6603145E-27)
#define KG_TO_AMU(amu)  ((amu) / 1.6603145E-27)

//
// variables
//

double mass           = AMU_TO_KG(4);
double T              = 298;
bool   mb_test_enable = false;
double KT;

//
// prototypes
//

char **extra_gnuplot_cmds(void);

double calc_rj(double f);
double calc_planck(double f);
double calc_mine1(double f);
double calc_mine2(double f);

void mb_init(void);
double mb_get_velocity(void);
double mb_get_energy(void);
double mb_probability(double velocity);
void mb_test(void);

// -----------------  MAIN  ----------------------------------------------------

int main(int argc, char **argv)
{
    // got options
    // -t <temp_deg_k> : temperature
    // -m <mass_amu>   : particle mass, this value does not affect the black body spectrum
    // -z              : display test plot of maxwell-boltzmann distribution
    while (true) {
        int c = getopt(argc, argv, "t:m:z");
        if (c == -1) break;
        switch (c) {
        case 't':   // temperature (degrees K)
            if (sscanf(optarg, "%lf", &T) != 1 || T < 1 || T > 20000) {
                printf("ERROR: invalid temperature, 1 to 20000 degrees K expected\n");
                exit(1);
            }
            break;
        case 'm': { // mass (AMU)
            double amu;
            if (sscanf(optarg, "%lf", &amu) != 1 || amu < 1 || amu > 1000) {
                printf("ERROR: invalid mass, 1 to 1000 AMU expected\n");
                exit(1);
            }
            mass = AMU_TO_KG(amu);
            break; }
        case 'z':
            mb_test_enable = true;
            break;
        default:
            exit(1);
        }
    }

    // init KT, global variable
    KT = K * T;

    // print params and KT value
    printf("T    = %0.1f degrees K\n", T);
    printf("mass = %0.1f AMU\n", KG_TO_AMU(mass));
    printf("KT   = %e joules\n", KT);
    printf("\n");

    // initialize
    mb_init();

    // if maxwell-boltzmann test enable then call mb_test; 
    // this test will plot the maxwell-boltzmann velocity distribution and
    // also call mb_get_velocity many times to confirm that mb_get_velocity
    // is returning random velocities in accordance to the maxwell-boltzmann
    // velocity probability distribution
    if (mb_test_enable) {
        mb_test();
        return 0;
    }

    // calculate the black-body energy density vs frequency using:
    // - Rayleigh–Jeans law
    // - Planck's Law
    // - My black body calculation, versions 1 and 2
    int max = 0;
    double logf;
    double log_freq[1000], rj[1000], planck[1000], mine1[1000], mine2[1000];;

    printf("black-body starting\n");
    uint64_t start = microsec_timer();
    for (logf = 10; logf <= 16; logf += .05) {  // .05 is good, use .5 for testing
        double f = pow(10, logf);

        log_freq[max] = logf;
        rj[max] = calc_rj(f);
        planck[max] = calc_planck(f);
        mine1[max] = calc_mine1(f);
        mine2[max] = calc_mine2(f);
        max++;
    }
    printf("black-body complete, %0.3f secs\n", (microsec_timer() - start) / 1000000.);
    printf("\n");

    // print results to file plot.dat, for gnuplot
    printf("plotting\n");
    FILE *fp = fopen("plot.dat", "w");
    for (int i = 0; i < max; i++) {
        double ratio = rj[i] / mine2[i];
        fprintf(fp, "%0.3f %10.3e %10.3e %10.3e %10.3e # %0.6f\n",
                log_freq[i], rj[i], planck[i], mine1[i], mine2[i], ratio);
    }
    fclose(fp);

    // run gnuplot
    double ymax;
    char yrange[100], title[100];
    sprintf(title, "Black Body - T=%0.1f K", T);

    ymax = max_array_val(max, mine1, planck, NULL);
    sprintf(yrange, "[0:%e]", 1.5*ymax);
    gnuplot(title, "plot.dat", 
            "Log Frequency", "[10:16]", 
            "Energy Density", yrange, 
            extra_gnuplot_cmds(),
            "Rayleigh–Jeans" , "1:2", "red",
            "Planck", "1:3", "purple",
            "Mine1", "1:4", "blue",
            NULL, NULL, NULL);

    ymax = max_array_val(max, mine2, planck, NULL);
    sprintf(yrange, "[0:%e]", 1.5*ymax);
    gnuplot(title, "plot.dat", 
            "Log Frequency", "[10:16]", 
            "Energy Density", yrange, 
            extra_gnuplot_cmds(),
            "Rayleigh–Jeans" , "1:2", "red",
            "Planck", "1:3", "purple",
            "Mine2", "1:5", "blue",
            NULL, NULL, NULL);

    // done
    return 0;
}

// return extra gnuplot cmds to display the location of the visible light 
// spectrum under the x-axis of the black-body plot;
// notes:
// - to enter utf8 using vi: ^vu<hex_code>
// - https://en.wikipedia.org/wiki/Block_Elements
// - the block element used below is U+2584 (lower half block)
char **extra_gnuplot_cmds(void)
{
    static char *extra_cmds[10];
    int max=0;

    #define ADD(wvlen_nm, color) \
        do { \
            char *p = calloc(200, sizeof(char)); \
            sprintf(p,  \
                "set label '▄' front at first %0.3f,0 center textcolor rgbcolor '%s'", \
                log10(C/((wvlen_nm)*1e-9)), \
                (color)); \
            extra_cmds[max++] = p; \
        } while (0)

    ADD(685, "red");
    ADD(605, "orange");
    ADD(580, "yellow");
    ADD(530, "green");
    ADD(475, "blue");
    ADD(420, "violet");
    extra_cmds[max] = NULL;

    return extra_cmds;
}

// -----------------  RAYLEIGH-JEANS BLACK-BODY  ----------------------

double calc_rj(double f)
{
    double mode_density, energy_density;

    mode_density = (8 * M_PI * (f*f)) / (C*C*C);
    energy_density = mode_density * KT;

    return energy_density;
}

// -----------------  PLANCK BLACK-BODY--------------------------------

double calc_planck(double f)
{
    return ((8 * M_PI * h * (f*f*f)) / (C*C*C)) * (1 / (exp((h * f) / KT) - 1));
}

// -----------------  MY BLACK-BODY VERSIONS 1 AND 2  -----------------

double calc_mine1(double f)
{
    double mode_density, energy_density, KT_quantized;;
    double hf = h * f;

    mode_density = (8 * M_PI * (f*f)) / (C*C*C);
    KT_quantized = floor(KT / hf) * hf;
    energy_density = mode_density * KT_quantized;

    return energy_density;
}

double calc_mine2(double f)
{
    #define MAX 500000
    #define SQRT_TWO_THIRDS 0.8164966

    double mode_density, energy_density;
    double hf, energy, sum_energy_quantized, avg_energy_quantized;

    mode_density = (8 * M_PI * (f*f)) / (C*C*C);

    hf = h * f;
    sum_energy_quantized = 0;
    for (int i = 0; i < MAX; i++) {
        energy = mb_get_energy() * SQRT_TWO_THIRDS;
        sum_energy_quantized += floor(energy / hf) * hf;
    }
    avg_energy_quantized = sum_energy_quantized / MAX;

    energy_density = mode_density * avg_energy_quantized * SQRT_TWO_THIRDS;

    return energy_density;
}

// -----------------  MAXWELL BOLTZMANN DISTRIBUTION  -----------------

// cum_mb: is the cumulative maxwell boltzmann probability array; indexed in delta_v units
double *cum_mb;      
int     max_cum_mb;
double  delta_v = 0.1;

void mb_init(void)
{
    double max_v = sqrt(2 * (15*KT) / mass);
    double sum_p = 0, p;

    max_cum_mb = max_v / delta_v;
    cum_mb = calloc(max_cum_mb, sizeof(double));

    for (int idx = 0; idx < max_cum_mb; idx++) {
        p = mb_probability(idx*delta_v);
        sum_p += p * delta_v;
        cum_mb[idx] = sum_p;
    }

#if 0
    printf("mb_init:\n");
    printf("  delt_v     = %f\n", delta_v);
    printf("  max_v      = %f\n", max_v);
    printf("  max_cum_mb = %d\n", max_cum_mb);
    printf("  sum_p      = %f\n", sum_p);
    //for (int idx = 0; idx < max_cum_mb; idx++) {
    //    printf("  %d  %0.18f\n", idx, cum_mb[idx]);
    //}
    printf("\n");
#endif

    assert(sum_p > 0.9999 && sum_p <= 1);
}

double mb_get_velocity(void)
{
    double rand_range_0_1;
    int idx;

    do {
        rand_range_0_1 = (double)random() / ((double)RAND_MAX + 1);
    } while (rand_range_0_1 >= cum_mb[max_cum_mb-1]);

    idx = binary_search(rand_range_0_1, cum_mb, max_cum_mb);
    if (idx < 0 || idx >= max_cum_mb-1) {
        printf("ERROR mb_get_velocity, bad idx %d\n", idx);
        exit(1);
    }

    return idx * delta_v;
}

double mb_get_energy(void)
{
    double v = mb_get_velocity();
    return 0.5 * mass * (v * v);
}

double mb_probability(double velocity)
{
    double velocity_squared = velocity * velocity;
    double probability;

    probability = pow(mass / (2*M_PI*KT), 1.5) * 
                  (4*M_PI) * velocity_squared *
                  exp(-(mass*velocity_squared) / (2*KT));
    return probability;
}

void mb_test(void)
{
    int i, velocity, max_plot_velocity, *histogram, histogram_cnt=0;
    FILE *fp;
    char title[100];
    double sum_p1=0, sum_p2=0;
    uint64_t start = microsec_timer();

    #define MAX_TEST 10000000

    printf("test starting\n");

    max_plot_velocity = sqrt(2 * (10*KT) / mass);
    histogram = calloc(max_plot_velocity, sizeof(int));
    printf("  max_plot_velocity = %d\n", max_plot_velocity);

    for (i = 0; i < MAX_TEST; i++) {
        velocity = (int)mb_get_velocity();
        if (velocity >= max_plot_velocity) {
            //printf("  velocity=%d m/s, out of range\n", velocity);
            continue;
        }
        histogram[velocity]++;
        histogram_cnt++;
    }
    printf("  MAX_TEST          = %d\n", MAX_TEST);
    printf("  histogram_cnt     = %d\n", histogram_cnt);

    fp = fopen("test.dat", "w");
    for (velocity = 0; velocity < max_plot_velocity; velocity++) {
        double p1 = (double)histogram[velocity] / histogram_cnt;
        double p2 = mb_probability(velocity);
        fprintf(fp, "%d %0.9f %0.9f\n", velocity, p1, p2);
        sum_p1 += p1;
        sum_p2 += p2;
    }
    fclose(fp);

    if (sum_p1 < 0.999 || sum_p1 > 1.001 || sum_p2 < 0.999 || sum_p2 > 1.001) {
        printf("  ERROR: sum_p1=%0.20f sum_p2=%0.6f\n", sum_p1, sum_p2);
        exit(1);
    }

    double kt_velocity = sqrt(2 * KT / mass);
    char cmd[100];
    char *extra_gnuplot_cmds[2] = {cmd, NULL};
    printf("  kt_velocity       = %f\n", kt_velocity);
    sprintf(cmd, "set label 'KT' front at first %0.3f,0 center textcolor rgbcolor 'blue'", kt_velocity);

    sprintf(title, "Maxwell Boltzmann Distribution - Mass=%0.0f AMU, T=%0.1f K", KG_TO_AMU(mass), T);
    gnuplot(title, "test.dat", 
            "Speed m/s", "[*:*]", 
            "Probability", "[*:*]", 
            extra_gnuplot_cmds,
            "", "1:2", "green", 
            "", "1:3", "blue",
            NULL, NULL, NULL);

    printf("  test done, %0.3f secs\n", (microsec_timer() - start) / 1000000.);
    printf("\n");
}
