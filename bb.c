// xxx
// mark KT on MB
// xxx optimize

#include "common.h"

//
// defines
//

#define   K    1.380649E-23
#define   C    299792458.
#define   h    6.62607004e-34

#define AMU_TO_KG(amu)  ((amu) * 1.6603145E-27)
#define KG_TO_AMU(amu)  ((amu) / 1.6603145E-27)

//
// variables
//

double mass        = AMU_TO_KG(4);
double T           = 298;
bool   mb_test_enable = false;
double KT;

//
// prototypes
//

char **extra_gnuplot_cmds(void);

double calc_rj(double f);
double calc_planck(double f);
double calc_mine(double f);

double mb_get_energy(void);
double mb_get_velocity(void);
void mb_init(void);
void mb_test(void);
double mb_probability(double velocity);

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

    // initialize
    mb_init();

    if (mb_test_enable) {
        mb_test();
    }

    // calculate the black-body energy density vs frequency using:
    // - Rayleigh–Jeans law
    // - Planck's Law
    // - My code
    int max = 0;
    double logf;
    double log_freq[1000], rj[1000], planck[1000], mine[1000];

    printf("black-body starting\n");
    uint64_t start = microsec_timer();
    for (logf = 10; logf <= 16; logf += .050) {  // .050 is good, use .5 for testing
        double f = pow(10, logf);

        log_freq[max] = logf;
        rj[max] = calc_rj(f);
        planck[max] = calc_planck(f);
        mine[max] = calc_mine(f);
        max++;
    }
    printf("black-body complete, %0.3f secs\n", (microsec_timer() - start) / 1000000.);

    // print results to file plot.dat, for gnuplot
    printf("plotting\n");
    FILE *fp = fopen("plot.dat", "w");
    for (int i = 0; i < max; i++) {
        double ratio_rj = mine[i] / rj[i];
        double ratio_planck = mine[i] / planck[i];
        fprintf(fp, "%0.3f %10.3e %10.3e %10.3e  # %10f %10f\n",
                log_freq[i], rj[i], planck[i], mine[i], ratio_rj, ratio_planck);
    }
    fclose(fp);

    // run gnuplot
    double ymax;
    char yrange[100];
    ymax = max_array_val(max, mine, planck, NULL);
    sprintf(yrange, "[0:%e]", 1.5*ymax);
    gnuplot("", "plot.dat", 
            "Frequency", "[*:*]", 
            "Energy Density", yrange, 
            extra_gnuplot_cmds(),
            "Rayleigh–Jeans" , "1:2", "red",
            "Planck", "1:3", "purple",
            "Mine", "1:4", "blue",
            NULL, NULL, NULL);

    // done
    return 0;
}

// return extra gnuplot cmds to display the location of the visible light 
// spectrum under the x-axis of the black-body plot;
// notes:
// - the '10' and '6' values in the code are the x-axis origin and span
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
                "set label '▄' front at graph %0.3f,0 center textcolor rgbcolor '%s'", \
                (log10(C/((wvlen_nm)*1e-9)) - 10) / 6, \
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

// -----------------  MY BLACK-BODY  ----------------------------------

// xxx need to understand the 2/3
// xxx and comments
double calc_mine(double f)
{
    #define MAX 500000

    double mode_density, energy_density, energy;
    double sum_energy_quantized, avg_energy_quantized;
    double hf = h * f * sqrt(2./3.);

    mode_density = (8 * M_PI * (f*f)) / (C*C*C);

    sum_energy_quantized = 0;
    for (int i = 0; i < MAX; i++) {
        energy = mb_get_energy() * (2./3.);
        sum_energy_quantized += floor(energy / hf) * hf;
    }
    avg_energy_quantized = sum_energy_quantized / MAX;

    energy_density = mode_density * avg_energy_quantized;

    return energy_density;
}

// -----------------  MAXWELL BOLTZMANN  ------------------------------

// cum_mb_p: is the cumulative maxwell boltzmann probability array;
//           indexed in delta_v units
double *cum_mb_p;      
int     max_cum_mb_p;
double  delta_v = 0.1;

void mb_init(void)
{
    double max_v = sqrt(2 * (15*KT) / mass);
    double sum_p = 0, p;

    max_cum_mb_p = max_v / delta_v;
    cum_mb_p = calloc(max_cum_mb_p, sizeof(double));

    for (int idx = 0; idx < max_cum_mb_p; idx++) {
        p = mb_probability(idx*delta_v);
        sum_p += p * delta_v;
        cum_mb_p[idx] = sum_p;
    }

    printf("mb_init:\n");
    printf("  delt_v       = %f\n", delta_v);
    printf("  max_v        = %f\n", max_v);
    printf("  sum_p        = %f\n", sum_p);  // xxx assert
    printf("  max_cum_mb_p = %d\n", max_cum_mb_p);

#if 0
    for (int idx = 0; idx < max_cum_mb_p; idx++) {
        printf("%d  %0.18f\n", idx, cum_mb_p[idx]);
    }
#endif
}

// xxx elim this ?
double mb_get_energy(void)
{
    double velocity, energy;

    velocity = mb_get_velocity();
    energy   = (.5 * mass) * (velocity * velocity);
    return energy;
}

double mb_get_velocity(void)
{
    double p;
    int idx;

    do {
        p = (double)random() / ((double)RAND_MAX + 1);
    } while (p >= cum_mb_p[max_cum_mb_p-1]);

    idx = binary_search(p, cum_mb_p, max_cum_mb_p);
    if (idx < 0 || idx >= max_cum_mb_p-1) {
        printf("ERROR bad idx %d\n", idx);
        exit(1);
    }

    return idx * delta_v;
}

double mb_probability(double velocity)
{
    double velocity_squared = velocity * velocity;
    double probability;

    probability = pow(mass / (2*M_PI*KT), 1.5) * 
                  (4*M_PI) * velocity_squared *
                  exp(-mass*(velocity_squared) / (2*KT));
    return probability;
}

// - - - - - - - - - - - - - - - - 

// xxx also print the maxwell boltzmen probability
void mb_test(void)
{
    int *histogram = calloc(max_cum_mb_p, sizeof(int));  // xxx calloc size
    int i, velocity, max_plot_velocity;
    FILE *fp;
    char title[100];
    double sum_p1=0, sum_p2=0;
//    double max = -INFINITY;

    #define MAX_TEST 10000000

    printf("test starting\n");

    for (i = 0; i < MAX_TEST; i++) {
        velocity = (int)mb_get_velocity();
        histogram[velocity]++;
    }

    max_plot_velocity = sqrt(2 * (10*KT) / mass);
    printf("max_plot_velocity = %d\n", max_plot_velocity);
    assert(max_plot_velocity <= max_cum_mb_p);  // xxx check this

    fp = fopen("test.dat", "w");
    for (velocity = 0; velocity < max_plot_velocity; velocity++) {
        double p1 = (double)histogram[velocity] / MAX_TEST;
        double p2 = mb_probability(velocity);
        fprintf(fp, "%d %0.9f %0.9f\n", velocity, p1, p2);
        sum_p1 += p1;
        sum_p2 += p2;
        //xxx if (p1 > max) max = p1;
        //xxx if (p2 > max) max = p2;
    }
    fclose(fp);
    printf("SUM_P1/2 %0.6f %0.6f\n", sum_p1, sum_p2);
    // asserts

    double ktv = sqrt(2 * KT / mass);
    char cmd[100];
    char *extra_gnuplot_cmds[2] = {cmd, NULL};
    printf("ktv %f\n", ktv);
    sprintf(cmd, "set label 'KT' front at graph %0.3f,0 center textcolor rgbcolor 'blue'", 
          ktv / 4000);  //xxx 4000
            

    sprintf(title, "Maxwell Boltzmann - m=%0.0f AMU, t=%0.1f K", KG_TO_AMU(mass), T);
    gnuplot(title, "test.dat", 
            "Speed m/s", "[*:*]", 
            "Probability", "[*:*]", 
            extra_gnuplot_cmds,
            "", "1:2", "green", 
            "", "1:3", "blue",
            NULL, NULL, NULL);

    printf("init test done\n");
}
