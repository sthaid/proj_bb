// show where light specturm is in plot
// labelling gnuplot 
// args for test plot and for Mass

// xxx optimize
// xxx label plot axis and title,  and temperature in plot
#include "common.h"

#define   K    1.380649E-23
#define   C    299792458.
#define   h    6.62607004e-34

#define AMU_TO_KG(amu)  ((amu) * 1.6603145E-27)
#define KG_TO_AMU(amu)  ((amu) / 1.6603145E-27)

double mass        = AMU_TO_KG(4);
double T           = 298;
bool   test_enable = false;
double KT;

char **extra_cmds(void);
double calc_rj(double f);
double calc_planck(double f);
double calc_mine(double f);
double get_energy(void);
double get_velocity(void);
void init(void);
void test(void); //xxx name
double maxwell_boltzmann_probability(double velocity);

// -----------------  MAIN  ----------------------------------------------------

int main(int argc, char **argv)
{
    int max, i;
    double logf, ymax;
    double log_freq[1000], rj[1000], planck[1000], mine[1000];
    char yrange[100];
    FILE *fp;

    // xxx
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
            test_enable = true;
            break;
        default:
            exit(1);
        }
    }

    printf("T    = %0.1f degrees K\n", T);
    printf("mass = %0.1f AMU\n", KG_TO_AMU(mass));
    KT = K * T;  // xxx print this ?

    // xxx
    init();
    if (test_enable) {
        test();
    }

    // calculate the energy density vs frequency using:
    // - Rayleigh–Jeans law
    // - Planck's Law
    // - xxx
    printf("starting\n");
    for (max = 0, logf = 10; logf <= 16; logf += .050) {  // xxx .050 is good  USE .5 for testing
        double f = pow(10, logf);

        log_freq[max] = logf;
        rj[max] = calc_rj(f);
        planck[max] = calc_planck(f);
        mine[max] = calc_mine(f);
        max++;
    }
    // xxx print how long it took

    // print results to file plot.dat, for gnuplot
    printf("plotting\n");
    fp = fopen("plot.dat", "w");
    for (i = 0; i < max; i++) {
        double ratio = mine[i] / rj[i];  // xxx temp
        double ratio2 = mine[i] / planck[i];  // xxx temp
        fprintf(fp, "%0.3f %10.3e %10.3e %10.3e  # %10f %10f\n",
                log_freq[i], rj[i], planck[i], mine[i], ratio, ratio2);
    }
    fclose(fp);

    // find largest value in planck[] or mine[] array, which will be used for gnuplot yrange
    ymax = 0;
    for (i = 0; i < max; i++) {
        if (mine[i] > ymax) ymax = mine[i];
        if (planck[i] > ymax) ymax = planck[i];
    }

    // run gnuplot
    sprintf(yrange, "[0:%e]", 1.5*ymax);
    gnuplot("", "plot.dat", 
            "Frequency", "[*:*]", 
            "Energy Density", yrange, 
            extra_cmds(),
            "Rayleigh–Jeans" , "1:2", "red",
            "Planck", "1:3", "purple",
            "Mine", "1:4", "blue",
            NULL, NULL, NULL);

    // done
    return 0;
}

char **extra_cmds(void)
{
    static char *extra_cmds[10];
    int max=0;

    #define ADD(wvlen_nm, color) \
        do { \
            char *p = calloc(200, sizeof(char)); \
            sprintf(p,  \
                "set label '▄' front at graph %0.3f,0 center textcolor rgbcolor '%s'", \
                (log10(C/((wvlen_nm)*1e-9)) - 10) / 6,  /* xxx comment */ \
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

// -----------------  RAYLEIGH-JEANS  ---------------------------------

double calc_rj(double f)
{
    double mode_density, energy_density;

    mode_density = (8 * M_PI * (f*f)) / (C*C*C);
    energy_density = mode_density * KT;

    return energy_density;
}

// -----------------  PLANCK  -----------------------------------------

double calc_planck(double f)
{
    return ((8 * M_PI * h * (f*f*f)) / (C*C*C)) * (1 / (exp((h * f) / KT) - 1));
}

// -----------------  MINE  -------------------------------------------

// xxx move prototypes here

double calc_mine(double f)
{
    //#define MAX 100000
    #define MAX 500000

    // xxx cleanup
    double mode_density, energy_density, energy;
    double sum_energy_quantized, avg_energy_quantized;
    double hf = h * f * sqrt(2./3.);

    mode_density = (8 * M_PI * (f*f)) / (C*C*C);

    sum_energy_quantized = 0;
    for (int i = 0; i < MAX; i++) {
        energy = get_energy() * (2./3.);
        sum_energy_quantized += floor(energy / hf) * hf;
    }
    avg_energy_quantized = sum_energy_quantized / MAX;

    energy_density = mode_density * avg_energy_quantized;

    return energy_density;
}

// - - - - - - - - - - - - - - - - 

double *array;  // xxx names
int    max_array;
double delta_v;

void init(void)
{
    double max_velocity = sqrt(2 * (15*KT) / mass);
    double sum_p = 0, p;
    int idx;

    delta_v = 0.1;
    max_array = max_velocity / delta_v;

    array = calloc(max_array, sizeof(double));

    for (idx = 0; idx < max_array; idx++) {
        p = maxwell_boltzmann_probability(idx*delta_v);
        sum_p += p * delta_v;
        array[idx] = sum_p;
    }

    printf("init:\n");
    printf("  sum_p        = %f\n", sum_p);
    printf("  max_array    = %d\n", max_array);
    printf("  max_velocity = %f\n", max_velocity);
    printf("  delt_v       = %f\n", delta_v);

#if 0
    for (idx = 0; idx < max_array; idx++) {
        printf("%d  %0.18f\n", idx, array[idx]);
    }
    exit(1);
#endif
}

// xxx get rid of this, and put all these in maxwel boltzman section
double get_energy(void)
{
    double velocity, energy;

    velocity = get_velocity();
    energy   = (.5 * mass) * (velocity * velocity);
    return energy;
}

double get_velocity(void)
{
    double p;
    int idx;

again:
    p = (double)random() / ((double)RAND_MAX + 1);
    if (p >= array[max_array-1]) goto again;

    // xxx can this be optimized
    idx = binary_search(p, array, max_array);

    if (idx < 0 || idx >= max_array-1) {
        printf("ERROR bad idx %d\n", idx);
        exit(1);
    }

    return idx * delta_v;
}

double maxwell_boltzmann_probability(double velocity)
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
void test(void)
{
    int *histogram = calloc(max_array, sizeof(int));  // xxx calloc size
    int i, velocity, max_plot_velocity;
    FILE *fp;
    char title[100];
    double sum_p1=0, sum_p2=0;
    double max = -INFINITY;

    #define MAX_TEST 10000000

    printf("test starting\n");

    for (i = 0; i < MAX_TEST; i++) {
        velocity = (int)get_velocity();
        histogram[velocity]++;
    }

    max_plot_velocity = sqrt(2 * (10*KT) / mass);
    printf("max_plot_velocity = %d\n", max_plot_velocity);
    assert(max_plot_velocity <= max_array);  // xxx check this

    fp = fopen("test.dat", "w");
    for (velocity = 0; velocity < max_plot_velocity; velocity++) {
        double p1 = (double)histogram[velocity] / MAX_TEST;
        double p2 = maxwell_boltzmann_probability(velocity);
        fprintf(fp, "%d %0.9f %0.9f\n", velocity, p1, p2);
        sum_p1 += p1;
        sum_p2 += p2;
        if (p1 > max) max = p1;
        if (p2 > max) max = p2;
    }
    fclose(fp);
    printf("SUM_P1/2 %0.6f %0.6f\n", sum_p1, sum_p2);

    sprintf(title, "Maxwell Boltzmann - m=%0.0f AMU, t=%0.1f K", KG_TO_AMU(mass), T);
    gnuplot(title, "test.dat", 
            "Speed m/s", "[*:*]", 
            "Probability", "[*:*]", 
            NULL,
            "", "1:2", "green", 
            "", "1:3", "blue",
            NULL, NULL, NULL);

    printf("init test done\n");
}
