// xxx label plot axis and title,  and temperature in plot
#include "common.h"

#define   K    1.380649E-23
#define   C    299792458.
#define   h    6.62607004e-34

#define MASS (4.002603 * 1.6603145E-27)  // 4He mass in kg xxx eliminate use of mass  or test with differnet mass
//#define MASS (1000 * 1.6603145E-27)  // 4He mass in kg xxx eliminate use of mass  or test with differnet mass

double T = 6000;
double KT;

double calc_rj(double f);
double calc_planck(double f);
double calc_mine(double f);
double get_energy(void);
double get_velocity(void);
void calc_mine_init(void);
void test(void); //xxx name
double maxwell_boltzmann_probability(double velocity);

// -----------------  MAIN  ----------------------------------------------------

int main(int argc, char **argv)
{
    int max, i, cnt;
    double logf, ymax;
    double log_freq[1000], rj[1000], planck[1000], mine[1000];
    char yrange[100];
    FILE *fp;

    time_t t = time(NULL);
    srandom(t);

    // get temperature from argv[1]
    if (argc > 1) {
        cnt = sscanf(argv[1], "%lf", &T);
        if (cnt != 1 || T < 1 || T > 10000) {
            printf("ERROR: T must be in range 1 to 10000 degrees K\n");
            exit(1);
        }
    }
    printf("T = %0.0f degrees K\n", T);
    KT = K * T;

    // xxx
    calc_mine_init();

    // calculate the energy density vs frequency using:
    // - Rayleigh–Jeans law
    // - Planck's Law
    for (max = 0, logf = 10; logf <= 16; logf += .10) {  // xxx .1 vs .01
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
    printf("yrange '%s'\n", yrange);
    gnuplot("plot.dat", "[*:*]", yrange, 3, 
            "1:2", "red",
            "1:3", "purple",
            "1:4", "blue");

    // done
    return 0;
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

double calc_mine(double f)
{
    #define MAX 1000000

    // xxx cleanup
    double mode_density, energy_density, energy;
    double sum_energy_quantized, avg_energy_quantized;
    double hf = sqrt(1.5) * h * f;

    mode_density = (8 * M_PI * (f*f)) / (C*C*C);

    sum_energy_quantized = 0;
    for (int i = 0; i < MAX; i++) {
        energy = get_energy();
        sum_energy_quantized += floor(energy / hf) * hf;
    }
    avg_energy_quantized = sum_energy_quantized / MAX / 1.5;

    energy_density = mode_density * avg_energy_quantized;

    return energy_density;
}

// - - - - - - - - - - - - - - - - 

double *array;  // xxx names
int    max_array;

void calc_mine_init(void)
{
    double max_velocity = sqrt(2 * (50*KT) / MASS);
    double sum_p = 0, velocity, p;

    // xxx determine max-array here
    array = calloc(max_velocity+10, sizeof(double));

    array[max_array++] = 0;
    for (velocity = 0; velocity < max_velocity; velocity++) {
        p = maxwell_boltzmann_probability(velocity);
        sum_p += p;
        array[max_array++] = sum_p;
    }
    array[max_array++] = 1;

    printf("calc_mine_init:\n");
    printf("  sum_p        = %f\n", sum_p);
    printf("  max array    = %d\n", max_array);
    printf("  max velocity = %f\n", max_velocity);

    // xxxx assert 

    test();
}

double get_energy(void)
{
    double velocity, energy;

    velocity = get_velocity();
    energy   = .5 * MASS * (velocity * velocity);
    return energy;
}

double get_velocity(void)
{
    double p;
    int velocity;

    p = (double)random() / ((double)RAND_MAX + 1);
    velocity = binary_search(p, array, max_array);

    if (velocity < 0 || velocity >= max_array-1) {
        printf("ERROR bad velocity %d\n", velocity);
        exit(1);
    }

    return velocity;
}

// - - - - - - - - - - - - - - - - 

void test(void)
{
    int *histogram = calloc(max_array, sizeof(int));
    int i, velocity, max_plot_velocity;
    FILE *fp;

    printf("test starting\n");

    for (i = 0; i < 10000000; i++) {
        velocity = get_velocity();
        histogram[velocity]++;
    }

    max_plot_velocity = sqrt(2 * (10*KT) / MASS);
    assert(max_plot_velocity <= max_array);

    fp = fopen("test.dat", "w");
    for (velocity = 0; velocity < max_plot_velocity; velocity++) {
        fprintf(fp, "%d %d\n", velocity, histogram[velocity]);
    }
    fclose(fp);

    gnuplot("test.dat", "[*:*]", "[*:*]", 1, "1:2", "green");

    printf("calc_mine_init test done\n");
}

double maxwell_boltzmann_probability(double velocity)
{
    double velocity_squared = velocity * velocity;
    double probability;

    // xxx try eliminating the mass, and just use energy
    probability = pow(MASS / (2*M_PI*KT), 1.5) * 
                 (4*M_PI) * velocity_squared *
                 exp(-MASS*(velocity_squared) / (2*KT));
    return probability;
}
