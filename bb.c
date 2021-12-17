// usage: ./bb [-t <temp_deg_k] [-z]
//   -t <temp_deg_k> : black body temperature
//   -z              : display test plot of maxwell-boltzmann distribution

#include "common.h"

//
// defines
//

#define K 1.38064852E-23   // boltzmann constant
#define C 299792458.       // speed of light
#define h 6.62607004e-34   // planck constant

#define AMU_TO_KG(amu)  ((amu) * 1.6603145E-27)
#define KG_TO_AMU(amu)  ((amu) / 1.6603145E-27)

//
// variables
//

double T = 5000;
double KT;

//
// prototypes
//

char **extra_gnuplot_cmds(void);

double calc_rj(double f);
double calc_planck(double f);
double calc_mine(double f);
void calc_mine_init(void);
void calc_mine_test(void);

// -----------------  MAIN  ----------------------------------------------------

int main(int argc, char **argv)
{
    bool test_enable = false;

    // got options
    // -t <temp_deg_k> : temperature
    // -z              : display test plot of maxwell-boltzmann distribution
    while (true) {
        int c = getopt(argc, argv, "t:z");
        if (c == -1) break;
        switch (c) {
        case 't':   // temperature (degrees K)
            if (sscanf(optarg, "%lf", &T) != 1 || T < 1 || T > 20000) {
                printf("ERROR: invalid temperature, 1 to 20000 degrees K expected\n");
                exit(1);
            }
            break;
        case 'z':
            test_enable = true;
            break;
        default:
            exit(1);
        }
    }

    // init KT, global variable
    KT = K * T;

    // print params and KT value
    printf("T  = %0.1f degrees K\n", T);
    printf("KT = %e joules\n", KT);
    printf("\n");

    // init my black body code
    calc_mine_init();
    if (test_enable) {
        calc_mine_test();
    }

    // check that both Planck and My black-body calculations agree with
    // Rayleigh-Jeans at low frequency (1 MHz); print warning if check fails
    double tst_freq   = 1000000;
    double tst_rj     = calc_rj(tst_freq);
    double tst_planck = calc_planck(tst_freq);
    double tst_mine   = calc_mine(tst_freq);
    if ((tst_rj/tst_planck < 0.995 || tst_rj/tst_planck > 1.005) ||
        (tst_rj/tst_mine < 0.995 || tst_rj/tst_mine > 1.005))
    {
        printf("WARNING at 1MHz: rj=%0.6e planck=%0.6e mine=%0.6e\n", tst_rj, tst_planck, tst_mine);
        printf("- rj/planck = %0.6f\n", tst_rj/tst_planck);
        printf("- rj/mine   = %0.6f\n", tst_rj/tst_mine);
        printf("\n");
    }

    // calculate the black-body energy density vs frequency using:
    // - Rayleigh–Jeans law
    // - Planck's Law
    // - My black body calculation, versions 1 and 2
    int max = 0;
    double logfreq;
    static double logf[10000], rj[10000], planck[10000], mine[10000];;

    printf("black-body starting\n");
    uint64_t start = microsec_timer();
    for (logfreq = 10; logfreq <= 16; logfreq += .05) {
        double f = pow(10, logfreq);
        logf[max]   = logfreq;
        rj[max]     = calc_rj(f);
        planck[max] = calc_planck(f);
        mine[max]   = calc_mine(f);
        max++;
    }
    printf("black-body complete, %0.3f secs\n", (microsec_timer() - start) / 1000000.);
    printf("\n");

    // print results to file plot.dat, for gnuplot
    printf("plotting\n");
    FILE *fp = fopen("plot.dat", "w");
    for (int i = 0; i < max; i++) {
        fprintf(fp, "%8.4f %10.3e %10.3e %10.3e # mine/rj=%0.6f mine/planck=%0.6f\n",
                logf[i], rj[i], planck[i], mine[i], 
                mine[i]/rj[i], mine[i]/planck[i]);
    }
    fclose(fp);

    // run gnuplot
    double ymax;
    char yrange[100], title[100];
    sprintf(title, "Black Body - T=%0.1f K", T);
    ymax = max_array_val(max, mine, planck, NULL);
    sprintf(yrange, "[0:%e]", 1.5*ymax);
    gnuplot(title, "plot.dat", 
            "Log Frequency", "[*:*]", 
            "Energy Density", yrange, 
            extra_gnuplot_cmds(),
            "Rayleigh–Jeans" , "1:2", "red",
            "Planck", "1:3", "purple",
            "mine", "1:4", "blue",
            NULL, NULL, NULL);
    printf("\n");

    // done
    printf("done\n");
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

#define MODE_DENSITY(wvlen) (8 * M_PI * pow((wvlen), -4))

double calc_rj(double f)
{
    return MODE_DENSITY(C/f) * KT;
}

// -----------------  PLANCK BLACK-BODY--------------------------------

double calc_planck(double f)
{
    return ((8 * M_PI * h * (f*f*f*f*f)) / (C*C*C*C)) * (1 / (exp((h * f) / KT) - 1));
}

// -----------------  MY BLACK-BODY -----------------------------------

void *hndl;
double mb_energy_probability(double energy);

double calc_mine(double f)
{
    #define MAX 1000000
    #define H_FUDGE_FACTOR .777

    double hf = (h*H_FUDGE_FACTOR) * f;
    double sum_energy_quantized = 0;
    double avg_energy_quantized, energy_density;

    for (int i = 0; i < MAX; i++) {
        double energy = probdist_get_value(hndl) * (2./3.);
        sum_energy_quantized += floor(energy / hf) * hf;
    }
    avg_energy_quantized = sum_energy_quantized / MAX;

    energy_density = MODE_DENSITY(C/f) * avg_energy_quantized;

    return energy_density;
}

void calc_mine_init(void)
{
    hndl = probdist_create(mb_energy_probability, 0, 20*KT);
}

// xxx comments
void calc_mine_test(void)
{
    printf("calc_mine_test ...\n");

    #define MAX_TEST 1000000
    double sum=0;
    for (int i = 0; i < MAX_TEST; i++) {
        sum += probdist_get_value(hndl);
    }
    double avg = sum / MAX_TEST;
    printf("- AVG   = %0.30f\n", avg);
    printf("- KT    = %0.30f\n", KT);
    printf("- RATIO = %0.6f\n", avg/KT);

    probdist_test(hndl);

    printf("\n");
}

double mb_energy_probability(double energy)
{
    return 2. / sqrt(M_PI) *
           pow(KT, -1.5) *
           sqrt(energy) *
           exp(-energy/KT);
}

// -----------------  NOT USED  ---------------------------------------

#if 0
double mb_velocity_probability(double velocity)
{
    double velocity_squared = velocity * velocity;
    double probability;

    probability = pow(mass / (2*M_PI*KT), 1.5) * 
                  (4*M_PI) * velocity_squared *
                  exp(-(mass*velocity_squared) / (2*KT));
    return probability;
}
#endif
