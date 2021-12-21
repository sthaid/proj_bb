// xxx
// usage: ./bb [-rpmz] [temp_deg_k]
//   -r: plot Rayleigh-Jeans 
//   -p: plot Planck
//   -m: plot Mine
//   -z: test plot of maxwell-boltzmann energy distribution

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

double T = 5800;   // black body temperature of the sun
double KT;

//
// prototypes
//

char **extra_gnuplot_cmds(void);

double calc_rj(double wvlen);
double calc_planck(double wvlen);
double calc_mine(double wvlen);
void calc_mine_init(void);
void calc_mine_test(void);

// -----------------  MAIN  ----------------------------------------------------

int main(int argc, char **argv)
{
    bool test_enable = false;
    bool plot_rj = false;
    bool plot_planck = false;
    bool plot_mine = false;

    // got options
    while (true) {
        int c = getopt(argc, argv, "rpmz");
        if (c == -1) break;
        switch (c) {
        case 'r': plot_rj     = true; break;
        case 'p': plot_planck = true; break;
        case 'm': plot_mine   = true; break;
        case 'z': test_enable = true; break;
        default:
            exit(1);
        }
    }

    // if temperature arg provided then scan it
    if (argc > optind) {
        if (sscanf(argv[optind], "%lf", &T) != 1 || T < 1 || T > 20000) {
            printf("ERROR: invalid temperature, 1 to 20000 degrees K expected\n");
            exit(1);
        }
    }

    // if no plots selected then default to all of them
    if (plot_rj == false && plot_planck == false && plot_mine == false) {
        plot_rj     = true;
        plot_planck = true;
        plot_mine   = true;
    }

    // must at a minimum choose to plot either mine or planck
    if (plot_planck == false && plot_mine == false) {
        printf("ERROR: must at a minimum choose to plot either mine or planck\n");
        exit(1);
    }

    // print opts and args
    printf("T           = %0.1f degrees K\n", T);
    printf("plot_rj     = %s\n", bool2str(plot_rj));
    printf("plot_planck = %s\n", bool2str(plot_planck));
    printf("plot_mine   = %s", bool2str(plot_mine));
    if (test_enable) {
        printf("   test_enable");
    }
    printf("\n");
    printf("\n");

    // init KT, global variable
    KT = K * T;

    // init my black body code
    calc_mine_init();
    if (test_enable) {
        calc_mine_test();
    }

    // check that both Planck and My black-body calculations agree with
    // Rayleigh-Jeans at low frequency (1 MHz); print warning if check fails
    double tst_freq   = 1000000;
    double tst_rj     = calc_rj(C/tst_freq);
    double tst_planck = calc_planck(C/tst_freq);
    double tst_mine   = calc_mine(C/tst_freq);
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
    #define START 10
    #define END   16
    #define INCR .05 
    uint64_t start = microsec_timer();
    for (logfreq = START; logfreq <= END; logfreq += INCR) {
        double f = pow(10, logfreq);
        double wvlen = C/f;
        logf[max] = logfreq;
        if (plot_rj)     rj[max]     = calc_rj(wvlen);
        if (plot_planck) planck[max] = calc_planck(wvlen);
        if (plot_mine)   mine[max]   = calc_mine(wvlen);
        max++;
    }
    printf("black-body complete, %0.3f secs\n", (microsec_timer() - start) / 1000000.);
    printf("\n");

    // if plot_planck is enabled then scan the planck spectrum to determine
    // the frequency range, and print it
    if (plot_planck) {
        double max_val = max_array_val(max, planck, NULL);
        double min_freq = -1, max_freq = -1;
        for (int i = 0; i < max; i++) {
            if (min_freq == -1 && planck[i] > max_val*.10) {
                min_freq = pow(10,logf[i]);
            }
            if (max_freq == -1 && min_freq != -1 && planck[i] < max_val*.10) {
                max_freq = pow(10,logf[i]);
                break;
            }
        }
        printf("PLANCK SPECTURM: T = %0.1lf   "
               "FREQ_RANGE_THZ = %0.1lf to %0.1lf   "
               "WAVELEN_RANGE_NM = %0.0lf to %0.0lf\n",
              T,
              min_freq/1e12, max_freq/1e12,
              C/max_freq*1e9, C/min_freq*1e9);
        printf("\n");
    }

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
            3,
            !plot_rj     ? NULL : "Rayleigh–Jeans" , "1:2", "red",
            !plot_planck ? NULL : "Planck", "1:3", "purple",
            !plot_mine   ? NULL : "Mine", "1:4", "blue");
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

// the "(8 * M_PI / pow((wvlen),4)" term is from here, and allowing for two polarizations
//   http://www.reading.ac.uk/physicsnet/units/3/3pha4/Lectures/l1.pdf
// the "(C / (4 * M_PI))" is for units conversion
#define MODE_DENSITY(wvlen) ((8 * M_PI / pow((wvlen),4)) * (C / (4 * M_PI)))

double calc_rj(double wvlen)
{
    return MODE_DENSITY(wvlen) * KT;
}

// -----------------  PLANCK BLACK-BODY--------------------------------

double calc_planck(double wvlen)
{
    return ((2 * h * C * C) / pow(wvlen, 5)) * (1 / (exp((h * C) / (wvlen * KT)) - 1));
}

// -----------------  MY BLACK-BODY -----------------------------------

void *hndl;
double mb_energy_probability(double energy);

double calc_mine(double wvlen)
{
    #define MAX 1000000
    #define H_FUDGE_FACTOR .777

    double hf = (h*H_FUDGE_FACTOR) * (C/wvlen);
    double sum_energy_quantized = 0;
    double avg_energy_quantized;

    for (int i = 0; i < MAX; i++) {
        double energy = probdist_get_value(hndl) * (2./3.);
        sum_energy_quantized += floor(energy / hf) * hf;
    }
    avg_energy_quantized = sum_energy_quantized / MAX;

    return MODE_DENSITY(wvlen) * avg_energy_quantized;
}

void calc_mine_init(void)
{
    hndl = probdist_create(mb_energy_probability, 0, 20*KT);
}

void calc_mine_test(void)
{
    printf("calc_mine_test ...\n");

    // the probdist_get_value returns a random energy value that conforms to the
    // maxwell-boltzmann energy probability distribution; averaging these values should
    // result in 3/2 KT
    #define MAX_TEST 1000000
    double sum=0;
    for (int i = 0; i < MAX_TEST; i++) {
        sum += probdist_get_value(hndl);
    }
    double avg_energy = sum / MAX_TEST;

    // print result, so user can confirm the avg_energy is 3/2 KT
    printf("- AVG ENERGY = %0.30f\n", avg_energy);
    printf("- 3/2 KT     = %0.30f\n", 1.5 * KT);
    printf("- RATIO      = %0.6f\n", avg_energy/(1.5*KT));
    printf("\n");

    // probdist_test will plot both the distribution and the random values
    // (those returned by probdist_get_value), which should conform to the distribtuion
    probdist_test(hndl);
}

double mb_energy_probability(double energy)
{
    return 2. / sqrt(M_PI) *
           pow(KT, -1.5) *
           sqrt(energy) *
           exp(-energy/KT);
}

#if 0   // no longer used
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
