// gcc -Wall -O2 -o t1 t1.c -lm; t1

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>

#define N 1000000

void random_vector(double magnitude, double * x, double * y, double * z);
double hypotenuse(double x, double y, double z);

// -----------------  MAIN  -------------------------------------

int main(int argc, char **argv)
{
    double x,y,z;
    double sum_x=0, sum_y=0, sum_z=0;
    double sum_x2=0, sum_y2=0, sum_z2=0;

    for (int i = 0; i < N; i++) {
        // generate vector with random direction and length of 1
        random_vector(1, &x, &y, &z);

        // sum the components
        sum_x += fabs(x);
        sum_y += fabs(y);
        sum_z += fabs(z);

        sum_x2 += x*x;
        sum_y2 += y*y;
        sum_z2 += z*z;
    }

    // print average x component
    printf("avg_x = %0.6f\n", sum_x/N);
    printf("avg_y = %0.6f\n", sum_y/N);
    printf("avg_z = %0.6f\n", sum_z/N);
    printf("\n");

    printf("avg_x2 = %0.6f\n", sum_x2/N);
    printf("avg_y2 = %0.6f\n", sum_y2/N);
    printf("avg_z2 = %0.6f\n", sum_z2/N);

    return 0;
}

// -----------------  SUPPORT------------------------------------

// returns a vector whose length equals 'magnitude' and with a random direction
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
    
#if 0
    // verification
    double magnitude_check = hypotenuse(*x, *y, *z);
    if (fabsf(magnitude_check-magnitude) > 2) {
        FATAL("magnitude=%f magnitude_check=%f, xyz=%f %f %f\n",
              magnitude, magnitude_check, *x, *y, *z);
    }
#endif
}

double hypotenuse(double x, double y, double z)
{
    return sqrt(x*x + y*y + z*z);
}

