#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <stdbool.h>
#define PI 3.1415926535897932
#define eps 0.000000000000000000001


double make_ng() {

    double sigma = 1.0;

    double x1 = (double)rand() / RAND_MAX;
    double x2 = (double)rand() / RAND_MAX;
    double ng = sigma * sqrt(-2.0 * log(x1 + eps)) * cos(2.0 * (double)PI * x2);

    return ng;
}


int main(int argc, char** argv) {

    srand(time(NULL));

    int n = 10000;

    for(int i=0; i<n; ++i) {
        double ng = make_ng();
        printf("%f\n", ng);
    }

    return 0;
}
