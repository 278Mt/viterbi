#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <stdbool.h>
#define PI 3.1415926535897932
#define eps 0.000000000000000000001


double mean(double* z, int len) {

    double sum = 0;
    for(int i=0; i<len; ++i) sum += z[i];

    return (sum / len);
}


double make_ng(double sigma) {

    double x1 = (double)rand() / (double)RAND_MAX;
    double x2 = (double)rand() / (double)RAND_MAX;
    double ng = sigma * sqrt(-2.0 * log(x1 + eps)) * cos(2.0 * (double)PI * x2);

    return ng;
}


void print_double_array(double* array, int len) {
    for(int i=0; i<len; ++i) printf(" %f", array[i]);
    printf("\n");
}


void print_int_array(int* array, int len) {
    for(int i=0; i<len; ++i) printf("%d", array[i]);
    printf("\n");
}


double prob_neq_array(int* corpus, int* recv, int len) {

    double err = 0.0;

    for(int i=0; i<len; ++i) err += (corpus[i] != recv[i] ? 1.0 : 0.0);

    return (err / (double)len);
}


double m_fn(int len, double sigma) {

    int* corpus;
    corpus = malloc(sizeof(int) * len);
    for(int i=0; i<len; ++i) corpus[i] = rand() & 1;

    double* d;
    d = malloc(sizeof(double) * len);
    for(int i=0; i<len; ++i) d[i] = (corpus[i] == 1 ? 1.0 : -1.0);

    double* x;
    x = malloc(sizeof(double) * len);
    for(int i=0; i<len; ++i) x[i] = d[i] + make_ng(sigma);

    int* recv;
    recv = malloc(sizeof(int) * len);
    for(int i=0; i<len; ++i) recv[i] = (x[i] > 0 ? 1 : 0);

    double res = prob_neq_array(corpus, recv, len);

    return res;

}

int main(int argc, char** argv) {

    srand(time(NULL));

    printf("db,avg\n");
    for(double db=1.0; db<=10.0; db += 0.5) {

        int n=10000;
        int len = 100;

        double sigma = 1.0;

        double s = 1.0;
        sigma = pow(10.0, -db/20.0)*sqrt(s / 2.0);

        double avg = 0.0;
        for(int i=0; i<n; ++i) avg += m_fn(len, sigma);
        //printf("db=%.6f\n", db);
        //printf("si=%.6f\n", sigma);
        //printf("pe=%.6f\n\n", avg / (double)n);
        printf("%.6f,%.6f\n", db, avg / (double)n);
    }

    return 0;
}
