#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#define BUF_MAX 1024
#define PI 3.1415926535897932
#define eps 0.000000000000000000001


// global variable
int** branch;          // branches of trellis
int** node;            // nodes of trellis
int state_len;         // state length
int number_memory;     // number of memory
int number_state;      // number of state
int series_len;        // length of series
int test_time;         // test time


// trellis diagram
typedef struct{

    int series[BUF_MAX];    // input series
    int weight;          //weight of node, hamming distance
    bool valid;          // validation of node

} Node;


// function for making ng
// input : sigma: double
// output: ng: double (gaussian noise)
double make_ng(double sigma) {

    double x1 = (double)rand() / (double)RAND_MAX;
    double x2 = (double)rand() / (double)RAND_MAX;
    double ng = sigma * sqrt(-2.0 * log(x1 + eps)) * cos(2.0 * (double)PI * x2);

    return ng;
}


// function for printing array by binary
// input : array: int[]
void print_bin_array(int array[]) {
    for(int i=0; i<series_len; ++i) printf(" %d%d", array[i*2+1], array[i*2]);
    printf("\n");
}


// function for printing array by integer
// input : array: int[]
void print_int_array(int array[]) {
    for(int i=0; i<series_len; ++i) printf("  %d", array[i]);
    printf("\n");
}


// function for probability of nonequal between two arrays
// input : corpus: int[]
//         recv: int[]
// output: err: double (probability of nonequal between two arrays)
double prob_neq_array(int corpus[], int recv[]) {

    double err = 0.0;

    for(int i=0; i<series_len; ++i) err += (corpus[i] != recv[i] ? 1.0 : 0.0);

    return (err / (double)series_len);
}


// function for convolution branched by state length
// input : pre_input: int[]
// output: res: int[] (convolutional function)
int* convolute(int pre_input[]) {

    int* input;

    input = malloc(sizeof(int) * (series_len+number_memory));
    for(int i=0; i<number_memory; ++i) input[i] = 0;
    for(int i=0; i<series_len; ++i) input[i+number_memory] = pre_input[i];
    int* res;
    res = malloc(sizeof(int) * series_len * 2);

    if(state_len == 3) {
        for(int i=0; i<series_len; ++i) {
            res[i*2+1] = input[i+2]              ^ input[i];
            res[i*2  ] = input[i+2] ^ input[i+1] ^ input[i];
        }
    }
    else if(state_len > 3) {
        for(int i=0; i<series_len; ++i) {
            res[i*2+1] = input[i+number_memory] ^ input[i+number_memory-1]              ^ input[i];
            res[i*2  ] = input[i+number_memory] ^ input[i+number_memory-1] ^ input[i+1] ^ input[i];
        }
    }
    else printf("ERROR: state_len=%d\n", state_len);

    return res;
}


// function for hamming distance between a and b
// input : a: int
//         b: int
// output: i: int (hamming distance between a and b)
int hamming_d(int a, int b) {
    int tmp = a ^ b;
    int dist;

    for(dist=0; tmp; ++dist) tmp &= tmp-1;

    return dist;
}


// function for convolution
// input : input: int (input signal)
//         current: int (current state)
// output: res: int* (array of storing result of state transision function
//              res[0]: int(current state)
//              res[1]: int(next state)
int* conv_state(int input, int current) {

    int* res;
    res = malloc(sizeof(int) * 2);

    for(int i=0; i<number_state; ++i) {
        if(i == current) {
            res[0] = (input == 0 ? branch[i][0] : branch[i][1]);
            res[1] = (input == 0 ? node[i][0] : node[i][1]);
            break;
        }
    }

    return res;
}


// function for whether node is updatable or not
// input : next_node: Node[]
//         next_in: int
//         weight: int (hamming distance)
// output: res: bool
bool updatable(Node next_node[], int next_in, int weight) {

    bool res0, res1;
    res0 = !next_node[next_in].valid;
    res1 = next_node[next_in].valid && (next_node[next_in].weight > weight);

    return (res0 || res1);
}


// function for inverse viterbi
// input : recv: int[] (received series)
// output: hatinput: int[] (estimation of input series)
int* viterbi(int recv[]) {

    // estimation of input series
    int* hatinput;
    hatinput = malloc(sizeof(int) * BUF_MAX);

    // initialize node
    Node current_node[number_state];
    Node next_node[number_state];
    for(int i=0; i<number_state; ++i) {
        current_node[i].weight = 0; current_node[i].valid = false;
        next_node[i].weight = 0; next_node[i].valid = false;
    }

    int node_len = sizeof(current_node)/sizeof(Node);

    // initialize of current_node[0]
    current_node[0].series[0] = 0;
    current_node[0].valid = true;

    int** trellis;
    trellis = malloc(sizeof(int *) * 2);
    for(int i=0; i<2; ++i) trellis[i] = malloc(sizeof(int) * 2);

    for(int tim=0; tim<series_len; ++tim) {
        for(int i=0; i<node_len; ++i) {
            if(! current_node[i].valid) continue;
            else if(series_len-tim > number_memory) {

                trellis[0] = conv_state(0, i);
                int weight = current_node[i].weight + hamming_d(recv[tim], trellis[0][0]);

                if(updatable(next_node, trellis[0][1], weight)) {
                    next_node[trellis[0][1]].weight = weight;

                    for(int j=0; j<tim; ++j) next_node[trellis[0][1]].series[j] = current_node[i].series[j];
                    next_node[trellis[0][1]].series[tim] = 0;
                }

                trellis[1] = conv_state(1, i);
                weight = current_node[i].weight + hamming_d(recv[tim], trellis[1][0]);

                if(updatable(next_node, trellis[1][1], weight)) {
                    next_node[trellis[1][1]].weight = weight;

                    for(int j=0; j<tim; ++j) next_node[trellis[1][1]].series[j] = current_node[i].series[j];
                    next_node[trellis[1][1]].series[tim] = 1;
                }

                next_node[trellis[0][1]].valid = next_node[trellis[1][1]].valid = true;
            }
            else {

                trellis[0] = conv_state(0, i);
                int weight = current_node[i].weight + hamming_d(recv[tim], trellis[0][0]);

                if(updatable(next_node, trellis[0][1], weight)) {
                    next_node[trellis[0][1]].weight = weight;

                    for(int j=0; j<tim; ++j) next_node[trellis[0][1]].series[j] = current_node[i].series[j];
                    next_node[trellis[0][1]].series[tim] = 0;
                }

                next_node[trellis[0][1]].valid = true;
            }
        }

        for(int i=0; i<node_len; ++i) {
            for(int j=0; j<BUF_MAX; ++j) current_node[i].series[j] = next_node[i].series[j];

            current_node[i].weight = next_node[i].weight;
            current_node[i].valid = next_node[i].valid;

            next_node[i].valid = false;
        }

    }

    for(int i=0; i<series_len; ++i) hatinput[i] = current_node[0].series[i];

    return hatinput;
}


// function for midmain because of calculating probability of error
// input : sigma: double (rootened variance)
// output: res: double (probability of error)
double m_fn(double sigma) {

    int* input;
    input = malloc(sizeof(int) * series_len);
    for(int i=0; i<series_len-number_memory; ++i) input[i] = rand() & 1;
    for(int i=0; i<number_memory; ++i) input[series_len-i-1] = 0;
    //printf("input     :");
    //print_int_array(input);

    // convoluted series
    int* conv;
    conv = malloc(sizeof(int) * series_len * 2);
    conv = convolute(input);
    //printf("convoluted:");
    //print_bin_array(conv);

    // convoluted series made of {-1, 1}
    double* d;
    d = malloc(sizeof(double) * series_len * 2);
    for(int i=0; i<series_len*2; ++i) d[i] = (conv[i] == 1 ? 1.0 : -1.0);

    // noised series
    double* x;
    x = malloc(sizeof(double) * series_len * 2);
    for(int i=0; i<series_len*2; ++i) x[i] = d[i] + make_ng(sigma);

    // received series
    int* recv;
    recv = malloc(sizeof(int) * series_len * 2);
    for(int i=0; i<series_len*2; ++i) recv[i] = (x[i] > 0 ? 1 : 0);
    //printf("received  :");
    //print_bin_array(recv);

    // received series by bielement
    int *post_recv;
    post_recv = malloc(sizeof(int) * series_len);
    for(int i=0; i<series_len; ++i) post_recv[i] = (recv[i*2+1] << 1) + recv[i*2];

    // estimated input
    int* hatinput;
    hatinput = viterbi(post_recv);
    //printf("hatinput  :");
    //print_int_array(hatinput);
    //printf("\n");

    double res = prob_neq_array(input, hatinput);

    return res;
}


// function for making branch of trellis
// output: branch[][]: int (data of branch of trellis)
int** make_branch() {

    int** branch;
    branch = malloc(sizeof(int *) * number_state);
    for(int i=0; i<number_state; ++i) branch[i] = malloc(sizeof(int) * 2);

    int quarter = number_state >> 2;
    int b0, b1, tmp;

    b0 = 0b00;
    b1 = 0b11;

    for(int i=0; i<quarter; ++i) {
        branch[i][0] = b0;
        branch[i][1] = b1;
        tmp = b0;
        b0 = b1;
        b1 = tmp;
    }
    b0 = 0b01;
    b1 = 0b10;
    for(int i=quarter; i<quarter*2; ++i) {
        branch[i][0] = b0;
        branch[i][1] = b1;
        tmp = b0;
        b0 = b1;
        b1 = tmp;
    }
    b0 = 0b11;
    b1 = 0b00;
    for(int i=quarter*2; i<quarter*3; ++i) {
        branch[i][0] = b0;
        branch[i][1] = b1;
        tmp = b0;
        b0 = b1;
        b1 = tmp;
    }
    b0 = 0b10;
    b1 = 0b01;
    for(int i=quarter*3; i<quarter*4; ++i) {
        branch[i][0] = b0;
        branch[i][1] = b1;
        tmp = b0;
        b0 = b1;
        b1 = tmp;
    }

    return branch;
}


// function for making node of trellis
// output: node[][]: int (data of node of trellis)
int** make_node() {

    int** node;
    node = malloc(sizeof(int *) * number_state);
    for(int i=0; i<number_state; ++i) node[i] = malloc(sizeof(int) * 2);
    for(int i=0; i<number_state; ++i) node[i/2][i%2] = node[(number_state>>1)+i/2][i%2] = i;

    return node;
}


int main(int argc, char **argv) {

    srand(time(NULL));

    printf("test_time      <- ");
    scanf("%d", &test_time);
    printf("series_len     <- ");
    scanf("%d", &series_len);
    printf("state_len      <- ");
    scanf("%d", &state_len);

    number_memory = state_len-1;
    number_state = (1 << number_memory);

    printf("number_memory  == %d\n", number_memory);
    printf("number_state   == %d\n", number_state);
    printf("\n");

    branch = malloc(sizeof(int *) * number_state);
    for(int i=0; i<number_state; ++i) branch[i] = malloc(sizeof(int) * 2);
    branch = make_branch();

    node = malloc(sizeof(int *) * number_state);
    for(int i=0; i<number_state; ++i) node[i] = malloc(sizeof(int) * 2);
    node = make_node();

    printf("db,avg\n");
    for(double db=1.0; db<=10.0; db += 0.5) {

        double sigma = 1.0;

        double s = 1.0;
        sigma = pow(10.0, -db/20.0)*sqrt(s / 2.0);

        double avg = 0.0;
        for(int i=0; i<test_time; ++i) avg += m_fn(sigma);
        //printf("db=%.6f\n", db);
        //printf("si=%.6f\n", sigma);
        //printf("pe=%.6f\n\n", avg / (double)n);
        printf("%.6f,%.6f\n", db, avg / (double)test_time);
    }

    return 0;
}
