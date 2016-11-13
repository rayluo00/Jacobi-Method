/*
 * Ray Weiming Luo
 * Benjamin Ellerby
 * CSCI 322
 *
 * Program computes 2048x2048 matrix with Jacobi iteration using multiple threads.
 * with multiple threads.
 */

#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

#include <sys/types.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <unistd.h>
#include <semaphore.h>

pthread_cond_t cond;
pthread_mutex_t mutex;
pthread_mutex_t maxmutex;
volatile int thds = 0;
volatile int broken = 0;
volatile double maxdiff = 0.0;

typedef struct {
    double **initial;
    double **result;
    unsigned num_threads;
    unsigned threadNum;
} args_t;

double max (double x, double y) {
    return (x > y) ? x : y;
}

// Clean up code to free any allocated memory used for 2d array.
void cleanup (double **initial, double **output, double **result) {
    for (int i = 0; i < 2049; ++i) {
        free(initial[i]);
        free(result[i]);
        free(output[i]);
    }
    free(initial);
    free(result);
    free(output);
}

// Computes the Jacobi iteration.
void* performMath (void *argv) {
    unsigned i, j, k, threadNum, num_threads;
    double **initial;
    double **result;

    args_t* args = (args_t*)(argv);
    initial = args->initial;
    result = args->result;
    num_threads = args->num_threads;
    threadNum = args->threadNum;

    while (1) {
        for (i = 1; i < 2048; ++i) {
            for (j = (threadNum*2048)/num_threads; j < ((threadNum+1) * 2048)/num_threads; ++j) {
                if (j > 0 && j < 2048)
                    result[i][j] = (initial[i-1][j] + initial[i+1][j] + initial[i][j-1] + initial[i][j+1]) * 0.25;
            }
        }

        pthread_mutex_lock(&mutex);
        int broken_state = broken;
        thds = thds + 1;
        if (thds == num_threads) {
            thds = 0;
            broken = !broken;
            pthread_cond_broadcast(&cond);
        }
        while (broken == broken_state) {
            pthread_cond_wait(&cond, &mutex);
        }
        pthread_mutex_unlock(&mutex);

        double maxdifflocal = 0.0;

        for (i = 1; i < 2048; ++i) {
            for (j = (threadNum*2048)/num_threads; j < ((threadNum+1) * 2048)/num_threads; ++j) {
                if (j > 0 && j < 2048)
                    maxdifflocal = max(maxdifflocal, fabs(result[i][j]-initial[i][j]));
            }
        }

        pthread_mutex_lock(&maxmutex);
        if (maxdifflocal > maxdiff) {
            maxdiff = maxdifflocal;
        }
        pthread_mutex_unlock(&maxmutex);

        pthread_mutex_lock(&mutex);
        broken_state = broken;
        thds = thds + 1;
        if (thds == num_threads) {
            thds = 0;
            broken = !broken;
            pthread_cond_broadcast(&cond);
        }
        while (broken == broken_state) {
            pthread_cond_wait(&cond, &mutex);
        }
        pthread_mutex_unlock(&mutex);
        printf("md = '%d'\n", maxdiff);

        if (maxdiff < 0.0001)
            break;

        for (i = 1; i < 2048; ++i) {
            for (j = (threadNum*2048)/num_threads; j < ((threadNum+1) * 2048)/num_threads; ++j) {
                if (j > 0 && j < 2048)
                    initial[i][j] = result[i][j];
            }
        }

        pthread_mutex_lock(&mutex);
        broken_state = broken;
        thds = thds + 1;
        if (thds == num_threads) {
            thds = 0;
            broken = !broken;
            pthread_cond_broadcast(&cond);
        }
        while (broken == broken_state) {
            pthread_cond_wait(&cond, &mutex);
        }
        pthread_mutex_unlock(&mutex);

        maxdiff = 0.0;


    }

    return NULL;
}

int main(int argc, char* argv[]) {
    int i, j;
    unsigned num_threads = 1;
    double **initial;
    double **result;
    double **output;

    initial = malloc(2049 * sizeof(double *));
    result = malloc(2049 * sizeof(double *));
    output = malloc(2049 * sizeof(double *));

    for (i = 0; i < 2049; i++) {
        result[i] = malloc(2049 * sizeof(double));
        initial[i] = malloc(2049 * sizeof(double));
        output[i] = malloc(2049 * sizeof(double));
    }

    if (argc > 1)
        num_threads = atoi(argv[1]);

    pthread_t threads[num_threads];
    args_t initArgs[num_threads];

    FILE *infp = fopen("/home/clausoa/public/input.mtx", "r");
    FILE *outfp = fopen("/home/clausoa/public/output.mtx", "r");
    char inbuf[15];
    char outbuf[15];
    for (i = 0; i < 2049; ++i) {
        for (j = 0; j < 2049; ++j) {
            fscanf(infp, "%s", inbuf);
            fscanf(outfp, "%s", outbuf);

            initial[i][j] = atof(inbuf);
            result[i][j] = initial[i][j];
            output[i][j] = atof(outbuf);
        }
    }
    fclose(infp);
    fclose(outfp);

    pthread_cond_init(&cond, NULL);
    pthread_mutex_init(&mutex, NULL);
    pthread_mutex_init(&maxmutex, NULL);

    for (i = 0; i < num_threads; ++i) {
        args_t args = {initial, result, num_threads, i};
        initArgs[i] = args;
        pthread_create(&threads[i], NULL, performMath, &initArgs[i]);
    }

    for(i = 0; i < num_threads; ++i) {
        pthread_join(threads[i], NULL);
    }

    FILE *fd;
    fd = fopen("jacobiout.txt", "w+");
    for (i = 0; i < 2048; ++i) {
        for (j = 0; j < 2048; ++j) { 
            fprintf(fd, "%.10f", result[i][j]);
            fprintf(fd, "%c", 32);
        }
    } 

    cleanup(initial, output, result);

    return EXIT_SUCCESS;
}
