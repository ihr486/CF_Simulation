#include <stdio.h>
#include <math.h>
#include <time.h>
#include <random>

#define TEST_SIZE (1024*1024*64)

extern "C" void dotproduct(const double *A, const double *B, double *C, int num);

int main(int argc, const char *argv[])
{
    double *A = (double *)aligned_alloc(32, sizeof(double) * TEST_SIZE * 4);
    double *B = (double *)aligned_alloc(32, sizeof(double) * TEST_SIZE * 4);
    double *C = (double *)aligned_alloc(32, sizeof(double) * TEST_SIZE);

    std::mt19937 engine(time(NULL));
    std::uniform_real_distribution<double> dist(0, 1);

    for (int i = 0; i < TEST_SIZE * 4; i++)
    {
        A[i] = dist(engine);
        B[i] = dist(engine);
    }

    double start_time = (double)clock() / CLOCKS_PER_SEC;

    dotproduct(A, B, C, TEST_SIZE);

    double finish_time = (double)clock() / CLOCKS_PER_SEC;

    double mean_error = 0;
    for (int i = 0; i < TEST_SIZE; i++)
    {
        double ref_dp = 0;
        for (int j = 0; j < 4; j++)
        {
            ref_dp += A[i * 4 + j] * B[i * 4 + j];
        }
        mean_error += (C[i] - ref_dp) / ref_dp;
    }
    mean_error /= TEST_SIZE;

    printf("Mean error = %le\n", mean_error);
    printf("Elapsed time = %lf[ms]\n", (finish_time - start_time) * 1E+3);
    return 0;
}
