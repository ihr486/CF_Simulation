#include <stdio.h>
#include <immintrin.h>
#include <random>
#include <algorithm>
#include <memory>
#include <math.h>

extern "C" void rsqrt(double *v);

int main(int argc, const char *argv[])
{
    std::random_device device;
    std::vector<uint32_t> seed_vector(10);
    std::generate(seed_vector.begin(), seed_vector.end(), std::ref(device));
    std::seed_seq seed(seed_vector.begin(), seed_vector.end());
    std::mt19937 engine(seed);
    std::uniform_real_distribution<double> dist(-4, 4);

    double *values = (double *)aligned_alloc(32, sizeof(double) * 256);

    for (int i = 0; i < 256; i++)
    {
        values[i] = pow(10.0, dist(engine));
    }

    double *ref_recip = (double *)malloc(sizeof(double) * 256);

    for (int i = 0; i < 256; i++)
    {
        ref_recip[i] = 1.0 / sqrt(values[i]);
    }

    rsqrt(values);

    double mean_error = 0.0;
    for (int i = 0; i < 256; i++)
    {
        printf("%lf = %lf\n", values[i], ref_recip[i]);
        mean_error += (values[i] - ref_recip[i]) / ref_recip[i];
    }
    mean_error /= 256;

    printf("Mean error = %lf\n", mean_error);

    return 0;
}
