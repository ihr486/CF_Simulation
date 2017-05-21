#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define THETA_DIV (64)  //Division of the coil current
#define PHI_DIV (64)    //Angular resolution of the magnetic field
#define R_DIV (32)      //Radial resolution of the magnetic field
#define L_DIV (32)      //Depth resolution of the magnetic field
#define PI (3.14159265359)  //Approximation of pi

typedef struct vector_tag
{
    float x, y, z;
} vector_t;

typedef struct coil_tag
{
    float R, P, N;
    float *sin_table, *cos_table;
} coil_t;

static void init_coil(coil_t *coil)
{
    coil->sin_table = aligned_alloc(32, sizeof(float) * THETA_DIV);
    coil->cos_table = aligned_alloc(32, sizeof(float) * THETA_DIV);
    for (int i = 0; i < THETA_DIV; i++)
    {
        float rad = 2.0 * PI * i / THETA_DIV;
        coil->sin_table[i] = sin(rad);
        coil->cos_table[i] = cos(rad);
    }
}

//extern vector_t biot_savart(const vector_t v, const coil_t *coil);

static vector_t vadd(const vector_t a, const vector_t b)
{
    return (vector_t){a.x + b.x, a.y + b.y, a.z + b.z};
}

static vector_t vsub(const vector_t a, const vector_t b)
{
    return (vector_t){a.x - b.x, a.y - b.y, a.z - b.z};
}

static vector_t vcross(const vector_t a, const vector_t b)
{
    //return (vector_t){a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
    return (vector_t){a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y};
}

static vector_t vscale(const vector_t v, float c)
{
    return (vector_t){v.x * c, v.y * c, v.z * c};
}

static float vsize(const vector_t v)
{
    return sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
}

static vector_t biot_savart(const vector_t v, const coil_t *coil)
{
    vector_t B = {0, 0, 0};
    float l = 2.0 * PI * coil->R / THETA_DIV;
    for (int i = 0; i < coil->N; i++)
    {
        for (int theta = 0; theta < THETA_DIV; theta++)
        {
            vector_t p = {coil->R * coil->cos_table[theta], coil->P * (i - 0.5 * (coil->N - 1)), coil->R * coil->sin_table[theta]};
            vector_t j = {-coil->sin_table[theta], 0, coil->cos_table[theta]};
            vector_t r = vsub(v, p);
            float d = vsize(r);
            B = vadd(B, vscale(vcross(j, r), 1E-7 * l / (d * d * d)));
        }
    }
    return B;
}

static float self_inductance(const coil_t *coil)
{
    float B_total = 0;
    for (int d = 0; d < L_DIV; d++)
    {
        for (int phi = 0; phi < PHI_DIV; phi++)
        {
            for (int r = 0; r < R_DIV; r++)
            {
                float dist = coil->R * (r + 0.5) / R_DIV;
                vector_t v = {dist * coil->cos_table[phi], 0.5 * d * coil->P * (coil->N - 1) / L_DIV, dist * coil->sin_table[phi]};
                float dS = PI * coil->R * coil->R / (R_DIV * R_DIV) * (2 * r + 1) / PHI_DIV;
                vector_t B = biot_savart(v, coil);
                B_total -= dS * (B.y + B.x * 1E-4 + B.z * 1E-4);
            }
        }
    }
    return B_total * coil->N / L_DIV;
}

int main(int argc, const char *argv[])
{
    coil_t tx = {30E-3, 0.22E-3, 100, NULL, NULL};
    init_coil(&tx);
    double start_time = (double)clock() / CLOCKS_PER_SEC;
    float L = self_inductance(&tx);
    double finish_time = (double)clock() / CLOCKS_PER_SEC;
    printf("Self inductance = %lf[uH]\n", L * 1E+6);
    printf("Elapsed time = %lf[ms]\n", (finish_time - start_time) * 1E+3);
    return 0;
}
