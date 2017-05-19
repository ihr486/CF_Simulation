#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define THETA_DIV (64)  //Division of the coil current
#define PHI_DIV (64)    //Angular resolution of the magnetic field
#define R_DIV (32)      //Radial resolution of the magnetic field
#define L_DIV (32)      //Depth resolution of the magnetic field
#define PI (3.14159265359)  //Approximation of pi

typedef struct vector_tag
{
    double x, y, z;
} vector_t;

typedef struct coil_tag
{
    double R, P, N;
    double *sin_table, *cos_table;
} coil_t;

static void init_coil(coil_t *coil)
{
    coil->sin_table = aligned_alloc(32, sizeof(double) * THETA_DIV);
    coil->cos_table = aligned_alloc(32, sizeof(double) * THETA_DIV);
    for (int i = 0; i < THETA_DIV; i++)
    {
        double rad = 2.0 * PI * i / THETA_DIV;
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
    return (vector_t){a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}

static vector_t vscale(const vector_t v, double c)
{
    return (vector_t){v.x * c, v.y * c, v.z * c};
}

static double vsize(const vector_t v)
{
    return sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
}

static vector_t biot_savart(const vector_t v, const coil_t *coil)
{
    vector_t B = {0, 0, 0};
    double l = 2.0 * PI * coil->R / THETA_DIV;
    for (int i = 0; i < coil->N; i++)
    {
        for (int theta = 0; theta < THETA_DIV; theta++)
        {
            vector_t p = {coil->R * coil->cos_table[theta], coil->P * (i - 0.5 * (coil->N - 1)), coil->R * coil->sin_table[theta]};
            vector_t j = {-coil->sin_table[theta], 0, coil->cos_table[theta]};
            vector_t r = vsub(v, p);
            double d = vsize(r);
            B = vadd(B, vscale(vcross(j, r), 1E-7 * l / (d * d * d)));
        }
    }
    return B;
}

static double self_inductance(const coil_t *coil)
{
    double B_total = 0;
    for (int d = 0; d < L_DIV; d++)
    {
        for (int phi = 0; phi < PHI_DIV; phi++)
        {
            double rad = 2.0 * PI * phi / PHI_DIV;
            for (int r = 0; r < R_DIV; r++)
            {
                double dist = coil->R * (r + 0.5) / R_DIV;
                vector_t v = {dist * cos(rad), 0.5 * d * coil->P * (coil->N - 1) / L_DIV, dist * sin(rad)};
                double dS = PI * coil->R * coil->R / (R_DIV * R_DIV) * (2 * r + 1) / PHI_DIV;
                B_total -= biot_savart(v, coil).y * dS;
            }
        }
    }
    return B_total * coil->N / L_DIV;
}

int main(int argc, const char *argv[])
{
    coil_t tx = {30E-3, 0.22E-3, 100, NULL, NULL};
    init_coil(&tx);
    printf("Self inductance = %lf[uH]\n", self_inductance(&tx) * 1E+6);
    return 0;
}
