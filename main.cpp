#include <stdio.h>
#include <math.h>

#define THETA_DIV (32)

static const double PI = 3.14159265358;
static const double mu0 = 4E-7 * PI;
static const double rho = 1.68E-8;

typedef struct
{
    double R;   //Coil radius [m]
    double r;   //Wire radius [m]
    double P;   //Winding pitch [m]
    int N;   //Number of turns
} coil_t;

typedef struct
{
    double x, y, z;
} vector_t;

static vector_t cross(const vector_t a, const vector_t b)
{
    return vector_t{a.y * b.z - a.z - b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}

static vector_t biot_savart(const vector_t v, const coil_t *coil)
{
    for (int theta = 0; theta < THETA_DIV; theta++)
    {
        double rad = 2.0 * PI * theta / THETA_DIV;
        vector_t p = {coil->R * cos(rad), 0, coil->R * sin(rad)};
        vector_t j = {-sin(rad), 0, cos(rad)};
        for (int i = 0; i < coil->N; i++)
        {
            p.y = coil->P * (i - 0.5 * coil->N);
            
        }
    }
}

int main(int argc, const char *argv[])
{
    return 0;
}
