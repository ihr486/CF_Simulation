#include <stdio.h>
#include <math.h>

#define THETA_DIV (32)  //Division of the coil current
#define PHI_DIV (32)    //Angular resolution of the magnetic field
#define R_DIV (16)      //Radial resolution of the magnetic field
#define L_DIV (16)      //Depth resolution of the magnetic field

static const double PI = 3.14159265358;
static const double mu0 = 4E-7 * PI;
static const double rho = 1.68E-8;

struct coil_t
{
    double R;   //Coil radius [m]
    double r;   //Wire radius [m]
    double P;   //Winding pitch [m]
    int N;   //Number of turns
};

struct vector_t {
    double x, y, z;
    vector_t &operator+=(const vector_t v)
    {
        x += v.x;
        y += v.y;
        z += v.z;
    }
    double size2()
    {
        return x * x + y * y + z * z;
    }
    double size()
    {
        return sqrt(x * x + y * y + z * z);
    }
};

static vector_t operator+(const vector_t a, const vector_t b)
{
    return vector_t{a.x + b.x, a.y + b.y, a.z + b.z};
}

static vector_t operator-(const vector_t a, const vector_t b)
{
    return vector_t{a.x - b.x, a.y - b.y, a.z - b.z};
}

static vector_t operator*(const vector_t a, double c)
{
    return vector_t{a.x * c, a.y * c, a.z * c};
}

static double dot(const vector_t a, const vector_t b)
{
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

static vector_t cross(const vector_t a, const vector_t b)
{
    return vector_t{a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}

static vector_t biot_savart(const vector_t v, const coil_t& coil)
{
    vector_t B{0, 0, 0};
    double l = 2.0 * PI * coil.R / THETA_DIV;
    for (int theta = 0; theta < THETA_DIV; theta++)
    {
        double rad = 2.0 * PI * theta / THETA_DIV;
        vector_t p{coil.R * cos(rad), 0, coil.R * sin(rad)};
        vector_t j{-sin(rad), 0, cos(rad)};
        for (int i = 0; i < coil.N; i++)
        {
            double r = (v - p).size();
            p.y = coil.P * (i - 0.5 * coil.N);
            B += cross(j, v - p) * (1E-7 * l / (r * r * r));
        }
    }
    return B;
}

static double self_inductance(const coil_t& coil)
{
    double Bs = 0;
    for (int d = 0; d < L_DIV; d++)
    {
        for (int phi = 0; phi < PHI_DIV; phi++)
        {
            double rad = 2.0 * PI * phi / PHI_DIV;
            for (int r = 1; r < R_DIV; r++)
            {
                double dist = coil.R * r / R_DIV;
                vector_t v{dist * cos(rad), 0.5 * d * coil.P * coil.N / L_DIV, dist * sin(rad)};
                double dS = PI * coil.R * coil.R / (R_DIV * R_DIV) * (2 * r + 1) / PHI_DIV;
                Bs -= biot_savart(v, coil).y * dS;
            }
        }
    }
    return Bs * coil.N / L_DIV;
}

int main(int argc, const char *argv[])
{
    coil_t tx{30E-3, 0.1E-3, 0.22E-3, 100};
    double L = self_inductance(tx);
    printf("Self inductance [Biot-Savart] = %lf[uH]\n", L * 1E+6);
    //printf("Self inductance [] = %lf\n", L1);
    return 0;
}
