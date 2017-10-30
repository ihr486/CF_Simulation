#include <stdio.h>
#include <math.h>

#define THETA_DIV (64)  //Division of the coil current
#define PHI_DIV (64)    //Angular resolution of the magnetic field
#define R_DIV (32)      //Radial resolution of the magnetic field
#define L_DIV (32)      //Depth resolution of the magnetic field

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
            p.y = coil.P * (i - 0.5 * (coil.N - 1));
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
            for (int r = 0; r < R_DIV; r++)
            {
                double dist = coil.R * (r + 0.5) / R_DIV;
                vector_t v{dist * cos(rad), 0.5 * d * coil.P * (coil.N - 1) / L_DIV, dist * sin(rad)};
                double dS = PI * coil.R * coil.R / (R_DIV * R_DIV) * (2 * r + 1) / PHI_DIV;
                Bs -= biot_savart(v, coil).y * dS;
            }
        }
    }
    return Bs * coil.N / L_DIV;
}

static double first_elliptic_integral(double k)
{
    double a = 1.0, b = sqrt(1.0 - k * k);
    while (b < 0.999 * a)
    {
        double a1 = 0.5 * (a + b);
        double b1 = sqrt(a * b);
        a = a1, b = b1;
    }
    return 0.5 * PI / a;
}

static double second_elliptic_integral_sub(double a, double b, double *I)
{
    double a1 = 0.5 * (a + b), b1 = sqrt(a * b);
    if (b1 < 0.999 * a1)
    {
        double J1 = second_elliptic_integral_sub(a1, b1, I);
        return 2.0 * J1 - a * b * (*I);
    }
    else
    {
        *I = 0.5 * PI / a1;
        return PI * (a1 - 0.5 * a * b / a1);
    }
}

static double second_elliptic_integral(double k)
{
    double I = 0.0;
    return second_elliptic_integral_sub(1.0, sqrt(1.0 - k * k), &I);
}

static double self_inductance_ref(const coil_t& coil)
{
    double a = 0.5 * coil.P * (coil.N - 1);
    double b = a / coil.R;
    double k2 = 1.0 / (1.0 + b * b);
    double k = sqrt(k2);
    double T = 4.0 / (3.0 * PI * b * k2 * k) * ((2.0 * k2 - 1) * first_elliptic_integral(k) + (1.0 - k2) * second_elliptic_integral(k) - k2 * k);
    return 4E-7 * PI * PI * coil.N * coil.N * coil.R * T / (2.0 * b);
}

static double self_inductance_ref2(const coil_t& coil)
{
    double d = coil.P * coil.N / coil.R;
    double c = 1.0 / (2.3 + 1.6 * d + 0.44 * d * d);
    return mu0 * coil.N * coil.N * coil.R * (c + log(1 + PI / d));
}

int main(int argc, const char *argv[])
{
    coil_t tx{30E-3, 0.1E-3, 0.22E-3, 100};
    double L = self_inductance(tx);
    printf("Self inductance [Biot-Savart] = %lf[uH]\n", L * 1E+6);
    double L2 = self_inductance_ref(tx);
    printf("Self inductance [1] = %lf[uH]\n", L2 * 1E+6);
    double L3 = self_inductance_ref2(tx);
    printf("Self inductance [2] = %lf[uH]\n", L3 * 1E+6);
    return 0;
}
