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
    float R, P;
    int N;
} coil_t;

//extern vector_t biot_savart(const vector_t v, const coil_t *coil);

static const float cos_table[PHI_DIV] = {
    1.0,
    0.995184726672,
    0.980785280403,
    0.956940335732,
    0.923879532511,
    0.881921264348,
    0.831469612303,
    0.773010453363,
    0.707106781187,
    0.634393284164,
    0.55557023302,
    0.471396736826,
    0.382683432365,
    0.290284677254,
    0.195090322016,
    0.0980171403296,
    6.12323399574e-17,
    -0.0980171403296,
    -0.195090322016,
    -0.290284677254,
    -0.382683432365,
    -0.471396736826,
    -0.55557023302,
    -0.634393284164,
    -0.707106781187,
    -0.773010453363,
    -0.831469612303,
    -0.881921264348,
    -0.923879532511,
    -0.956940335732,
    -0.980785280403,
    -0.995184726672,
    -1.0,
    -0.995184726672,
    -0.980785280403,
    -0.956940335732,
    -0.923879532511,
    -0.881921264348,
    -0.831469612303,
    -0.773010453363,
    -0.707106781187,
    -0.634393284164,
    -0.55557023302,
    -0.471396736826,
    -0.382683432365,
    -0.290284677254,
    -0.195090322016,
    -0.0980171403296,
    -1.83697019872e-16,
    0.0980171403296,
    0.195090322016,
    0.290284677254,
    0.382683432365,
    0.471396736826,
    0.55557023302,
    0.634393284164,
    0.707106781187,
    0.773010453363,
    0.831469612303,
    0.881921264348,
    0.923879532511,
    0.956940335732,
    0.980785280403,
    0.995184726672
};
static const float sin_table[PHI_DIV] = {
    0.0,
    0.0980171403296,
    0.195090322016,
    0.290284677254,
    0.382683432365,
    0.471396736826,
    0.55557023302,
    0.634393284164,
    0.707106781187,
    0.773010453363,
    0.831469612303,
    0.881921264348,
    0.923879532511,
    0.956940335732,
    0.980785280403,
    0.995184726672,
    1.0,
    0.995184726672,
    0.980785280403,
    0.956940335732,
    0.923879532511,
    0.881921264348,
    0.831469612303,
    0.773010453363,
    0.707106781187,
    0.634393284164,
    0.55557023302,
    0.471396736826,
    0.382683432365,
    0.290284677254,
    0.195090322016,
    0.0980171403296,
    1.22464679915e-16,
    -0.0980171403296,
    -0.195090322016,
    -0.290284677254,
    -0.382683432365,
    -0.471396736826,
    -0.55557023302,
    -0.634393284164,
    -0.707106781187,
    -0.773010453363,
    -0.831469612303,
    -0.881921264348,
    -0.923879532511,
    -0.956940335732,
    -0.980785280403,
    -0.995184726672,
    -1.0,
    -0.995184726672,
    -0.980785280403,
    -0.956940335732,
    -0.923879532511,
    -0.881921264348,
    -0.831469612303,
    -0.773010453363,
    -0.707106781187,
    -0.634393284164,
    -0.55557023302,
    -0.471396736826,
    -0.382683432365,
    -0.290284677254,
    -0.195090322016,
    -0.0980171403296
};

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
            vector_t p = {coil->R * cos_table[theta], coil->P * (i - 0.5 * (coil->N - 1)), coil->R * sin_table[theta]};
            vector_t j = {-sin_table[theta], 0, cos_table[theta]};
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
                vector_t v = {dist * cos_table[phi], 0.5 * d * coil->P * (coil->N - 1) / L_DIV, dist * sin_table[phi]};
                float dS = PI * coil->R * coil->R / (R_DIV * R_DIV) * (2 * r + 1) / PHI_DIV;
                vector_t B = biot_savart(v, coil);
                B_total -= dS * B.y;
            }
        }
    }
    return B_total * coil->N / L_DIV;
}

int main(int argc, const char *argv[])
{
    coil_t tx = {30E-3, 0.22E-3, 100};
    double start_time = (double)clock() / CLOCKS_PER_SEC;
    float L = self_inductance(&tx);
    double finish_time = (double)clock() / CLOCKS_PER_SEC;
    printf("Self inductance = %lf[uH]\n", L * 1E+6);
    printf("Elapsed time = %lf[ms]\n", (finish_time - start_time) * 1E+3);
    return 0;
}
