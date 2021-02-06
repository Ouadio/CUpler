#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include "math.h"
#include <time.h>
#include "sphereUtils.h"
#include "utils.h"

int main(int argc, char **argv)
{

    size_t N = 5;
    if (argc > 1)
    {
        sscanf(argv[1], "%u", &N);
    }

    size_t maxR = 10000;
    size_t minR = 1;

    printf("N = %u \n", N);

    double *phis = (double *)malloc(N * sizeof(double));
    double *thetas = (double *)malloc(N * sizeof(double));
    double *uv1 = (double *)malloc(N * sizeof(double));
    double *uv2 = (double *)malloc(N * sizeof(double));
    double *uv3 = (double *)malloc(N * sizeof(double));

    double *r = (double *)malloc(N * sizeof(double));
    double *w = (double *)malloc(N * sizeof(double));

    double distance = 0;

    generatePlaneSpheres(PI / 2,       //Theta
                         PI / 6,       //Phi
                         PI / 6,       //Gamma
                         N,            //Number of spheres
                         2 * PI,       //Maximal angular speed
                         2 * PI / 100, //Minimum angular speed
                         maxR,         //Maximal radius
                         minR,         //Minimal radius
                                       //OUTPUT
                         phis,         //Phi angles
                         thetas,       //Theta angles
                         uv1,          //Movement directions x
                         uv2,          //Movement directions y
                         uv3,          //Movement directions z
                         r,            //Radiuses
                         w);           //Angular speeds

    Sphere *spheres = (Sphere *)malloc(N * sizeof(Sphere));

    for (size_t i = 0; i != N; i++)
    {
        spheres[i].phi0 = phis[i];
        spheres[i].theta0 = thetas[i];
        spheres[i].r = r[i];
        spheres[i].w = w[i];
        spheres[i].u1 = uv1[i];
        spheres[i].u2 = uv2[i];
        spheres[i].u3 = uv3[i];
    }

    free(phis);
    free(thetas);

    free(uv1);
    free(uv2);
    free(uv3);

    free(r);
    free(w);

    double start, end;

    printf("\n-----------starting computations------------ \n\n");
    start = now();

    double minDistance = maxR - minR + 1; //Rmax - Rmin +1
    size_t idx1, idx2;
    double minD;

    for (size_t i = 0; i != N; i++)
    {
        computeRefAngles(spheres + i);
        computeRelativePhi0(spheres + i);
        computeMvmtDirection(spheres + i);
    }

    for (size_t i = 0; i != (N - 1); i++)
    {
        for (size_t j = i + 1; j != (N); j++)
        {
            minD = minimalDistance(spheres + i, spheres + j);

            if (minD <= minDistance)
            {
                idx1 = i;
                idx2 = j;
                minDistance = minD;
            }
        }
    }

    Sphere spX1, spX2;
    spX1 = spheres[idx1];
    spX2 = spheres[idx2];

    double t = tProximity(&spX1, &spX2);

    end = now();
    printf("Closest planets : (%u , %u) at %lf \n", idx1, idx2, t);
    printf("Corresponding Distance = %lf \n", minDistance);

    printf("\n---------------------------------------------- \n");
    printf("Sphere %u :\n\t r = %lf, w = %lf \n\t theta0 = %lf pi, phi0 = %lf pi  \n\t uv=[%lf, %lf, %lf] \n\t Swept angle = %lf \n", idx1, spX1.r, spX1.w, spX1.theta0 / PI, spX1.phi0 / PI, spX1.u1, spX1.u2, spX1.u3, fabs(spX1.w) * t);
    printf("Sphere %u :\n\t r = %lf, w = %lf \n\t theta0 = %lf pi, phi0 = %lf pi  \n\t uv=[%lf, %lf, %lf] \n\t Swept angle = %lf \n", idx2, spX2.r, spX2.w, spX2.theta0 / PI, spX2.phi0 / PI, spX2.u1, spX2.u2, spX2.u3, fabs(spX2.w) * t);
    printf("---------------------------------------------- \n");
    printf("Elapsed time : %lf s\n", end - start);

    free(spheres);

    return 0;
}
