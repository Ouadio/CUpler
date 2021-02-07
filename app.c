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

    size_t maxR = 10000000;
    size_t minR = 1;

    printf("N = %u \n", N);

    double distance = 0;

    Sphere *spheres = (Sphere *)malloc(N * sizeof(Sphere));

    generatePlaneSpheres(PI / 2,        //Theta
                         PI / 6,        //Phi
                         PI / 6,        //Gamma
                         N,             //Number of spheres
                         2 * PI,        //Maximal angular speed
                         2 * PI / 1000, //Minimum angular speed
                         maxR,          //Maximal radius
                         minR,          //Minimal radius
                         spheres);      //Angular speeds

    double start, end;

    printf("\n-----------starting computations------------ \n\n");

    start = now();

    //Temporary distance
    double minD;
    //Minimal distance between closest spheres
    double minDistance = maxR - minR + 1; //Rmax - Rmin +1
    //Indexes of closest spheres
    size_t idx1, idx2;
    
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
