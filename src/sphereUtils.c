#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include "math.h"
#include <time.h>
#include "sphereUtils.h"

#define PI 3.14159265358979323846

void genSphere(Sphere *sphere, double r, double phi0, double theta0, double w, double gamma)
{

    sphere->r = r;
    sphere->phi0 = phi0;
    sphere->theta0 = theta0;
    sphere->w = w;

    sphere->u1 = -cos(phi0) * cos(theta0) * sin(gamma) - sin(phi0) * cos(gamma);
    sphere->u2 = -sin(phi0) * cos(theta0) * sin(gamma) + cos(phi0) * cos(gamma);
    sphere->u3 = sin(theta0) * sin(gamma);

    return;
}

void getNormal(Sphere *sphere, double *n)
{

    n[0] = -sphere->u2 * cos(sphere->theta0) + sphere->u3 * sin(sphere->theta0) * sin(sphere->phi0);
    n[1] = -sphere->u3 * cos(sphere->phi0) * sin(sphere->theta0) + sphere->u1 * cos(sphere->theta0);
    n[2] = sin(sphere->theta0) * (-sphere->u1 * sin(sphere->phi0) + sphere->u2 * cos(sphere->phi0));
    return;
}

void computeRelativePhi0(Sphere *sp)
{
    double phi0Prime = acos(sin(sp->theta0) * cos(sp->phi0 - ((sp->beta < PI) * sp->beta + (sp->beta >= PI) * (sp->beta - PI))));
    phi0Prime += (sp->beta - PI >= 0.0 && roundl(sin(sp->theta0) * 1000) / 1000 != 0 && roundf((phi0Prime - PI) * 100000) < 0.0 ? PI : 0);
    phi0Prime -= sp->beta * (roundf(fabs(sp->alpha - PI) * 1000) == 0.0 || roundf(fabs(sp->alpha) * 1000) == 0.0);

    sp->phi0relative = phi0Prime;
}

void computeRefAngles(Sphere *sp)
{
    double *n = (double *)malloc(3 * sizeof(double));
    getNormal(sp, n);

    double beta = 0;
    double alpha = acos(n[2]);

    if (fabs(n[2]) == 1.0)
    {

        if (n[2] > 0)
        {
            beta = sp->phi0;
        }
        else
        {
            beta = PI + sp->phi0;
        }
    }

    else
    {
        beta = acos(-n[1] / sqrt(pow(n[0], 2) + pow(n[1], 2)));
        if ((roundl(n[0] * 1e5)) < 0)
        {
            beta = 2 * PI - beta;
        }
    }

    sp->alpha = alpha;
    sp->beta = beta;

    free(n);
    return;
}

void computeMvmtDirection(Sphere *sp)
{
    double alpha = sp->alpha;
    double beta = sp->beta;

    int8_t direction = 1;

    if (alpha > PI / 2)
    {
        direction = -1;
    }
    if (roundl(1000 * (alpha - PI / 2)) == 0.0 && roundl(10000 * (beta - PI) >= 0))
    {
        direction = -1;
    }
    sp->w *= direction;
}

uint8_t neverMeet(Sphere *sp1, Sphere *sp2)
{

    uint8_t neverMeet = 0;

    if (roundl((sp1->w - sp2->w) * 1e10) == 0.0)
    {
        neverMeet = 1;
    }
    return (neverMeet);
}

double minimalDistance(Sphere *sp1, Sphere *sp2)
{
    //Base case : w1 and w2 are different OR w1 and w2 are the same but movements are opposite
    double distance = fabsl(sp2->r - sp1->r);

    double alpha1 = sp1->alpha;
    double alpha2 = sp2->alpha;
    double beta1 = sp1->beta;
    double beta2 = sp2->beta;

    double w1 = sp1->w;
    double w2 = sp2->w;

    //Same velocity and same direction
    if (roundl((w1 - w2) * 1e10) == 0.0)
    {
        distance = pow(sp1->r * cos(sp1->theta0) - sp2->r * cos(sp2->theta0), 2) +
                   pow(sp1->r * sin(sp1->theta0) * cos(sp1->phi0) - sp2->r * sin(sp2->theta0) * cos(sp2->phi0), 2),
        pow(sp1->r * sin(sp1->theta0) * sin(sp1->phi0) - sp2->r * sin(sp2->theta0) * sin(sp2->phi0), 2);
        distance = sqrtl(distance);
    }
    return (distance);
}

double tProximity(Sphere *sp1, Sphere *sp2)
{
    //Computing relative initial Phi0 to the movement plane
    double phi0Relat1 = sp1->phi0relative;
    double phi0Relat2 = sp2->phi0relative;

    int8_t direction1, direction2;

    double timeProximity;

    uint8_t dontMeet = neverMeet(sp1, sp2);

    if (dontMeet)
    {
        timeProximity = 0;
    }
    else
    {
        if ((phi0Relat2 - phi0Relat1) * (sp1->w - sp2->w) >= 0)
        {
            timeProximity = (phi0Relat2 - phi0Relat1) / (sp1->w - sp2->w);
        }
        else
        {
            if ((sp1->w - sp2->w) < 0)
            {
                timeProximity = (phi0Relat2 - phi0Relat1 - 2 * PI) / (sp1->w - sp2->w);
            }

            else
            {
                timeProximity = (phi0Relat2 - phi0Relat1 + 2 * PI) / (sp1->w - sp2->w);
            }
        }
    }
    return (timeProximity);
}

void generatePlaneSpheres(double theta0,
                          double phi0,
                          double gamma,
                          uint64_t N,
                          double maxW,
                          double minW,
                          uint64_t maxR,
                          uint64_t minR,
                          double *phis,
                          double *thetas,
                          double *uv1,
                          double *uv2,
                          double *uv3,
                          double *r,
                          double *w)
{
    srand(time(0));

    //Reference sphere
    Sphere sp;
    genSphere(&sp, 1, phi0, theta0, 1, gamma);

    double *n = (double *)malloc(3 * sizeof(double));

    getNormal(&sp, n);

    computeRefAngles(&sp);

    double phi0Prime = acos(sin(theta0) * cos(phi0 - ((sp.beta < PI) * sp.beta + (sp.beta >= PI) * (sp.beta - PI))));
    phi0Prime += (sp.beta - PI >= 0.0 && roundl(sin(theta0) * 1000) / 1000 != 0 && roundf((phi0Prime - PI) * 100000) < 0.0 ? PI : 0);

    double *randomAngles = (double *)malloc(N * sizeof(double));
    int8_t *randomSigns = (int8_t *)malloc(N * sizeof(int8_t));

    for (size_t i = 0; i != N; i++)
    {
        randomAngles[i] = (2 * PI) * (double)rand() / RAND_MAX;
        randomSigns[i] = 1 - 2 * (rand() & 1); //1 or -1
        r[i] = ((double)rand() / RAND_MAX) * (maxR - minR) + minR;
        w[i] = ((double)rand() / RAND_MAX) * (maxW - minW) + minW;
    }

    double *x, *y, *z;

    x = (double *)malloc(N * sizeof(double));
    y = (double *)malloc(N * sizeof(double));
    z = (double *)malloc(N * sizeof(double));

    for (size_t i = 0; i != N; i++)
    {

        //Cartesian coordinates
        x[i] = (cos(sp.beta) * cos(phi0Prime + randomAngles[i]) - sin(sp.beta) * cos(sp.alpha) * sin(phi0Prime + randomAngles[i]));
        y[i] = (sin(sp.beta) * cos(phi0Prime + randomAngles[i]) + cos(sp.beta) * cos(sp.alpha) * sin(phi0Prime + randomAngles[i]));
        z[i] = sin(sp.alpha) * sin(phi0Prime + randomAngles[i]);

        //Movement direction using a unit vector
        uv1[i] = n[1] * z[i] - n[2] * y[i];
        uv2[i] = n[2] * x[i] - n[0] * z[i];
        uv3[i] = n[0] * y[i] - n[1] * x[i];

        //Randomly flipping movement direction (remain in the same plan)
        uv1[i] *= randomSigns[i];
        uv2[i] *= randomSigns[i];
        uv3[i] *= randomSigns[i];
    }

    for (size_t i = 0; i != N; i++)
    {
        thetas[i] = acos(z[i]);
        if (roundf(thetas[i] * 100000) == 0.0 || roundf((thetas[i] - PI) * 100000) == 0.0)
        {
            phis[i] = 0.0;
        }
        else
        {
            phis[i] = y[i] < 0 ? 2 * PI - acos(x[i] / sin(thetas[i])) : acos(x[i] / sin(thetas[i]));
        }
    }

    free(x);
    free(y);
    free(z);
    free(n);

    return;
}
