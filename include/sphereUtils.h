#ifndef SPHERE_UTILS_H
#define SPHERE_UTILS_H

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include "math.h"
#include <time.h>

#define PI 3.14159265358979323846

typedef struct Sphere
{
    double r;
    double phi0;
    double theta0;
    double w;
    double u1;
    double u2;
    double u3;

    double alpha;
    double beta;
    double phi0relative;

} Sphere;

void getAngles(double *n, double phi0, double *angles);

void genSphere(Sphere *sphere, double r, double phi0, double theta0, double w, double gamma);

void getNormal(Sphere *sphere, double *n);

void computeRefAngles(Sphere *sp);

void computeRelativePhi0(Sphere *sp);

void computeMvmtDirection(Sphere *sp);

double dotProd(double *x, double *y, size_t len);

double minimalDistance(Sphere *sp1, Sphere *sp2);

uint8_t neverMeet(Sphere *sp1, Sphere *sp2);

double tProximity(Sphere *sp1, Sphere *sp2);

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
                          double *w);

#endif