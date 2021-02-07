#ifndef SPHERE_UTILS_H
#define SPHERE_UTILS_H

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include "math.h"
#include <time.h>

#define PI 3.14159265358979323846

//Sphere/Planet as a structure
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

//Normal vector to the directed movement plane
void getNormal(Sphere *sphere, double *n);

//Sphere initialization
void genSphere(Sphere *sphere, double r, double phi0, double theta0, double w, double gamma);

//Reference angles (alpha & beta) defining the directed movement plan
void computeRefAngles(Sphere *sp);

//Relative phi0 in the movement circle
void computeRelativePhi0(Sphere *sp);

//Movement direction
void computeMvmtDirection(Sphere *sp);

//Minimal distance between two spheres
double minimalDistance(Sphere *sp1, Sphere *sp2);

//Whether two spheres can get aligned with the center of rotation
uint8_t neverMeet(Sphere *sp1, Sphere *sp2);

//First time two spheres get the closest to each other
double tProximity(Sphere *sp1, Sphere *sp2);

//Generating a random set of spheres data (radius, velocity, direction etc..)
//moving in the SAME plan as the initial sphere defined by its spherical
//coordinates (theta0, phi0) and its motion direction (gamma)
void generatePlaneSpheresData(double theta0,
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

//Generating a random set of spheres moving in the SAME plan as the initial
//sphere defined by its spherical coordinates (theta0, phi0) and its
//motion direction (gamma)
void generatePlaneSpheres(double theta0,
                          double phi0,
                          double gamma,
                          uint64_t N,
                          double maxW,
                          double minW,
                          uint64_t maxR,
                          uint64_t minR,
                          Sphere *spheres);

#endif