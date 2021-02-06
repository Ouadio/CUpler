#include "sys/time.h"
#include "time.h"
#include "utils.h"

double now()
{
    struct timeval t;
    double f_t;
    gettimeofday(&t, NULL);
    f_t = t.tv_usec;
    f_t = f_t / ((double)1E6);
    f_t += t.tv_sec;
    return (f_t);
}

double dotProd(double *x, double *y, size_t len)
{
    double result = 0;

    for (size_t i = 0; i != len; i++)
    {
        result += x[i] * y[i];
    }
    return (result);
}