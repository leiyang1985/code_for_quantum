#include <stdio.h> //如需使用printf()，需引入该头文件
#include <math.h>
#include <gsl/gsl_integration.h>

double f(double x, void *params)
{
    double sigma = *(double *)params;
    double f = exp(-2 * x * x / sigma / sigma);
    return f;
}

int main(void)
{
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(1000);

    double expected = sqrt(M_PI/2);//严格值

    double sigma = 1.0;

    gsl_function F;
    F.function = &f;
    F.params = &sigma;

    double result, error;
    gsl_integration_qagi(&F, 1e-3, 1e-6, 1000, w, &result, &error);

    printf("result          = % .18f\n", result);
    printf("exact result    = % .18f\n", expected);
    printf("estimated error = % .18f\n", error);
    printf("actual error    = % .18f\n", result - expected);
    printf("intervals       = %zu\n", w->size);

    gsl_integration_workspace_free(w);

    return 0;
}