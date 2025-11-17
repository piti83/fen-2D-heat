#ifndef SRC_INTEGRALS_C_
#define SRC_INTEGRALS_C_

#include <mesh.h>

double Gauss1D2P(double (*f)(double));
double Gauss1D3P(double (*f)(double));
double Gauss2D2P(double (*f)(double, double));
double Gauss2D3P(double (*f)(double, double));

void CalcUniversalVals(UniversalVals* uni_vals);
void CalcJacobians(Grid* grid, UniversalVals* uni_vals);
void CalcHMatrix(Grid* grid, GlobalData* glob_data);
void CalcHbcMatrix(Grid* grid, UniversalVals* uni_vals, GlobalData* glob_data);

#endif // SRC_INTEGRALS_C_
