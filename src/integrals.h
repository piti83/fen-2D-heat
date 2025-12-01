#ifndef SRC_INTEGRALS_C_
#define SRC_INTEGRALS_C_

#include <mesh.h>

double N1(double ksi, double eta);
double N2(double ksi, double eta);
double N3(double ksi, double eta);
double N4(double ksi, double eta);

double Gauss1D2P(double (*f)(double));
double Gauss1D3P(double (*f)(double));
double Gauss2D2P(double (*f)(double, double));
double Gauss2D3P(double (*f)(double, double));

void CalcUniversalVals(UniversalVals* uni_vals, int nip);
void CalcJacobians(Grid* grid, UniversalVals* uni_vals, GlobalData* glob_data);
void CalcHMatrices(Grid* grid, GlobalData* glob_data);
void CalcHbcMatrices(Grid* grid, UniversalVals* uni_vals, GlobalData* glob_data);
void CalcPVectors(Grid* grid, UniversalVals* uni_vals, GlobalData* glob_data);
void CalcCMatrices(Grid* grid, UniversalVals* uni_vals, GlobalData* glob_data);

#endif // SRC_INTEGRALS_C_
