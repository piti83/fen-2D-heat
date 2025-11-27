#ifndef SRC_EQUATION_H_
#define SRC_EQUATION_H_

#include "mesh.h"

typedef struct {
  int nn;
  double* t;
  double** hg;
  double* pg;
} Equation;

void InitEquation(const GlobalData* glob_data, Equation* equation);
void AgregatePVectors(const GlobalData* glob_data, const Grid* grid, Equation* equation);
void SolveEquation(const GlobalData* glob_data, Equation* equation);
void EquationCleanup(Equation* equation);

#endif
