#ifndef SRC_EQUATION_H_
#define SRC_EQUATION_H_

#include "mesh.h"

typedef struct {
  int nn;
  double* t;
  double** hg;
  double* pg;
} Equation;

void InitEquation(GlobalData* glob_data, Equation* equation);
void SolveEquation(GlobalData* glob_data, Equation* equation);
void EquationCleanup(Equation* equation);

#endif
