#include "equation.h"

#include <stdlib.h>

void InitEquation(GlobalData* glob_data, Equation* equation) {
  equation->nn = glob_data->n_nodes;

  equation->t = (double*)calloc(equation->nn, sizeof(double));

  equation->hg = (double**)malloc(equation->nn * sizeof(double*));
  for (int i = 0; i < equation->nn; ++i) {
    equation->hg[i] = (double*)calloc(equation->nn, sizeof(double));
  }

  equation->pg = (double*)calloc(equation->nn, sizeof(double));
}

void SolveEquation(GlobalData* glob_data, Equation* equation) {
  ;
}

void EquationCleanup(Equation* equation) {
  free(equation->t);
  for (int i = 0; i < equation->nn; ++i) {
    free(equation->hg[i]);
  }
  free(equation->hg);
  free(equation->pg);
}
