#include "equation.h"

#include <stdlib.h>

void InitEquation(const GlobalData* glob_data, Equation* equation) {
  equation->nn = glob_data->n_nodes;

  equation->t = (double*)calloc(equation->nn, sizeof(double));

  equation->hg = (double**)malloc(equation->nn * sizeof(double*));
  for (int i = 0; i < equation->nn; ++i) {
    equation->hg[i] = (double*)calloc(equation->nn, sizeof(double));
  }

  equation->pg = (double*)calloc(equation->nn, sizeof(double));
}

void AgregatePVectors(const GlobalData* glob_data, const Grid* grid,
                      Equation* equation) {

  for (int i = 0; i < glob_data->n_elements; ++i) {
    Element* e = &grid->elements[i];
    for (int local_id = 0; local_id < 4; ++local_id) {
      int global_id = e->nodes[local_id] - 1;
      equation->pg[global_id] += e->p_vector[local_id];
    }
  }
}

void SolveEquation(const GlobalData* glob_data, Equation* equation) { ; }

void EquationCleanup(Equation* equation) {
  free(equation->t);
  for (int i = 0; i < equation->nn; ++i) {
    free(equation->hg[i]);
  }
  free(equation->hg);
  free(equation->pg);
}
