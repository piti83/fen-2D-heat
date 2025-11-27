#include "constants.h"
#include "equation.h"
#include "h_matrix.h"
#include "integrals.h"
#include "io.h"
#include "mesh.h"

int main() {
  InitConstants();

  GlobalData data;
  Grid grid;

  ReadFile("run/test.txt", &data, &grid);
  PrintInfo(&data, &grid);

  UniversalVals uni_vals;
  CalcUniversalVals(&uni_vals);
  CalcJacobians(&grid, &uni_vals);
  CalcHMatrix(&grid, &data);

#ifdef DEBUG
  ExportJacobianData(&grid, &uni_vals);
#endif  // DEBUG

  GlobHMatrix h_matrix;
  InitHMatrix(&h_matrix, &data);
  CalcGlobalHMatrix(&h_matrix, &grid);

#ifdef DEBUG
  PrintGlobalH(&h_matrix);
#endif  // DEBUG

  CalcHbcMatrix(&grid, &uni_vals, &data);
#ifdef DEBUG
  ExportHbcMatrices(&grid);
#endif

  CalcPVector(&grid, &uni_vals, &data);

  Equation eq;
  InitEquation(&data, &eq);
  AgregatePVectors(&data, &grid, &eq);

  SolveEquation(&data, &eq);

  EquationCleanup(&eq);
  HMatrixCleanup(&h_matrix);
  GridCleanup(&grid);

  return 0;
}
