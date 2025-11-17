#include "constants.h"
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

  GridCleanup(&grid);
  HMatrixCleanup(&h_matrix);

  return 0;
}
