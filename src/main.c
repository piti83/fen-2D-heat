#include <stdio.h>
#include <stdlib.h>

#include "constants.h"
#include "equation.h"
#include "integrals.h"
#include "io.h"
#include "mesh.h"

int main(int argc, char* argv[]) {
  InitConstants();

  GlobalData data;
  Grid grid;

  int nip_elem = 2;
  int nip_bc = 2;

  if (argc >= 3) {
    nip_elem = atoi(argv[1]);
    nip_bc = atoi(argv[2]);
    printf("NIP set from cl arguments: Element=%d, Surface=%d\n", nip_elem, nip_bc);
  } else {
    printf("Usage: %s <nip_element> <nip_surface>\n", argv[0]);
    printf("Using default values: 2 2\n");
  }

  if (nip_elem < 1 || nip_elem > 4 || nip_bc < 1 || nip_bc > 4) {
    printf("Error: NIP has to be in range [1, 4]\n");
    return 1;
  }

  ReadFile("run/test.txt", &data, &grid);

  data.nip_elem = nip_elem;
  data.nip_bc = nip_bc;

#ifdef DEBUG
  PrintInfo(&data, &grid);
#endif

  UniversalVals uni_vals;
  CalcUniversalVals(&uni_vals, data.nip_elem);
  CalcJacobians(&grid, &uni_vals, &data);
  CalcHMatrices(&grid, &data);

#ifdef DEBUG
  ExportJacobianData(&grid, &uni_vals);
#endif

  CalcHbcMatrices(&grid, &uni_vals, &data);

#ifdef DEBUG
  ExportHbcMatrices(&grid);
#endif

  CalcPVectors(&grid, &uni_vals, &data);

#ifdef DEBUG
  ExportPVectors(&grid);
#endif

  CalcCMatrices(&grid, &uni_vals, &data);

#ifdef DEBUG
  ExportCMatrices(&grid);
#endif

  Equation eq;
  InitEquation(&data, &eq);

  AggregatePVector(&data, &grid, &eq);
  AggregateHMatrix(&data, &grid, &eq);
  AggregateCMatrix(&data, &grid, &eq);

#ifdef DEBUG
  ExportGlobalP(&eq);
  ExportGlobalH(&eq);
  ExportGlobalC(&eq);
#endif

  SolveEquation(&data, &eq);

  EquationCleanup(&eq);
  GridCleanup(&grid);

  return 0;
}
