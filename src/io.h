#ifndef SRC_INPUT_H_
#define SRC_INPUT_H_

#include "equation.h"
#include "mesh.h"

void ReadFile(const char* path, GlobalData* glob_data, Grid* grid);

void PrintInfo(const GlobalData* glob_data, const Grid* grid);

void ExportJacobianData(const Grid* grid, const UniversalVals* uni_vals);

void ExportGlobalH(const Equation* equation);
void ExportGlobalP(const Equation* equation);
void ExportGlobalC(const Equation* equation);

void ExportHbcMatrices(const Grid* grid);
void ExportPVectors(const Grid* grid);
void ExportCMatrices(const Grid* grid);

void ExportTemperatureVector(const Equation* equation);

void InitNonStationaryExport(const GlobalData* data);
void ExportTempSnapshot(const GlobalData* data, const Equation* equation, double current_time);

#endif  // SRC_INPUT_H_
