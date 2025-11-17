#include <stdio.h>

#include "constants.h"
#include "integrals.h"

double Gauss1D2P(double (*f)(double)) {
  return kWeightsN2[0] * f(kPointsN2[0]) + kWeightsN2[1] * f(kPointsN2[1]);
}

double Gauss1D3P(double (*f)(double)) {
  double result = 0.0;
  for (int i = 0; i < 3; ++i) {
    result += kWeightsN3[i] * f(kPointsN3[i]);
  }
  return result;
}

double Gauss2D2P(double (*f)(double, double)) {
  double result = 0.0;
  for (int i = 0; i < 2; ++i) {
    for (int j = 0; j < 2; ++j) {
      result += kWeightsN2[i] * kWeightsN2[j] * f(kPointsN2[i], kPointsN2[j]);
    }
  }
  return result;
}

double Gauss2D3P(double (*f)(double, double)) {
  double result = 0.0;
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      result += kWeightsN3[i] * kWeightsN3[j] * f(kPointsN3[i], kPointsN3[j]);
    }
  }
  return result;
}

void CalcUniversalVals(UniversalVals* uni_vals) {
  for (int i = 0; i < 4; ++i) {
    uni_vals->dn_dksi[i][0] = -0.25 * (1 - kPointsN2D2[i].y);
    uni_vals->dn_dksi[i][1] = 0.25 * (1 - kPointsN2D2[i].y);
    uni_vals->dn_dksi[i][2] = 0.25 * (1 + kPointsN2D2[i].y);
    uni_vals->dn_dksi[i][3] = -0.25 * (1 + kPointsN2D2[i].y);
  }
  for (int i = 0; i < 4; ++i) {
    uni_vals->dn_deta[i][0] = -0.25 * (1 - kPointsN2D2[i].x);
    uni_vals->dn_deta[i][1] = -0.25 * (1 + kPointsN2D2[i].x);
    uni_vals->dn_deta[i][2] = 0.25 * (1 + kPointsN2D2[i].x);
    uni_vals->dn_deta[i][3] = 0.25 * (1 - kPointsN2D2[i].x);
  }
}

void CalcJacobians(Grid* grid, UniversalVals* uni_vals) {
  for (int i = 0; i < grid->n_elements; ++i) {
    for (int ip = 0; ip < 4; ++ip) {
      int id1 = grid->elements[i].nodes[0] - 1;
      int id2 = grid->elements[i].nodes[1] - 1;
      int id3 = grid->elements[i].nodes[2] - 1;
      int id4 = grid->elements[i].nodes[3] - 1;

      double dx_dksi = uni_vals->dn_dksi[ip][0] * grid->nodes[id1].x +
                       uni_vals->dn_dksi[ip][1] * grid->nodes[id2].x +
                       uni_vals->dn_dksi[ip][2] * grid->nodes[id3].x +
                       uni_vals->dn_dksi[ip][3] * grid->nodes[id4].x;

      double dx_deta = uni_vals->dn_deta[ip][0] * grid->nodes[id1].x +
                       uni_vals->dn_deta[ip][1] * grid->nodes[id2].x +
                       uni_vals->dn_deta[ip][2] * grid->nodes[id3].x +
                       uni_vals->dn_deta[ip][3] * grid->nodes[id4].x;

      double dy_dksi = uni_vals->dn_dksi[ip][0] * grid->nodes[id1].y +
                       uni_vals->dn_dksi[ip][1] * grid->nodes[id2].y +
                       uni_vals->dn_dksi[ip][2] * grid->nodes[id3].y +
                       uni_vals->dn_dksi[ip][3] * grid->nodes[id4].y;

      double dy_deta = uni_vals->dn_deta[ip][0] * grid->nodes[id1].y +
                       uni_vals->dn_deta[ip][1] * grid->nodes[id2].y +
                       uni_vals->dn_deta[ip][2] * grid->nodes[id3].y +
                       uni_vals->dn_deta[ip][3] * grid->nodes[id4].y;

      grid->elements[i].jacobian[ip].j[0][0] = dx_dksi;
      grid->elements[i].jacobian[ip].j[0][1] = dy_dksi;
      grid->elements[i].jacobian[ip].j[1][0] = dx_deta;
      grid->elements[i].jacobian[ip].j[1][1] = dy_deta;

      double det_j = (dx_dksi * dy_deta) - (dy_dksi * dx_deta);
      grid->elements[i].jacobian[ip].det_j = det_j;

      double j1_00 = dy_deta / det_j;
      double j1_01 = -dy_dksi / det_j;
      double j1_10 = -dx_deta / det_j;
      double j1_11 = dx_dksi / det_j;

      grid->elements[i].jacobian[ip].j1[0][0] = j1_00;
      grid->elements[i].jacobian[ip].j1[0][1] = j1_01;
      grid->elements[i].jacobian[ip].j1[1][0] = j1_10;
      grid->elements[i].jacobian[ip].j1[1][1] = j1_11;

      for (int k = 0; k < 4; ++k) {
        double dN_dksi_k = uni_vals->dn_dksi[ip][k];
        double dN_deta_k = uni_vals->dn_deta[ip][k];

        grid->elements[i].jacobian[ip].dN_dx[k] =
            j1_00 * dN_dksi_k + j1_01 * dN_deta_k;

        grid->elements[i].jacobian[ip].dN_dy[k] =
            j1_10 * dN_dksi_k + j1_11 * dN_deta_k;
      }
    }
  }
  grid->jacobians_calculated = 1;
}

void CalcHMatrix(Grid* grid, GlobalData* glob_data) {
  if (!grid->jacobians_calculated) {
    printf("Call CalcJacobians before calculating H\n");
    return;
  }

  for (int e = 0; e < grid->n_elements; ++e) {
    for (int ip = 0; ip < 4; ++ip) {
      Jacobian* c_j = &grid->elements[e].jacobian[ip];
      double (*c_h)[4] = grid->elements[e].h_matrix;

      for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
          c_h[i][j] +=
              glob_data->conductivity *
              (c_j->dN_dx[i] * c_j->dN_dx[j] + c_j->dN_dy[i] * c_j->dN_dy[j]) *
              c_j->det_j;
        }
      }
    }
  }
}

void CalcHbcMatrix(Grid* grid, UniversalVals* uni_vals, GlobalData* glob_data) {
  ;
}
