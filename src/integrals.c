#include "integrals.h"

#include <math.h>
#include <stdio.h>

#include "constants.h"
#include "mesh.h"

double N1(double ksi, double eta) { return 0.25 * (1 - ksi) * (1 - eta); }

double N2(double ksi, double eta) { return 0.25 * (1 + ksi) * (1 - eta); }

double N3(double ksi, double eta) { return 0.25 * (1 + ksi) * (1 + eta); }

double N4(double ksi, double eta) { return 0.25 * (1 - ksi) * (1 + eta); }

double Gauss1D2P(double (*f)(double)) {
  return kWeightsN2[0] * f(kPointsN2[0]) + kWeightsN2[1] * f(kPointsN2[1]);
}

void CalcUniversalVals(UniversalVals* uni_vals, int nip) {
  double* points = GetPoints(nip);

  if (!points) {
    printf("Error. Unsupported NIP: %d\n", nip);
    return;
  }

  int idx = 0;

  for (int i = 0; i < nip; ++i) {
    for (int j = 0; j < nip; ++j) {
      double eta = points[i];
      double ksi = points[j];

      uni_vals->dn_dksi[idx][0] = -0.25 * (1 - eta);
      uni_vals->dn_dksi[idx][1] = 0.25 * (1 - eta);
      uni_vals->dn_dksi[idx][2] = 0.25 * (1 + eta);
      uni_vals->dn_dksi[idx][3] = -0.25 * (1 + eta);

      uni_vals->dn_deta[idx][0] = -0.25 * (1 - ksi);
      uni_vals->dn_deta[idx][1] = -0.25 * (1 + ksi);
      uni_vals->dn_deta[idx][2] = 0.25 * (1 + ksi);
      uni_vals->dn_deta[idx][3] = 0.25 * (1 - ksi);

      idx++;
    }
  }
}

void CalcJacobians(Grid* grid, UniversalVals* uni_vals, GlobalData* glob_data) {
  int nip = glob_data->nip_elem;
  double* weights = GetWeights(nip);

  for (int i = 0; i < grid->n_elements; ++i) {
    int idx = 0;
    for (int y = 0; y < nip; ++y) {
      for (int x = 0; x < nip; ++x) {
        int ip = idx;
        grid->elements[i].jacobian[ip].weight = weights[y] * weights[x];

        double dx_dksi = 0.0, dx_deta = 0.0, dy_dksi = 0.0, dy_deta = 0.0;

        for (int k = 0; k < 4; ++k) {
          double x_node = grid->nodes[grid->elements[i].nodes[k] - 1].x;
          double y_node = grid->nodes[grid->elements[i].nodes[k] - 1].y;

          dx_dksi += uni_vals->dn_dksi[ip][k] * x_node;
          dx_deta += uni_vals->dn_deta[ip][k] * x_node;
          dy_dksi += uni_vals->dn_dksi[ip][k] * y_node;
          dy_deta += uni_vals->dn_deta[ip][k] * y_node;
        }

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

          grid->elements[i].jacobian[ip].dN_dx[k] = j1_00 * dN_dksi_k + j1_01 * dN_deta_k;

          grid->elements[i].jacobian[ip].dN_dy[k] = j1_10 * dN_dksi_k + j1_11 * dN_deta_k;
        }
        idx++;
      }
    }
  }
  grid->jacobians_calculated = 1;
}

void CalcHMatrices(Grid* grid, GlobalData* glob_data) {
  if (!grid->jacobians_calculated) {
    printf("Call CalcJacobians before calculating H\n");
    return;
  }

  int nip = glob_data->nip_elem;
  int total_points = nip * nip;

  for (int e = 0; e < grid->n_elements; ++e) {
    for (int i = 0; i < 4; ++i)
      for (int j = 0; j < 4; ++j) grid->elements[e].h_matrix[i][j] = 0.0;

    for (int ip = 0; ip < total_points; ++ip) {
      Jacobian* c_j = &grid->elements[e].jacobian[ip];
      double (*c_h)[4] = grid->elements[e].h_matrix;

      double factor = glob_data->conductivity * c_j->det_j * c_j->weight;

      for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
          c_h[i][j] += factor * (c_j->dN_dx[i] * c_j->dN_dx[j] + c_j->dN_dy[i] * c_j->dN_dy[j]);
        }
      }
    }
  }
}

void CalcHbcMatrices(Grid* grid, UniversalVals* uni_vals, GlobalData* glob_data) {
  int nip = glob_data->nip_bc;
  double* points = GetPoints(nip);
  double* weights = GetWeights(nip);

  for (int i = 0; i < grid->n_elements; ++i) {
    Element* e = &grid->elements[i];

    for (int hbc_i = 0; hbc_i < 4; ++hbc_i) {
      for (int hbc_j = 0; hbc_j < 4; ++hbc_j) {
        e->hbc_matrix[hbc_i][hbc_j] = 0.0;
      }
    }

    int id1 = e->nodes[0] - 1;
    int id2 = e->nodes[1] - 1;
    int id3 = e->nodes[2] - 1;
    int id4 = e->nodes[3] - 1;

    int side_nodes[4][2] = {{id1, id2}, {id2, id3}, {id3, id4}, {id4, id1}};

    for (int side = 0; side < 4; ++side) {
      Node* n_a = &grid->nodes[side_nodes[side][0]];
      Node* n_b = &grid->nodes[side_nodes[side][1]];

      if (n_a->bc == false || n_b->bc == false) {
        continue;
      }

      double det_j_side = sqrt(pow(n_b->x - n_a->x, 2) + pow(n_b->y - n_a->y, 2)) / 2.0;

      for (int ip = 0; ip < nip; ++ip) {
        double ip_val = points[ip];
        double weight = weights[ip];

        double ksi, eta;

        switch (side) {
          case 0:
            ksi = ip_val;
            eta = -1.0;
            break;
          case 1:
            ksi = 1.0;
            eta = ip_val;
            break;
          case 2:
            ksi = -ip_val;
            eta = 1.0;
            break;  // Uwaga: kierunek? Dla symetrii Gaussa ok
          case 3:
            ksi = -1.0;
            eta = -ip_val;
            break;
        }

        double N[4];
        N[0] = N1(ksi, eta);
        N[1] = N2(ksi, eta);
        N[2] = N3(ksi, eta);
        N[3] = N4(ksi, eta);

        // Zapisz funkcje kształtu dla ewentualnego debugowania (jeśli Surface obsłuży >2 pkt)
        for (int i_surface = 0; i_surface < 4; ++i_surface) {
          if (ip < MAX_NIP_1D) e->surface[side].n[ip][i_surface] = N[i_surface];
        }

        double coefficient = glob_data->alfa * weight * det_j_side;

        for (int i_hbc = 0; i_hbc < 4; ++i_hbc) {
          for (int j_hbc = 0; j_hbc < 4; ++j_hbc) {
            e->hbc_matrix[i_hbc][j_hbc] += coefficient * N[i_hbc] * N[j_hbc];
          }
        }
      }
    }
  }
}

void CalcPVectors(Grid* grid, UniversalVals* uni_vals, GlobalData* glob_data) {
  int nip = glob_data->nip_bc;
  double* points = GetPoints(nip);
  double* weights = GetWeights(nip);

  for (int i = 0; i < grid->n_elements; ++i) {
    for (int j = 0; j < 4; ++j) {
      grid->elements[i].p_vector[j] = 0;
    }
  }

  for (int i = 0; i < grid->n_elements; ++i) {
    Element* e = &grid->elements[i];

    int id1 = e->nodes[0] - 1;
    int id2 = e->nodes[1] - 1;
    int id3 = e->nodes[2] - 1;
    int id4 = e->nodes[3] - 1;

    int side_nodes[4][2] = {{id1, id2}, {id2, id3}, {id3, id4}, {id4, id1}};

    for (int side = 0; side < 4; ++side) {
      Node* n_a = &grid->nodes[side_nodes[side][0]];
      Node* n_b = &grid->nodes[side_nodes[side][1]];

      if (n_a->bc == false || n_b->bc == false) {
        continue;
      }

      double length = sqrt(pow(n_b->x - n_a->x, 2) + pow(n_b->y - n_a->y, 2));
      double det_j_side = length / 2.0;

      for (int ip = 0; ip < nip; ++ip) {
        double weight = weights[ip];
        double point = points[ip];

        double ksi = 0.0;
        double eta = 0.0;

        switch (side) {
          case 0:
            ksi = point;
            eta = -1.0;
            break;
          case 1:
            ksi = 1.0;
            eta = point;
            break;
          case 2:
            ksi = -point;
            eta = 1.0;
            break;
          case 3:
            ksi = -1.0;
            eta = -point;
            break;
        }

        double N[4];
        N[0] = N1(ksi, eta);
        N[1] = N2(ksi, eta);
        N[2] = N3(ksi, eta);
        N[3] = N4(ksi, eta);

        double coefficient = weight * glob_data->alfa * glob_data->tot * det_j_side;

        for (int k = 0; k < 4; ++k) {
          e->p_vector[k] += N[k] * coefficient;
        }
      }
    }
  }
}

void CalcCMatrices(Grid* grid, UniversalVals* uni_vals, GlobalData* glob_data) {

}
