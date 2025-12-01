#ifndef SRC_MESH_H_
#define SRC_MESH_H_

#include <stdint.h>
#include <constants.h>
#include <stdbool.h>

#define MAX_NIP_1D 4
#define MAX_NIP_2D (MAX_NIP_1D * MAX_NIP_1D)

typedef struct {
  double x;
  double y;
  bool bc;
} Node;

typedef struct {
  double j[2][2];
  double j1[2][2];
  double det_j;
  double dN_dx[4];
  double dN_dy[4];
  double weight;
} Jacobian;

typedef struct {
  double n[MAX_NIP_1D][4];
} Surface;

typedef struct {
  unsigned int nodes[4];
  Jacobian jacobian[MAX_NIP_2D];
  double h_matrix[4][4];
  double hbc_matrix[4][4];
  Surface surface[4];
  double p_vector[4];
  double c_matrix[4][4];
} Element;

typedef struct {
  unsigned int n_nodes;
  unsigned int n_elements;
  Node* nodes;
  Element* elements;
  int16_t jacobians_calculated;
} Grid;

typedef struct {
  double dn_dksi[MAX_NIP_2D][4];
  double dn_deta[MAX_NIP_2D][4];
} UniversalVals;

typedef struct {
  double sim_time;
  double sim_step_time;
  double conductivity;
  double alfa;
  double tot;
  double init_temp;
  double density;
  double spec_heat;
  int n_nodes;
  int n_elements;
  int nip_elem;
  int nip_bc;
} GlobalData;

void GridCleanup(Grid* grid);

#endif  // SRC_MESH_H_
