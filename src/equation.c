#include "equation.h"

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "io.h"

void InitEquation(const GlobalData* glob_data, Equation* equation) {
  equation->nn = glob_data->n_nodes;

  equation->t = (double*)calloc(equation->nn, sizeof(double));

  equation->hg = (double**)malloc(equation->nn * sizeof(double*));
  for (int i = 0; i < equation->nn; ++i) {
    equation->hg[i] = (double*)calloc(equation->nn, sizeof(double));
  }

  equation->pg = (double*)calloc(equation->nn, sizeof(double));

  equation->c = (double**)malloc(equation->nn * sizeof(double*));
  for (int i = 0; i < equation->nn; ++i) {
    equation->c[i] = (double*)calloc(equation->nn, sizeof(double));
  }
}

void AggregatePVector(const GlobalData* glob_data, const Grid* grid, Equation* equation) {
  for (int i = 0; i < glob_data->n_elements; ++i) {
    Element* e = &grid->elements[i];
    for (int local_id = 0; local_id < 4; ++local_id) {
      int global_id = e->nodes[local_id] - 1;
      equation->pg[global_id] += e->p_vector[local_id];
    }
  }
}

void AggregateHMatrix(const GlobalData* glob_data, const Grid* grid, Equation* equation) {
  for (int i = 0; i < glob_data->n_elements; ++i) {
    Element* e = &grid->elements[i];
    for (int row = 0; row < 4; ++row) {
      int global_row = e->nodes[row] - 1;
      for (int col = 0; col < 4; ++col) {
        int global_col = e->nodes[col] - 1;
        double value = e->h_matrix[row][col] + e->hbc_matrix[row][col];
        equation->hg[global_row][global_col] += value;
      }
    }
  }
}

void AggregateCMatrix(const GlobalData* glob_data, const Grid* grid, Equation* equation) {
  for (int i = 0; i < glob_data->n_elements; ++i) {
    Element* e = &grid->elements[i];
    for (int row = 0; row < 4; ++row) {
      int global_row = e->nodes[row] - 1;
      for (int col = 0; col < 4; ++col) {
        int global_col = e->nodes[col] - 1;
        equation->c[global_row][global_col] += e->c_matrix[row][col];
      }
    }
  }
}

int GaussElimination(int n, double** a, double* b, double* x) {
  for (int k = 0; k < n - 1; k++) {
    int pivot_row = k;
    double max_val = fabs(a[k][k]);

    for (int i = k + 1; i < n; i++) {
      if (fabs(a[i][k]) > max_val) {
        max_val = fabs(a[i][k]);
        pivot_row = i;
      }
    }

    if (pivot_row != k) {
      double* temp_ptr = a[k];
      a[k] = a[pivot_row];
      a[pivot_row] = temp_ptr;

      double temp_val = b[k];
      b[k] = b[pivot_row];
      b[pivot_row] = temp_val;
    }

    if (fabs(a[k][k]) < 1e-14) {
      return -1;
    }

    for (int i = k + 1; i < n; i++) {
      double factor = a[i][k] / a[k][k];

      for (int j = k; j < n; j++) {
        a[i][j] -= factor * a[k][j];
      }
      b[i] -= factor * b[k];
    }
  }

  for (int i = n - 1; i >= 0; i--) {
    double sum = 0.0;
    for (int j = i + 1; j < n; j++) {
      sum += a[i][j] * x[j];
    }
    x[i] = (b[i] - sum) / a[i][i];
  }

  return 0;
}

void SolveSteadyState(const GlobalData* glob_data, Equation* equation) {
  int n = equation->nn;

  double** a = (double**)malloc(n * sizeof(double*));
  for (int i = 0; i < n; i++) {
    a[i] = (double*)malloc(n * sizeof(double));
  }
  double* b = (double*)malloc(n * sizeof(double));
  double* x = (double*)malloc(n * sizeof(double));

  for (int i = 0; i < n; i++) {
    b[i] = equation->pg[i];

    for (int j = 0; j < n; j++) {
      a[i][j] = equation->hg[i][j];
    }
  }

  int status = GaussElimination(n, a, b, x);

  if (status == 0) {
    for (int i = 0; i < n; i++) {
      equation->t[i] = x[i];
    }
  }

  for (int i = 0; i < n; i++) free(a[i]);
  free(a);
  free(b);
  free(x);
}

void SolveNonStationary(const GlobalData* glob_data, Equation* equation) {
  int nn = equation->nn;
  double dt = glob_data->sim_step_time;
  double total_time = glob_data->sim_time;

  for (int i = 0; i < nn; ++i) {
    equation->t[i] = glob_data->init_temp;
  }

  for (int i = 0; i < nn; ++i) {
    for (int j = 0; j < nn; ++j) {
      equation->hg[i][j] += equation->c[i][j] / dt;
    }
  }

  double* p_static = (double*)malloc(nn * sizeof(double));
  for (int i = 0; i < nn; ++i) {
    p_static[i] = equation->pg[i];
  }

  double current_time = 0.0;
  int step = 0;

  InitNonStationaryExport(glob_data);

  while (current_time < total_time) {
    current_time += dt;
    step++;

    for (int i = 0; i < nn; ++i) {
      double c_dt_t0 = 0.0;
      for (int j = 0; j < nn; ++j) {
        c_dt_t0 += (equation->c[i][j] / dt) * equation->t[j];
      }
      equation->pg[i] = p_static[i] + c_dt_t0;
    }

    SolveSteadyState(glob_data, equation);

    double min_t = equation->t[0];
    double max_t = equation->t[0];
    for (int i = 1; i < nn; ++i) {
      if (equation->t[i] < min_t) min_t = equation->t[i];
      if (equation->t[i] > max_t) max_t = equation->t[i];
    }
    ExportTempSnapshot(glob_data, equation, current_time);
  }

  free(p_static);
}

void EquationCleanup(Equation* equation) {
  free(equation->t);
  for (int i = 0; i < equation->nn; ++i) {
    free(equation->hg[i]);
  }
  free(equation->hg);
  free(equation->pg);
  for (int i = 0; i < equation->nn; ++i) {
    free(equation->c[i]);
  }
  free(equation->c);
}
