#include "io.h"

#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void ReadFile(const char* path, GlobalData* glob_data, Grid* grid) {
  FILE* f_ptr = fopen(path, "r");
  if (!f_ptr) {
    printf("Failed to open file: %s\n", path);
    return;
  }

  char line[256];
  // 0=None, 1=Node, 2=Element, 3=BC
  int current_section = 0;

  while (fgets(line, sizeof(line), f_ptr)) {
    if (strstr(line, "*Node")) {
      current_section = 1;
      continue;
    }
    if (strstr(line, "*Element")) {
      current_section = 2;
      continue;
    }
    if (strstr(line, "*BC")) {
      current_section = 3;
      continue;
    }

    if (strstr(line, "SimulationTime"))
      sscanf(line, "%*s %lf", &glob_data->sim_time);
    else if (strstr(line, "SimulationStepTime"))
      sscanf(line, "%*s %lf", &glob_data->sim_step_time);
    else if (strstr(line, "Conductivity"))
      sscanf(line, "%*s %lf", &glob_data->conductivity);
    else if (strstr(line, "Alfa"))
      sscanf(line, "%*s %lf", &glob_data->alfa);
    else if (strstr(line, "Tot"))
      sscanf(line, "%*s %lf", &glob_data->tot);
    else if (strstr(line, "InitialTemp"))
      sscanf(line, "%*s %lf", &glob_data->init_temp);
    else if (strstr(line, "Density"))
      sscanf(line, "%*s %lf", &glob_data->density);
    else if (strstr(line, "SpecificHeat"))
      sscanf(line, "%*s %lf", &glob_data->spec_heat);

    else if (strstr(line, "Nodes number")) {
      sscanf(line, "%*s %*s %d", &grid->n_nodes);
      glob_data->n_nodes = grid->n_nodes;
      grid->nodes = (Node*)calloc(grid->n_nodes, sizeof(Node));
    } else if (strstr(line, "Elements number")) {
      sscanf(line, "%*s %*s %d", &grid->n_elements);
      glob_data->n_elements = grid->n_elements;
      grid->elements = (Element*)calloc(grid->n_elements, sizeof(Element));
    }

    else if (current_section == 1 && grid->nodes) {
      int id;
      double x, y;
      if (sscanf(line, "%d, %lf, %lf", &id, &x, &y) == 3) {
        grid->nodes[id - 1].x = x;
        grid->nodes[id - 1].y = y;
        grid->nodes[id - 1].bc = false;
      }
    } else if (current_section == 2 && grid->elements) {
      int id, n1, n2, n3, n4;
      if (sscanf(line, "%d, %d, %d, %d, %d", &id, &n1, &n2, &n3, &n4) == 5) {
        grid->elements[id - 1].nodes[0] = n1;
        grid->elements[id - 1].nodes[1] = n2;
        grid->elements[id - 1].nodes[2] = n3;
        grid->elements[id - 1].nodes[3] = n4;
      }
    } else if (current_section == 3 && grid->nodes) {
      char* token = strtok(line, ", ");
      while (token) {
        int id = atoi(token);
        if (id > 0 && id <= grid->n_nodes) {
          grid->nodes[id - 1].bc = true;
        }
        token = strtok(NULL, ", ");
      }
    }
  }

  fclose(f_ptr);
}

void PrintInfo(const GlobalData* glob_data, const Grid* grid) {
  printf("=============== GLOBAL DATA ===============\n\n");

  printf("- Simulation Time      :  %f\n", glob_data->sim_time);
  printf("- Simulation Step Time :  %f\n", glob_data->sim_step_time);
  printf("- Conductivity         :  %f\n", glob_data->conductivity);
  printf("- Alfa                 :  %f\n", glob_data->alfa);
  printf("- Tot                  :  %f\n", glob_data->tot);
  printf("- Initial Temperature  :  %f\n", glob_data->init_temp);
  printf("- Density              :  %f\n", glob_data->density);
  printf("- Specific Heat        :  %f\n", glob_data->spec_heat);
  printf("- Nodes                :  %d\n", glob_data->n_nodes);
  printf("- Elements             :  %d\n", glob_data->n_elements);
  printf("- NIP (Elements)       :  %d\n", glob_data->nip_elem);
  printf("- NIP (Surface)        :  %d\n\n", glob_data->nip_bc);

  printf("===========================================\n\n");

  printf("                   NODES                   \n");
  printf("+-----------------------------------------+\n");
  for (int i = 0; i < grid->n_nodes; ++i) {
    printf("|%-7d|%16.10lf|%16.10lf|\n", i + 1, grid->nodes[i].x, grid->nodes[i].y);
  }
  printf("+-----------------------------------------+\n\n");

  printf("                  ELEMENTS                 \n");
  printf("+-----------------------------------------+\n");
  for (int i = 0; i < grid->n_elements; ++i) {
    printf("|%-9d|%7d|%7d|%7d|%7d|\n", i + 1, grid->elements[i].nodes[0],
           grid->elements[i].nodes[1], grid->elements[i].nodes[2], grid->elements[i].nodes[3]);
  }
  printf("+-----------------------------------------+\n");

  printf("                     BC                    \n");
  printf("+-----------------------------------------+\n");
  for (int i = 0; i < grid->n_nodes; ++i) {
    if (grid->nodes[i].bc) {
      printf("%d ", i + 1);
    }
  }
  printf("\n");
  printf("+-----------------------------------------+\n");
}

void ExportJacobianData(const Grid* grid, const UniversalVals* uni_vals) {
  FILE* fptr;
  fptr = fopen("out/jacobian.txt", "w");

  fprintf(fptr, "==================== CALCULATED CONSTANTS ====================\n\n");
  fprintf(fptr, "dN/dKsi\n");
  for (int i = 0; i < 4; ++i) {
    fprintf(fptr, "\t%12.8f %12.8f %12.8f %12.8f\n", uni_vals->dn_dksi[i][0],
            uni_vals->dn_dksi[i][1], uni_vals->dn_dksi[i][2], uni_vals->dn_dksi[i][3]);
  }
  fprintf(fptr, "dN/dEta\n");
  for (int i = 0; i < 4; ++i) {
    fprintf(fptr, "\t%12.8f %12.8f %12.8f %12.8f\n", uni_vals->dn_deta[i][0],
            uni_vals->dn_deta[i][1], uni_vals->dn_deta[i][2], uni_vals->dn_deta[i][3]);
  }

  fprintf(fptr, "\n\n============== RESULTS FOR POINTS ==============\n\n");

  for (int i = 0; i < grid->n_elements; ++i) {
    fprintf(fptr, "Element nr %d\n", i + 1);
    fprintf(fptr, "\t Macierz H:\n");
    for (int h_i = 0; h_i < 4; ++h_i) {
      fprintf(fptr, "\t\t");
      for (int h_j = 0; h_j < 4; ++h_j) {
        fprintf(fptr, "%5.4lf ", grid->elements[i].h_matrix[h_i][h_j]);
      }
      fprintf(fptr, "\n");
    }
    fprintf(fptr, "\n");
    for (int ip = 0; ip < 4; ++ip) {
      Jacobian* j = &grid->elements[i].jacobian[ip];

      fprintf(fptr, "\tJakobi matrix for %d integration point\n", ip + 1);

      fprintf(fptr, "\t\t%12.8lf  %12.8lf\n", j->j[0][0], j->j[0][1]);
      fprintf(fptr, "\t\t%12.8lf  %12.8lf\n", j->j[1][0], j->j[1][1]);

      fprintf(fptr, "\t\tdet[J] = %-12.8lf\n\n", j->det_j);

      fprintf(fptr, "\t\tdN/dx\n");
      fprintf(fptr, "\t\t\t%12.5lf %12.5lf %12.5lf %12.5lf\n\n", j->dN_dx[0], j->dN_dx[1],
              j->dN_dx[2], j->dN_dx[3]);

      fprintf(fptr, "\t\tdN/dy\n");
      fprintf(fptr, "\t\t\t%12.5lf %12.5lf %12.5lf %12.5lf\n\n", j->dN_dy[0], j->dN_dy[1],
              j->dN_dy[2], j->dN_dy[3]);
    }
  }
  fprintf(fptr, "\n");
  fclose(fptr);
}

void ExportGlobalH(const Equation* equation) {
  FILE* fptr = fopen("out/global_h_matrix.txt", "w");
  fprintf(fptr, "[GLOBAL H MATRIX]\n\n");
  for (int i = 0; i < equation->nn; ++i) {
    for (int j = 0; j < equation->nn; ++j) {
      fprintf(fptr, "%6.2lf  ", equation->hg[i][j]);
    }
    fprintf(fptr, "\n");
  }
  fclose(fptr);
}

void ExportGlobalP(const Equation* equation) {
  FILE* fptr = fopen("out/global_p_vector.txt", "w");
  fprintf(fptr, "[GLOBAL P VECTOR]\n\n");
  for (int i = 0; i < equation->nn; ++i) {
    fprintf(fptr, "%6.2lf  ", equation->pg[i]);
  }
  fprintf(fptr, "\n");
  fclose(fptr);
}

void ExportGlobalC(const Equation* equation) {
  FILE* fptr = fopen("out/global_c_matrix.txt", "w");

  fprintf(fptr, "[GLOBAL C MATRIX]\n\n");
  for (int i = 0; i < equation->nn; ++i) {
    for (int j = 0; j < equation->nn; ++j) {
      fprintf(fptr, "%6.2lf  ", equation->c[i][j]);
    }
    fprintf(fptr, "\n");
  }

  fclose(fptr);
}

void ExportHbcMatrices(const Grid* grid) {
  FILE* fptr;
  fptr = fopen("out/hbc.txt", "w");

  fprintf(fptr, "[HBC MATRICES]\n\n");

  for (int i = 0; i < grid->n_elements; ++i) {
    Element* e = &grid->elements[i];
    fprintf(fptr, "\tElement %d:\n\n", i + 1);
    for (int hbc_i = 0; hbc_i < 4; ++hbc_i) {
      fprintf(fptr, "\t\t");
      for (int hbc_j = 0; hbc_j < 4; ++hbc_j) {
        fprintf(fptr, "%6.3lf   ", e->hbc_matrix[hbc_i][hbc_j]);
      }
      fprintf(fptr, "\n");
    }
    fprintf(fptr, "\n");
  }

  fclose(fptr);
}

void ExportPVectors(const Grid* grid) {
  FILE* fptr = fopen("out/local_p_vectors.txt", "w");

  fprintf(fptr, "[LOCAL P VECTORS]\n\n");
  for (int i = 0; i < grid->n_elements; ++i) {
    Element* e = &grid->elements[i];
    fprintf(fptr, "\tElement %d:\n\n\t\t", i + 1);
    for (int p_i = 0; p_i < 4; ++p_i) {
      fprintf(fptr, "%6.3lf   ", e->p_vector[p_i]);
    }
    fprintf(fptr, "\n");
  }

  fclose(fptr);
}

void ExportCMatrices(const Grid* grid) {
  FILE* fptr = fopen("out/local_c_matrices.txt", "w");

  fprintf(fptr, "[LOCAL C MATRICES]\n\n");

  for (int i = 0; i < grid->n_elements; ++i) {
    Element* e = &grid->elements[i];
    fprintf(fptr, "\tElement %d:\n\n", i + 1);
    for (int c_i = 0; c_i < 4; ++c_i) {
      fprintf(fptr, "\t\t");
      for (int c_j = 0; c_j < 4; ++c_j) {
        fprintf(fptr, "%6.3lf   ", e->c_matrix[c_i][c_j]);
      }
      fprintf(fptr, "\n");
    }
    fprintf(fptr, "\n");
  }

  fclose(fptr);
}

void ExportTemperatureVector(const Equation* equation) {
  FILE* fptr = fopen("out/temperature_vector.txt", "w");
  fprintf(fptr, "[TEMPERATURE VECTOR]\n\n");
  for (int i = 0; i < equation->nn; ++i) {
    fprintf(fptr, "%6.2lf  ", equation->t[i]);
  }
  fprintf(fptr, "\n");
  fclose(fptr);
}

void InitNonStationaryExport(const GlobalData* data) {
  FILE* fptr = fopen("out/non_stat_temperatures.txt", "w");
  fprintf(fptr, "[NON STATIONARY SOLUTION]\n\n");
  fprintf(fptr, "Time      |");
  for (int i = 0; i < data->n_nodes; ++i) {
    fprintf(fptr, " n%-9d|", i + 1);
  }
  fprintf(fptr, "\n");
  fclose(fptr);
}

void ExportTempSnapshot(const GlobalData* data, const Equation* equation, double current_time) {
  FILE* fptr = fopen("out/non_stat_temperatures.txt", "a");
  fprintf(fptr, "%-10.2lf|", current_time);
  for (int i = 0; i < data->n_nodes; ++i) {
    fprintf(fptr, " %10.4lf|", equation->t[i]);
  }
  fprintf(fptr, "\n");
  fclose(fptr);
}
