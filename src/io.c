#include "io.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "constants.h"

void ReadFile(const char* path, GlobalData* glob_data, Grid* grid) {
  FILE* f_ptr = fopen(path, "r");

  if (f_ptr == NULL) {
    printf("Failed to open file: %s. Aborting.", path);
    return;
  }

  if (fscanf(f_ptr, " %*[^0-9.-]%f", &glob_data->sim_time) != 1) {
    printf("Error reading sim_time.\n");
  }
  if (fscanf(f_ptr, " %*[^0-9.-]%f", &glob_data->sim_step_time) != 1) {
    printf("Error reading sim_step_time.\n");
  }
  if (fscanf(f_ptr, " %*[^0-9.-]%f", &glob_data->conductivity) != 1) {
    printf("Error reading conductivity.\n");
  }
  if (fscanf(f_ptr, " %*[^0-9.-]%f", &glob_data->alfa) != 1) {
    printf("Error reading alfa.\n");
  }
  if (fscanf(f_ptr, " %*[^0-9.-]%f", &glob_data->tot) != 1) {
    printf("Error reading tot.\n");
  }
  if (fscanf(f_ptr, " %*[^0-9.-]%f", &glob_data->init_temp) != 1) {
    printf("Error reading init_temp.\n");
  }
  if (fscanf(f_ptr, " %*[^0-9.-]%f", &glob_data->density) != 1) {
    printf("Error reading density.\n");
  }
  if (fscanf(f_ptr, " %*[^0-9.-]%f", &glob_data->spec_heat) != 1) {
    printf("Error reading spec_heat.\n");
  }
  if (fscanf(f_ptr, " %*[^0-9]%d", &glob_data->n_nodes) != 1) {
    printf("Error reading n_nodes.\n");
  }
  if (fscanf(f_ptr, " %*[^0-9]%d", &glob_data->n_elements) != 1) {
    printf("Error reading n_elements.\n");
  }

  grid->n_nodes = glob_data->n_nodes;
  grid->n_elements = glob_data->n_elements;

  grid->nodes = (Node*)calloc(grid->n_nodes, sizeof(Node));
  grid->elements = (Element*)calloc(grid->n_elements, sizeof(Element));

  if (!grid->nodes) {
    printf("Failed to allocate memory for %d nodes.", grid->n_nodes);
  }
  if (!grid->elements) {
    printf("Failed to allocate memory for %d elements.", grid->n_elements);
  }

  char line_buffer[100];
  if (!fgets(line_buffer, sizeof line_buffer, f_ptr)) {
    printf("Error during fgets()");
  }
  if (!fgets(line_buffer, sizeof line_buffer, f_ptr)) {
    printf("Error during fgets()");
  }
  for (int i = 0; i < grid->n_nodes; ++i) {
    if (!fgets(line_buffer, sizeof(line_buffer), f_ptr)) {
      printf("Error during fgets()");
    }
    sscanf(line_buffer, "%*d, %lf, %lf", &grid->nodes[i].x, &grid->nodes[i].y);
  }

  if (!fgets(line_buffer, sizeof line_buffer, f_ptr)) {
    printf("Error during fgets()");
  }
  for (int i = 0; i < grid->n_elements; ++i) {
    if (!fgets(line_buffer, sizeof(line_buffer), f_ptr)) {
      printf("Error during fgets()");
    }
    sscanf(line_buffer, "%*d, %d, %d, %d, %d", &grid->elements[i].nodes[0],
           &grid->elements[i].nodes[1], &grid->elements[i].nodes[2],
           &grid->elements[i].nodes[3]);
  }

  if (!fgets(line_buffer, sizeof(line_buffer), f_ptr)) {
    printf("Error during fgets()");
  }
  if (!fgets(line_buffer, sizeof(line_buffer), f_ptr)) {
    printf("Error during fgets()");
  }

  for (int i = 0; i < grid->n_nodes; ++i) {
    grid->nodes[i].bc = false;
  }

  char* token = strtok(line_buffer, ",");
  while (token) {
    int id = atoi(token);
    grid->nodes[id-1].bc = true;
    token = strtok(NULL, ",");
  }

  fclose(f_ptr);

  glob_data->nip = NIP;
  grid->jacobians_calculated = 0;
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
  printf("- NIP                  :  %d\n\n", glob_data->nip);

  printf("===========================================\n\n");

  printf("                   NODES                   \n");
  printf("+-----------------------------------------+\n");
  for (int i = 0; i < grid->n_nodes; ++i) {
    printf("|%-7d|%16.10lf|%16.10lf|\n", i + 1, grid->nodes[i].x,
           grid->nodes[i].y);
  }
  printf("+-----------------------------------------+\n\n");

  printf("                  ELEMENTS                 \n");
  printf("+-----------------------------------------+\n");
  for (int i = 0; i < grid->n_elements; ++i) {
    printf("|%-9d|%7d|%7d|%7d|%7d|\n", i + 1, grid->elements[i].nodes[0],
           grid->elements[i].nodes[1], grid->elements[i].nodes[2],
           grid->elements[i].nodes[3]);
  }
  printf("+-----------------------------------------+\n");

  printf("                     BC                    \n");
  printf("+-----------------------------------------+\n");
  for (int i = 0; i < grid->n_nodes; ++i) {
    if (grid->nodes[i].bc) {
      printf("%d ", i+1);
    }
  }
  printf("\n");
  printf("+-----------------------------------------+\n");
}

void ExportJacobianData(const Grid* grid, const UniversalVals* uni_vals) {
  FILE* fptr;
  fptr = fopen("out/jacobian.txt", "w");

  fprintf(fptr,
          "==================== OBLICZONE STAŁE ====================\n\n");
  fprintf(fptr, "dN/dKsi\n");
  for (int i = 0; i < 4; ++i) {
    fprintf(fptr, "\t%12.8f %12.8f %12.8f %12.8f\n", uni_vals->dn_dksi[i][0],
            uni_vals->dn_dksi[i][1], uni_vals->dn_dksi[i][2],
            uni_vals->dn_dksi[i][3]);
  }
  fprintf(fptr, "dN/dEta\n");
  for (int i = 0; i < 4; ++i) {
    fprintf(fptr, "\t%12.8f %12.8f %12.8f %12.8f\n", uni_vals->dn_deta[i][0],
            uni_vals->dn_deta[i][1], uni_vals->dn_deta[i][2],
            uni_vals->dn_deta[i][3]);
  }

  fprintf(fptr,
          "\n\n============== WYNIKI OBLICZEŃ DLA PUNKTÓW ==============\n\n");

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

      fprintf(fptr, "\tMacierz Jakobiego dla punktu całkowania %d\n", ip + 1);

      fprintf(fptr, "\t\t%12.8lf  %12.8lf\n", j->j[0][0], j->j[0][1]);
      fprintf(fptr, "\t\t%12.8lf  %12.8lf\n", j->j[1][0], j->j[1][1]);

      fprintf(fptr, "\t\tdet[J] = %-12.8lf\n\n", j->det_j);

      fprintf(fptr, "\t\tWartość dN/dx\n");
      fprintf(fptr, "\t\t\t%12.5lf %12.5lf %12.5lf %12.5lf\n\n", j->dN_dx[0],
              j->dN_dx[1], j->dN_dx[2], j->dN_dx[3]);

      fprintf(fptr, "\t\tWartość dN/dy\n");
      fprintf(fptr, "\t\t\t%12.5lf %12.5lf %12.5lf %12.5lf\n\n", j->dN_dy[0],
              j->dN_dy[1], j->dN_dy[2], j->dN_dy[3]);
    }
  }
  fprintf(fptr, "\n");
  fclose(fptr);
}

void PrintGlobalH(const GlobHMatrix* h_matrix) {
  printf("\n\n");
  printf("                                                  H MATRIX\n");
  printf("+");
  for (int i = 0; i < h_matrix->n * 7 - 1; ++i) {
    printf("-");
  }
  printf("+\n");
  for (int i = 0; i < h_matrix->n; ++i) {
    printf("|");
    for (int j = 0; j < h_matrix->n; ++j) {
      printf("%6.3lf|", h_matrix->mat[i][j]);
    }
    printf("\n");
  }
  printf("+");
  for (int i = 0; i < h_matrix->n * 7 - 1; ++i) {
    printf("-");
  }
  printf("+\n");
}
