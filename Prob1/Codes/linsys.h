#ifndef LINSYS_H
#define LINSYS_H

#include <stdio.h>
#include <stdbool.h>

// To exit after print error message

void raise(char *msg);

#define Error_msg_input_file_format "The input file format invalid.\n"
#define Error_msg_input_file_open "Error on opening the input file.\n"
#define Error_msg_output_file_open "Error on opening the output file.\n"
#define Error_msg_num_cells "num_cells should be power of 2\n"

// Grid structure

typedef struct _Grid {
    int num_cells;
    double a;
    double b;
    double *x;
    double h;
} Grid;

Grid* new_grid(int num_cells, double a, double b);
void del_grid(Grid *grid_c);
Grid* new_grid_import(FILE *fp);
void grid_export(Grid *grid, FILE *fp);

// Var structure

typedef struct _Var {
    Grid *grid;
    double *val;
    bool left_closed;
    bool right_closed;
    int num_vals;
} Var;

void del_var(Var *var_c);
Var *new_var_move(Grid *grid, bool left_closed, bool right_closed, Var *other_c);
Var* new_var_func(Grid *grid, bool left_closed, bool right_closed, double func(double, double*), double *param);
Var* new_var_fill(Grid *grid, bool left_closed, bool right_closed, double val);
void var_fill(Var *var, double val);
Var* new_var_import(Grid *grid, bool left_closed, bool right_closed, FILE *fp);
void var_set_from_coarser(Var *var, Var *coarser_var);
void var_set_from_finer(Var *var, Var *finer_var);
void var_add(Var *var, Var *operand);
double var_get_rms(Var *var);
void var_export(Var *var, FILE *fp);

#define VAR_ITER(var, i) for (int i = (var->left_closed ? 0 : 1); i < var->grid->num_cells + (var->right_closed ? 1 : 0); i++)

// Linear System structure

typedef struct _System {
    Grid *grid;
    Var *sol; // numerical solution
    Var *src; // source term
    Var *res; // residual
    Var *temp;
} System;

System* new_system(Grid *grid_c, Var *sol_c, Var *src_c);
void del_system(System *system_c);
void system_relax(System *system, int num_steps);
void system_set_residual(System *system);
double system_get_log10_rms_residual_h2(System *system);

// Multigrid structure

typedef struct _Multigrid {
    int num_lvs;
    System **system;
    System *top;
    int num_pre_iters;
    int num_post_iters;
    int order;
} Multigrid;

Multigrid* new_multigrid(Grid *grid_c, Var *sol_c, Var *src_c, int num_pre_iters, int num_post_iters, int order);
void del_multigrid(Multigrid *multigrid_c);
void multigrid_cycle(Multigrid *multigrid);

#endif