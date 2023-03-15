#include <stdlib.h>
#include <math.h>
#include "linsys.h"

void raise(char *msg)
{
    fprintf(stderr, "%s", msg);
    exit(EXIT_FAILURE);
}

// Grid structure

Grid* new_grid(int num_cells, double a, double b)
{
    Grid *grid = malloc(sizeof(Grid));
    grid->num_cells = num_cells;
    grid->a = a;
    grid->b = b;
    grid->h = (b - a) / (double)grid->num_cells;    
    grid->x = malloc((num_cells + 1) * sizeof(double));    
    grid->x[0] = a;
    for (int i = 1; i < grid->num_cells; i++) {
        grid->x[i] = a + i * grid->h;
    }
    grid->x[grid->num_cells] = b;
    return grid;
}

void del_grid(Grid *grid_c)
{
    free(grid_c->x);
    free(grid_c);
}

Grid* new_grid_import(FILE *fp)
{
    int num_cells;
    double a, b;
    if (EOF == fscanf(fp, "%d %lf %lf", &num_cells, &a, &b)) raise(Error_msg_input_file_format);
    return new_grid(num_cells, a, b);
}

void grid_export(Grid *grid, FILE *fp)
{
    fprintf(fp, "\n%d %.16e %.16e", grid->num_cells, grid->a, grid->b);
}

// Var structure

Var* _new_var(Grid *grid, bool left_closed, bool right_closed)
{
    Var *var = malloc(sizeof(Var));
    var->grid = grid;
    var->left_closed = left_closed;
    var->right_closed = right_closed;
    var->num_vals = var->grid->num_cells - 1 + (left_closed ? 1 : 0) + (right_closed ? 1 : 0);
    return var;
}

Var* _new_var_alloc(Grid *grid, bool left_closed, bool right_closed)
{
    Var *var = _new_var(grid, left_closed, right_closed);
    var->val = malloc((grid->num_cells + 1) * sizeof(double));
    return var;
}

void del_var(Var *var_c)
{
    free(var_c->val);
    free(var_c);
}

Var *new_var_move(Grid *grid, bool left_closed, bool right_closed, Var *other_c)
{
    Var *var = _new_var(grid, left_closed, right_closed);
    var->val = other_c->val;
    // consume other
    free(other_c);
    return var;
}

Var* new_var_func(Grid *grid, bool left_closed, bool right_closed, double func(double, double*), double *param)
{
    Var *var = _new_var_alloc(grid, left_closed, right_closed);
    VAR_ITER(var, i) {
        var->val[i] = func(grid->x[i], param);
    }
    return var;
}

Var* new_var_fill(Grid *grid, bool left_closed, bool right_closed, double val)
{
    Var *var = _new_var_alloc(grid, left_closed, right_closed);
    var_fill(var, val);
    return var;
}

void var_fill(Var *var, double val)
{
    VAR_ITER(var, i) {
        var->val[i] = val;
    }
}

Var* new_var_import(Grid *grid, bool left_closed, bool right_closed, FILE *fp)
{
    Var *var = _new_var_alloc(grid, left_closed, right_closed);
    VAR_ITER(var, i) {
        if (EOF == fscanf(fp, "%lf", var->val + i)) raise(Error_msg_input_file_format);
    }
    return var;
}

void var_add(Var *var, Var *operand)
{
    VAR_ITER(var, i) {
        var->val[i] += operand->val[i];
    }
}

double var_get_rms(Var *var)
{
    double sum = 0.;
    VAR_ITER(var, i) {
        sum += var->val[i] * var->val[i];
    }
    return sqrt(sum / var->num_vals);
}

void var_export(Var *var, FILE *fp)
{
    fprintf(fp, "\n");
    VAR_ITER(var, i) {
        fprintf(fp, "\t%.16e", var->val[i]);
    }
}

// Linear System structure

System* new_system(Grid *grid_c, Var *sol_c, Var *src_c)
{
    System *system = malloc(sizeof(System));
    system->grid = grid_c;
    system->sol = new_var_move(system->grid, true, true, sol_c);
    system->src = new_var_move(system->grid, false, false, src_c);
    system->res = new_var_fill(system->grid, false, false, 0.);
    system->temp = new_var_fill(system->grid, true, true, 0.);
    return system;
}

void del_system(System *system_c)
{
    del_var(system_c->sol);
    del_var(system_c->src);
    del_var(system_c->res);
    del_var(system_c->temp);
    del_grid(system_c->grid);
    free(system_c);
}

double system_get_log10_rms_residual_h2(System *system)
{
    system_set_residual(system);
    return log10(var_get_rms(system->res) * system->grid->h * system->grid->h);
}

// Multigrid structure

Multigrid* new_multigrid(Grid *grid_c, Var *sol_c, Var *src_c, int num_pre_iters, int num_post_iters, int order)
{
    Multigrid *multigrid = malloc(sizeof(Multigrid));
    int num_lvs = (int)log2(grid_c->num_cells);
    multigrid->num_lvs = num_lvs;
    multigrid->system = malloc(num_lvs * sizeof(System*));
    for (int i = 0; i < num_lvs - 1; i++) {
        Grid *sub_grid = new_grid(2 << i, grid_c->a, grid_c->b);
        Var *sub_sol = new_var_fill(sub_grid, true, true, 0.);
        Var *sub_src = new_var_fill(sub_grid, false, false, 0.);
        multigrid->system[i] = new_system(sub_grid, sub_sol, sub_src);
    }
    multigrid->system[num_lvs - 1] = new_system(grid_c, sol_c, src_c);
    multigrid->top = multigrid->system[num_lvs - 1];
    multigrid->num_pre_iters = num_pre_iters;
    multigrid->num_post_iters = num_post_iters;
    multigrid->order = order;

    return multigrid;
}

void del_multigrid(Multigrid *multigrid_c)
{
    for (int i = 0; i < multigrid_c->num_lvs; i++)
        del_system(multigrid_c->system[i]);
    free(multigrid_c->system);
    free(multigrid_c);
}