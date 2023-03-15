#include "linsys.h"
#include <stdlib.h>

void var_set_from_coarser(Var *var, Var *coarser_var)
{
    //2번 함수
}

void var_set_from_finer(Var *var, Var *finer_var)
{
    //3번 함수
}

void system_relax(System *system, int num_steps)
{
    //1번 함수
    for (int step = 0; step < num_steps; step++) {
    }
}

void system_set_residual(System *system)
{
    //1번 함수
}
void multigrid_cycle(Multigrid *multigrid) // 얘 빼고 복붙
{
}

double initial_data(double x, double *param)
{
    double a = param[0];
    double b = param[1];
    double u = param[2];
    double v = param[3];
    return 0.;
}

int main()
{
    FILE *fp = fopen("input4.txt", "r");
    if (NULL == fp) raise(Error_msg_input_file_open);

    // read boundary values
    double u, v;
    if (EOF == fscanf(fp, "%lf %lf", &u, &v)) raise(Error_msg_input_file_format);

    // read data
    Grid *grid = new_grid_import(fp);
    if (grid->num_cells <= 0 || grid->num_cells & (grid->num_cells - 1) != 0) raise (Error_msg_num_cells);
    Var *src = new_var_import(grid, false, false, fp);

    fclose(fp);

    double param[4] = {grid->a, grid->b, u, v};
    Var *sol = new_var_func(grid, true, true, initial_data, param);
    sol->val[0] = u;
    sol->val[grid->num_cells] = v;
    Multigrid *multigrid = new_multigrid(grid, sol, src, 1, 1, 1);

    double val;
    do {
        multigrid_cycle(multigrid);
    } while (system_get_log10_rms_residual_h2(multigrid->top) > -15.);
    
    // output
    fp = fopen("output4.txt", "w");
    if (NULL == fp) raise(Error_msg_output_file_open);
    grid_export(multigrid->top->grid, fp);
    var_export(multigrid->top->sol, fp);
    fclose(fp);

    del_multigrid(multigrid);

    return 0;
}