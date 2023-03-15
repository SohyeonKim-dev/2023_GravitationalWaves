#include "linsys.h"
#include <stdlib.h>

void var_set_from_coarser(Var *var, Var *coarser_var)
{
}

void var_set_from_finer(Var *var, Var *finer_var)
{
}

void system_relax(System *system, int num_steps)
{
    for (int step = 0; step < num_steps; step++) {
    }
}

void system_set_residual(System *system)
{
}

void multigrid_cycle(Multigrid *multigrid)
{
}

int main()
{
    FILE *fp = fopen("input5.txt", "r");
    if (NULL == fp) raise(Error_msg_input_file_open);

    // read data
    Grid *grid = new_grid_import(fp);
    
    fclose(fp);

    Var *sol = new_var_fill(grid, true, true, 1. + grid->a / grid->b);

    Var *src = new_var_fill(grid, true, true, 0);
    Multigrid *multigrid = new_multigrid(grid, sol, src, 1000, 1000, 1);

    double val;
    do {
        multigrid_cycle(multigrid);
    } while (system_get_log10_rms_residual_h2(multigrid->top) > -15.);
    
    // output
    fp = fopen("output5.txt", "w");
    if (NULL == fp) raise(Error_msg_output_file_open);
    grid_export(multigrid->top->grid, fp);
    var_export(multigrid->top->sol, fp);
    fclose(fp);

    del_multigrid(multigrid);

    return 0;
}