#include "linsys.h"
#include <stdlib.h>

void system_relax(System *system, int num_steps)
{
    for (int step = 0; step < num_steps; step++) {  // gaussian : 2.7772743701934814
        VAR_ITER(system->src, i) {
            system->sol->val[i] = (
                + system->sol->val[i - 1] 
                + system->sol->val[i + 1]
                - system->src->val[i] * (system->grid->h) * (system->grid->h) // gird -> h
            ) / 2. ;
        }
    }   
}

void system_set_residual(System *system) // ro_i 군요.......
{  
    VAR_ITER(system->res, i) {
        system->res->val[i] = 
            system->src->val[i] - (
                + system->sol->val[i-1]
                - system->sol->val[i] * 2.
                + system->sol->val[i + 1]
        ) / (system->grid->h * system->grid->h) ;
    }
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
    FILE *fp = fopen("input1.txt", "r");
    if (NULL == fp) raise(Error_msg_input_file_open);

    // read boundary values
    double u, v;
    if (EOF == fscanf(fp, "%lf %lf", &u, &v)) raise(Error_msg_input_file_format);

    // read data
    Grid *grid = new_grid_import(fp);
    Var *src = new_var_import(grid, false, false, fp);

    fclose(fp);

    double param[4] = {grid->a, grid->b, u, v};
    Var *sol = new_var_func(grid, true, true, initial_data, param);
    sol->val[0] = u;
    sol->val[grid->num_cells] = v;
    System *system = new_system(grid, sol, src);

    do {
        system_relax(system, 1000);
    } while (system_get_log10_rms_residual_h2(system) > -15.1);
    
    // output
    fp = fopen("output1.txt", "w");
    if (NULL == fp) raise(Error_msg_output_file_open);
    grid_export(system->grid, fp);
    var_export(system->sol, fp);
    fclose(fp);

    del_system(system);

    return 0;
}