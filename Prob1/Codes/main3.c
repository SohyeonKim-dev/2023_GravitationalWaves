#include "linsys.h"
#include <stdlib.h>

void var_set_from_finer(Var *var, Var *finer_var)
{
    var->val[0] = finer_var->val[0];
    var->val[var->num_vals-1] = finer_var->val[finer_var->num_vals-1];
    int half_index = 0;

    for (int i = 0; i < finer_var->num_vals - 1; i++) {
        half_index = i / 2;
        var->val[half_index] = (finer_var->val[i-1]/4.) 
                            + (finer_var->val[i]/2.)
                            + (finer_var->val[i+1]/4.);

    } 
    // var->val[-1] = finer_var->val[-1];
}

int main()
{
    FILE *fp = fopen("input3.txt", "r");
    if (NULL == fp) raise(Error_msg_input_file_open);

    // read data
    Grid *finer_grid = new_grid_import(fp);
    Var *finer_var = new_var_import(finer_grid, true, true, fp);

    fclose(fp);

    Grid *coarser_grid = new_grid(finer_grid->num_cells / 2, finer_grid->a, finer_grid->b);
    Var *coarser_var = new_var_fill(coarser_grid, true, true, 0.);
    var_set_from_finer(coarser_var, finer_var);

    fp = fopen("output3.txt", "w");
    if (NULL == fp) raise(Error_msg_output_file_open);

    // output
    grid_export(coarser_grid, fp);
    var_export(coarser_var, fp);
    
    fclose(fp);

    free(finer_grid);
    free(finer_var);
    free(coarser_grid);
    free(coarser_var);
    
    return 0;
}