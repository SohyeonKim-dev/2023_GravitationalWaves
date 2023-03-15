#include "linsys.h"
#include <stdlib.h>

void var_set_from_coarser(Var *var, Var *coarser_var)
{
    var->val[0] = coarser_var->val[0];
    //var->val[var->num_vals-1] = coarser_var->val[coarser_var->num_vals];
    var->val[-1] = coarser_var->val[-1];
    // int num = coarser_var->num_vals;
    int num = var->num_vals;
    int count=0;

    for (int i = 2; i < num; i++) {
        if (i % 2 == 0) {
            count++;
            // var->grid[i] = coarser_var->grid[i];
            var->val[i] = coarser_var->val[count];
        }
    }

    
    for (int j = 1; j < num; j++) {
        if (j % 2 == 1) {
            var->val[j] = (var->val[j-1] + var->val[j+1]) / 2;
        }
    }
}

int main()
{
    FILE *fp = fopen("input2.txt", "r");
    if (NULL == fp) raise(Error_msg_input_file_open);

    // read data
    Grid *coarser_grid = new_grid_import(fp);
    Var *coarser_var = new_var_import(coarser_grid, true, true, fp);

    fclose(fp);

    Grid *finer_grid = new_grid(coarser_grid->num_cells * 2, coarser_grid->a, coarser_grid->b);
    Var *finer_var = new_var_fill(finer_grid, true, true, 0.);
    var_set_from_coarser(finer_var, coarser_var);

    fp = fopen("output2.txt", "w");
    if (NULL == fp) raise(Error_msg_output_file_open);

    // output
    grid_export(finer_grid, fp);
    var_export(finer_var, fp);

    fclose(fp);

    free(coarser_grid);
    free(coarser_var);
    free(finer_grid);
    free(finer_var);

    return 0;
}