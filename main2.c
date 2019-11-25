#include <stdio.h> /* needed for input, output, ... */
#include <stdlib.h> /* needed for EXIT_SUCCESS, ... */
#include <glpk.h> /* the linear programming toolkit */

int main(void) {
    glp_prob *lp; /* the linear programming problem */
    int index[5]; /* indices to define constraint coefficients */
    double row[5]; /* values to define constraint coefficients */
    int i; /* loop variable to run over LP variables */

    lp = glp_create_prob(); /* create LP problem instance */
    glp_set_obj_dir(lp, GLP_MAX); /* set maximisation as objective */

    glp_add_cols(lp, 5); /* set number of variables to 5 */
    /* set bounds for the 5 variables */
    for ( i=1; i<=5; i++ ) { /* counting starts at 1, not 0! */
        glp_set_col_bnds(lp, i, GLP_LO, 0.0, 0.0); /* lower bound 0.0 set */
    }
    /* set objective function 2x_1+x_2(+0x_3+0x_4+0_x5) */
    glp_set_obj_coef(lp, 1, 2.0); /* 2x_1 */
    glp_set_obj_coef(lp, 2, 1.0); /* 1x_2 */
    glp_set_obj_coef(lp, 3, 0.0); /* 0x_3 */
    glp_set_obj_coef(lp, 4, 0.0); /* 0x_4 */
    glp_set_obj_coef(lp, 5, 0.0); /* 0x_5 */

    

    glp_add_rows(lp, 3); /* set number of constraints to 3 (3 constraints) */
    /* set constraint 1.5x_1-1x_2+1x_3+0x_4+0x_5=0 */
    glp_set_row_bnds(lp, 1, GLP_FX, 0.0, 0.0); /* set right-hand side to 0 */
    index[1]=1; row[1] = 1.5; /* 1.5x_1 */
    index[2]=2; row[2] = -1.0; /* -1x_2 */
    index[3]=3; row[3] = 1.0; /* 1x_3 */
    glp_set_mat_row(lp, 1, 3, index, row);
    /* set constraint 1x_1+1x_2+0x_3+1x_4+0x_5=4680 */
    glp_set_row_bnds(lp, 2, GLP_FX, 4680.0, 0.0); /* set right-hand side to 4680 */
    index[1]=1; row[1] = 1.0; /* 1x_1 */
    index[2]=2; row[2] = 1.0; /* 1x_2 */
    index[3]=4; row[3] = 1.0; /* 1x_4 */
    glp_set_mat_row(lp, 2, 3, index, row);
    /* set constraint 0.2x_1+1x_2+0x_3+0x_4+1x_5=3900 */
    glp_set_row_bnds(lp, 3, GLP_FX, 3900.0, 0.0); /* set right-hand side to 3900 */
    index[1]=1; row[1] = 0.2; /* 0.2x_1 */
    index[2]=2; row[2] = 1.0; /* 1x_2 */
    index[3]=5; row[3] = 1.0; /* 1x_5 */
    glp_set_mat_row(lp, 3, 3, index, row);

    /* now solve LP */
    glp_term_out(0); /* switch off debug output from GLPK */
    glp_simplex(lp, NULL); /* solve the LP */
    printf("Optimal solution has value %f (%f minutes reading, %f minutes programming).\n",
           glp_get_obj_val(lp), glp_get_col_prim(lp, 1), glp_get_col_prim(lp, 2));

    glp_delete_prob(lp); /* release memory used for LP */
    return EXIT_SUCCESS; /* signal all went okay */
}