#include <stdio.h> /* needed for input, output, ... */
#include <stdlib.h> /* needed for EXIT_SUCCESS, ... */
#include <math.h> /* nedded for myFloor() and myCeil() */
#include <glpk.h> /* the linear programming toolkit */

#define EPSILON 0.00001 /* small value to deal with rounding issues */

/* global variables -- YOU SHOULD NOT CHANGE THESE! */
/* (You are allowed to add your own if you want.) */
int size; /* number of rows/columns in the input */
double *input; /* array that contains the input */
double *solution; /* array that contains the solution */
int debug; /* flag for debug mode; 1 means debug mode, 0 means debug off */

/* prototypes of functions -- YOU SHOULD NOT CHANGE THESE! */
/* (Feel free to add your own as you like.) */
int readInput(char *filename); /* reads graph from file */
/* readInput creates/fills the data structures (global variables) as needed */
/* it returns 0 if all is okay and 1 otherwise */
int computeSolution(void); /* computes the solution, stores it in input */
/* the return value is 1 if a solution exists, 0 otherwise */
void checkSolution(void); /* checks if a solution is valid */
double myFloor(double x); /* slightly 'hacky' myFloor function */
double myCeil(double x); /* slightly 'hacky' myCeil function */
void printMatrix(char *title, double *matrix); /* print a matrix */
void fillRValues(double *, int);
void fillCValues(double *, int);

void fillRValues(double* rValues, int numR){
  int lastJPos = 0;
  int newSize = size;
  for (int i = 0; i< numR; ++i){
    double totalOfRows = 0;
    for (int j = lastJPos; j<newSize; ++j){
      totalOfRows += input[j];
      lastJPos = j+1;
    }
    newSize = newSize + lastJPos;
    rValues[i] = totalOfRows;
  }
}

void fillCValues(double *cValues, int numC){
  int lastJPos = 0;
  int newSize = size;
  for (int i = 0; i<size; ++i){
    double totalOfRows = 0;
    for (int j = 0; j<numC; ++j){
      int index = j*size+i;
      totalOfRows += input[index];
      lastJPos = j+1;
    }
    newSize = newSize + lastJPos;
    cValues[i] = totalOfRows;
  }
}

void fillIndexArray(int* index) {
  for (int i = 1; i < size + 1; ++i) {
    index[i] = i;
  }
}

/* This is the function that actually solves the problem. */
/* It is essentially empty and not functional. */
/* Your own implementation needs to go in here. */
int computeSolution(void) {
  int i, j; // For counting in for loops
  const int numX = size * size; // Number of X_nm variables
  const int numC = size; // Number of C_m variables
  const int numR = size; // Number of R_n variables
  const int numE = size * size; // Number of e_t variables
  const int numVariables = numC + numR + 2*numX + numE; // this is the number of columns in the matrix for GLPK. 2*numX because of Xnm' and Xnm'' values.
  double cPrimeValues[numC];
  double rPrimeValues[numR];
  glp_prob *lp;
  int index[numVariables + 1];
  memset(index, 0, sizeof(index));
  double row[numVariables + 1];
  memset(row, 0.0, sizeof(row));

  lp = glp_create_prob();
  glp_set_obj_dir(lp, GLP_MAX);

  glp_add_cols(lp, numVariables);

  fillRValues(rPrimeValues, numR);
  fillCValues(cPrimeValues, numC);

/* ADD COLUMN BOUNDS STARTS HERE *********************************************/
  // Add C_m' column variables
  int maxLoopSize = numC;
  for (i=1; i<=maxLoopSize; i++ ) {
    //*2 as C_m' max and min flow is 2 times greater than C_m because the outgoing of C_m is greater by 2
    const double lowerBound = 2*myFloor(cPrimeValues[i-1]);
    const double upperBound = 2*myCeil(cPrimeValues[i-1]);
    if (lowerBound != upperBound) {
      glp_set_col_bnds(lp, i, GLP_DB, lowerBound, upperBound);
    } else {
      glp_set_col_bnds(lp, i, GLP_FX, lowerBound, upperBound);
    }
  }

  // Add R_m' column variables
  maxLoopSize += numR;
  for (j = 0, i = numC; i<= maxLoopSize; ++j ,++i){
    //*2 as R_n' max and min flow is 2 times greater than max and min of R_n because the outgoing of node R_n is greater by a magnitude of 2 than R_n
    const double lowerBound = 2*myFloor(rPrimeValues[j]);
    const double upperBound = 2*myCeil(rPrimeValues[j]);
    if (lowerBound != upperBound) {
      glp_set_col_bnds(lp, i, GLP_DB, lowerBound, upperBound);
    } else {
      glp_set_col_bnds(lp, i, GLP_FX, lowerBound, upperBound);
    }
  }

  // Add X_nm' column variables numC and numR along in the column.
  maxLoopSize += numX;
  for (j = 0, i = numC + numR + 1; i<=maxLoopSize; ++j, ++i){
    const double lowerBound = myFloor(input[j]);
    const double upperBound = myCeil(input[j]);
    if (lowerBound != upperBound) {
      glp_set_col_bnds(lp, i, GLP_DB, lowerBound, upperBound);
    } else {
      glp_set_col_bnds(lp, i, GLP_FX, lowerBound, upperBound);
    }
  }

  // Add X_nm'' column variables numC, numR, and numX along in the column.
  maxLoopSize += numX;
  for (j = 0, i = numC + numR + numX + 1; i<=maxLoopSize; ++j, ++i){
    const double lowerBound = myFloor(input[j]);
    const double upperBound = myCeil(input[j]);
    if (lowerBound != upperBound) {
      glp_set_col_bnds(lp, i, GLP_DB, lowerBound, upperBound);
    } else {
      glp_set_col_bnds(lp, i, GLP_FX, lowerBound, upperBound);
    }
  }

  // Add e_t column variables numC, numR, and 2numX along in the column
  maxLoopSize += numE;
  for (j = 0, i = numC + numR + 2*numX + 1; i<=maxLoopSize; ++j, ++i){
    const double lowerBound = myFloor(2*input[j]);
    const double upperBound = myCeil(2*input[j]);
    if (lowerBound != upperBound) {
      glp_set_col_bnds(lp, i, GLP_DB, lowerBound, upperBound);
    } else {
      glp_set_col_bnds(lp, i, GLP_FX, lowerBound, upperBound);
    }
  }
/* ADD COLUMN BOUNDS ENDS HERE *********************************************/

/* ADD ROW COEFFICIENTS STARTS HERE *******************************************************/
  // Add C_m coefficients and R_n: C_1' + C_2' + .... + C_m' + R_1 + R_2 + .... + R_n'
  maxLoopSize = numC + numR;
  for (i=1; i<= maxLoopSize; i++) {
    glp_set_obj_coef(lp, i, 1.0);
  }

/* ADD ROW COEFFICIENTS ENDS HERE *******************************************************/

/* ADD ROWS BOUNDS STARTS HERE *******************************************************/

  const int numRows = numC + numR + numX + numE;
  glp_add_rows(lp, numRows);

  // Set bounds for C_m
  maxLoopSize = numC;
  for (i = 1; i<=maxLoopSize; ++i){
    const double lowerBound = 2*myFloor(cPrimeValues[i-1]);
    const double upperBound = 2*myCeil(cPrimeValues[i-1]);
    if (lowerBound != upperBound) {
      glp_set_row_bnds(lp, i, GLP_DB, lowerBound, upperBound);
    } else {
      glp_set_row_bnds(lp, i, GLP_FX, lowerBound, upperBound);
    }
  }

  // Set bounds for R_n
  maxLoopSize += numR;
  for (j=0, i = numC + 1; i<=maxLoopSize; ++i, ++j){
    const double lowerBound = 2*myFloor(rPrimeValues[j]);
    const double upperBound = 2*myCeil(rPrimeValues[j]);
    if (lowerBound != upperBound) {
      glp_set_row_bnds(lp, i, GLP_DB, lowerBound, upperBound);
    } else {
      glp_set_row_bnds(lp, i, GLP_FX, lowerBound, upperBound);
    }
  }

  // Set bounds for X_{n,m}
  maxLoopSize += numX;
  for (j=0, i=numC + numR + 1; i<= maxLoopSize; ++i, ++j){
    const double lowerBound = 2*myFloor(input[j]);
    const double upperBound = 2*myCeil(input[j]);
    if (lowerBound != upperBound) {
      glp_set_row_bnds(lp, i, GLP_DB, lowerBound, upperBound);
    } else {
      glp_set_row_bnds(lp, i, GLP_FX, lowerBound, upperBound);
    }
  }

  // Set bounds for Z_t = 0
  maxLoopSize += numE;
  for (i=numC + numR + numX + 1; i<=maxLoopSize; ++i){
    glp_set_row_bnds(lp, i, GLP_FX, 0.0, 0.0);
  }

/* ADD ROWS BOUNDS ENDS HERE *******************************************************/

/* ADD THE DATA FOR EACH ROW STARTS HERE ************************************************/

  // Set C_m constraints
  maxLoopSize = numC;
  for(int c = 0, i = 1; i<=maxLoopSize; ++i, c+= 3){
    index[1] = i, row[i] = 1;
    // maxSize is 1 + size because there will be 1 C_m' and size number of x_nm' variables
    for (j = 2; j<=1 + size; ++j){
      const int indexValue = numC + numR + j -1 + c;
      index[j] = indexValue, row[indexValue] = -1;
    }
    glp_set_mat_row(lp, i, 1 + size, index, row);
  }

  // Set R_n constraints
  maxLoopSize = numC + numR;
  for (int ii = 1, i = numC + 1; i<= maxLoopSize; i++, ii++) {
    index[1] = i, row[i] = 1;
    for (int c = 0, k = 2, j = 1 + numC; j<= numC + size; ++j, ++k, c+=3){
      const int indexValue = numC + numR + numX + c + ii;
      index[k] = indexValue, row[indexValue] = -1;
    }
    glp_set_mat_row(lp, i, 1 + size, index, row);
  }

  // Set X_{n,m} constraints
  maxLoopSize = numC + numR + numX;
  for (int i = numC + numR + 1; i <= maxLoopSize; i++){
    int indexValue = i;
    index[1] = indexValue, row[indexValue] = 1;
    indexValue = i + numX;
    index[2] = indexValue, row[indexValue] = 1;
    indexValue = i + numX + numX;
    index[3] = indexValue, row[indexValue] = -1;
    glp_set_mat_row(lp, i, 3, index, row);
  }

  // Set Z_t constraints
  maxLoopSize = numC + numR + numX + numE;
  for (int i = 1; i<=numE; ++i){
    int indexValue = i + numC + numR;
    index[1] = indexValue, row[indexValue] = 1;
    indexValue = i + numC + numR + numX;
    index[2] = indexValue, row[indexValue] = -1;
    glp_set_mat_row(lp, i, 2, index, row);
  }

/* ADD THE DATA FOR EACH ROW ENDS HERE **************************************************/

  glp_term_out(0); // switch off debug output from GLPK 
  glp_simplex(lp, NULL);

  // Fill in solution
  double perfectValue = glp_get_obj_val(lp);
  for (j = 0, i = numC + numR + 1; i<=numC + numR + numX; ++i, ++j) {
      solution[j] = glp_get_col_prim(lp, i);
  }

  // glp_delete_prob(lp); // release memory used for LP

  return 1; /* This is not always correct, of course, and needs to be changed. */
}

/* YOU SHOULD NOT CHANGE ANYTHING BELOW THIS LINE! */

int main(int argc, char **argv) {
    int arg; /* used to run over the command line parameters */

    if ( argc<2 ) { /* no command line parameter given */
        fprintf(stderr, "Usage: %s [file1] [file2] [file3] [...]\n"
                "Where each [file] indicates the name of a file with an input.\n",
                argv[0]);
        exit(EXIT_FAILURE);
    }

    if ( argv[1][0]=='-' && argv[1][1]=='d' && argv[1][2]==0 ) {
        /* If the first parameter is -d we activate debug mode. */
        debug=1; /* switch debug mode on */
        fprintf(stdout, "DEBUG: Debug mode activated\n"); /* tell us */
    } else {
        debug=0; /* switch debug mode off */
    }

    for ( arg=1+debug; arg<argc; arg++ ) { /* go over remaining parameters */
        if ( readInput(argv[arg]) ) { /* try to read file */
            /* readInput returned with error message */
            fprintf(stderr, "%s: Cannot read input with filename %s. Skipping.\n",
                    argv[0], argv[arg]);
        } else {
            if ( computeSolution() ) {
                fprintf(stdout, "%s: Input %s has a solution.\n",
                        argv[0], argv[arg]);
                checkSolution();
                printMatrix("Input", input);
                printMatrix("Solution", solution);
            } else {
                fprintf(stdout, "%s: Input %s does not have a solution.\n",
                        argv[0], argv[arg]);
            }
            /* free memory for next input */
            free(input);
            free(solution);
        }
    }
    return EXIT_SUCCESS;
}

/* The following function prints a matrix including the row and column sums */
void printMatrix(char *title, double *matrix) {
  int	i, j; /* looping over rows and columns */
  double sum; /* to compute the sum */

  fprintf(stdout, "%s:\n", title);
  /* print the matrix and compute the row sums on the fly */
  for ( i=0; i<size; i++ ) {
    for ( j=0, sum=0.0; j<size; j++ ) {
      sum += matrix[i*size+j];
      fprintf(stdout,"%8.2f ", matrix[i*size+j]);
    }
  fprintf(stdout, "(row sum: %8.2f)\n", sum);
  }
  /* print separating line */
  for ( j=0; j<size; j++ ) {
    fprintf(stdout, "---------");
  }
  fprintf(stdout,"\n");
  /* now compute the column sums and print them */
  for ( j=0; j<size; j++ ) { /* we consistently use j for columns */
    for ( i=0, sum=0.0; i<size; i++ ) {
      sum += matrix[i*size+j];
    }
    fprintf(stdout,"%8.2f ", sum);
  }
  fprintf(stdout,"(column sums)\n");
}

/* The following function checks of a solution is valid. */
void checkSolution(void) {
    int	i, j; /* to run over the arrays */
    double sum1, sum2; /* to compute the sums over the rows and columns */

  /* check rows and that all numbers are integers and rounded */
    for ( i=0; i<size; i++ ) {
      for ( j=0, sum1=sum2=0.0; j<size; j++ ) {
        sum1 += input[i*size+j];
        sum2 += solution[i*size+j];
        if ( myFloor(solution[i*size+j]) != solution[i*size+j] ) {
          fprintf(stdout, "Error: %lf is not an integer (%d/%d).\n",
          solution[i*size+j], i, j);
          return;
        }
        if ( (myFloor(input[i*size+j])!=solution[i*size+j]) &&
            (myCeil(input[i*size+j])!=solution[i*size+j]) ) {
          fprintf(stdout, "Error: %lf is not rounded from %lf (%d/%d).\n",
          solution[i*size+j], input[i*size+j], i, j);
          return;
        }
      }
      if ( (myFloor(sum1)!=sum2) && (myCeil(sum1)!=sum2) ) {
        fprintf(stdout, "Error: Row sum for row %d not valid.\n", i);
        return;
      }
    }
    /* check columns*/
    for ( j=0; j<size; j++ ) {
      for ( i=0, sum1=sum2=0.0; i<size; i++ ) {
        sum1 += input[i*size+j];
        sum2 += solution[i*size+j];
      }
      if ( (myFloor(sum1)!=sum2) && (myCeil(sum1)!=sum2) ) {
        fprintf(stdout, "Error: Row sum for row %d not valid.\n", i);
        return;
      }
    }
}

int readInput(char *filename) {
    FILE *fh; /* file handle to read input */
    int i, j; /* variables to run over array */
    double value; /* value read from input file */

    /* try to open the file for reading */
    if ( ( fh = fopen(filename, "rt") ) == NULL ) {
        if ( debug ) {
            fprintf(stdout, "DEBUG: Unable to open file %s for reading.\n",
             filename);
        }
        return 1; /* unable to open file, flag failure */
    }
    /* read the first integer, the number of columns/rows */
    if ( fscanf(fh, "%d", &size)!= 1) {
        if ( debug ) {
            fprintf(stdout, "DEBUG: Unable to read input size.\n");
        }
        fclose(fh); /* close file to avoid ununsed open files */
        return 1; /* flag failure */
    }

    if ( size<2 ) {
        if ( debug ) {
            fprintf(stdout, "DEBUG: Received %d as input size.\n", size);
        }
        fclose(fh); /* close file to avoid unused open files */
        return 1; /* flag failure */
    }
    /* allocate the memory for the input */
    if ( ( input = (double *)malloc(sizeof(double)*size*size) ) == NULL ) {
        if ( debug ) {
            fprintf(stdout, "DEBUG: Unable to allocate %ld bytes.\n",
                    sizeof(int)*size*size);
        }
        fclose(fh); /* close file to avoid unused open files */
        return 1; /* flag failure */
    }
  /* allocate the memory for the solution */
    if ( ( solution = (double *)malloc(sizeof(double)*size*size) ) == NULL ) {
        if ( debug ) {
            fprintf(stdout, "DEBUG: Unable to allocate %ld bytes.\n",
                    sizeof(int)*size*size);
        }
        fclose(fh); /* close file to avoid unused open files */
        return 1; /* flag failure */
    }

    /* read the actual values */
    for ( i=0; i<size; i++ ) {
        for ( j=0; j<size; j++ ) {
            /* attempt to read next value */
            if ( fscanf(fh, "%lf", &value)!= 1) {
                if ( debug ) {
                    fprintf(stdout, "DEBUG: Unable to read input value (%d/%d).\n",
                            i, j);
                }
                /* free memory and close file */
                free(input);
                free(solution);
                fclose(fh);
                return 1; /* flag failure */
            }
            input[i*size+j]=value;
        }
    }
  if ( debug ) {
    fprintf(stdout,"Read the following input with size %d x %d\n", size, size);
    for ( i=0; i<size; i++ ) {
      for ( j=0; j<size; j++ ) {
        fprintf(stdout,"%5.2f ", input[i*size+j]);
      }
      fprintf(stdout,"\n");
    }
  }
  return 0; /* all okay */
}

double myFloor(double x) {
  return floor(x+EPSILON);
}

double myCeil(double x) {
  return ceil(x-EPSILON);
}
