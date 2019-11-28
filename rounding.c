#include <stdio.h> /* needed for input, output, ... */
#include <stdlib.h> /* needed for EXIT_SUCCESS, ... */
#include <math.h> /* nedded for myFloor() and myCeil() */
#include <glpk.h> /* the linear programming toolkit */
#include <string.h>

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
  long newSize = size;
  long j = 0, i;
  for (i = 0; i < numR; ++i){
    double totalOfRows = 0;
    for (; j<newSize; ++j){
      totalOfRows += input[j];
    }
    newSize = newSize + numR;
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

double sumArray(double arr[], int n)  {  
    double sum = 0;  
    for (int i = 0; i < n; i++){
      sum += arr[i];  
    }
    return sum;  
}

/* This is the function that actually solves the problem. */
/* It is essentially empty and not functional. */
/* Your own implementation needs to go in here. */
int computeSolution(void) {
  int i, j; // For counting in for loops
  const int numX = size * size; // Number of X_nm variables
  const int numC = size; // Number of C_m variables
  const int numR = size; // Number of R_n variables
  const int numVariables = numC + numR + numX;
  double cValues[numC];
  double rValues[numR];
  glp_prob *lp;
  int index[numVariables + 1];
  memset(index, 0, sizeof(index));
  double row[numVariables + 1];
  memset(row, 0.0, sizeof(row));

  lp = glp_create_prob();
  glp_set_obj_dir(lp, GLP_MAX);

  glp_add_cols(lp, numVariables);

  fillRValues(rValues, numR);
  fillCValues(cValues, numC);

/* ADD COLUMN BOUNDS STARTS HERE *********************************************/
  // Add C_m' column variables
  int maxLoopSize = numC;
  for (i=1; i<=maxLoopSize; i++ ) {
    const double lowerBound = myFloor(cValues[i-1]);
    const double upperBound = myCeil(cValues[i-1]);
    if (lowerBound != upperBound) {
      glp_set_col_bnds(lp, i, GLP_DB, lowerBound, upperBound);
    } else {
      glp_set_col_bnds(lp, i, GLP_FX, lowerBound, upperBound);
    }
  }

  // Add R_m' column variables
  maxLoopSize += numR;
  for (j = 0, i = numC+1; i<= maxLoopSize; ++j ,++i){
    const double lowerBound = myFloor(rValues[j]);
    const double upperBound = myCeil(rValues[j]);
    if (lowerBound != upperBound) {
      glp_set_col_bnds(lp, i, GLP_DB, lowerBound, upperBound);
    } else {
      glp_set_col_bnds(lp, i, GLP_FX, lowerBound, upperBound);
    }
  }

  // Add X_nm column variables numC and numR along in the column.
  maxLoopSize = numC + numR + numX;
  for (j = 0, i = numC + numR + 1; i<=maxLoopSize; ++j, ++i){
    const double lowerBound = myFloor(input[j]);
    const double upperBound = myCeil(input[j]);
    if (lowerBound != upperBound) {
      glp_set_col_bnds(lp, i, GLP_DB, lowerBound, upperBound);
    } else {
      glp_set_col_bnds(lp, i, GLP_FX, lowerBound, upperBound);
    }
  }

  // Add the co-efficient function C_1' + C_2' + ... C_m'
  maxLoopSize = numC + numR + numX;
  for (i=1; i<= maxLoopSize; i++) {
    glp_set_obj_coef(lp, i, 0.0);
  }

  maxLoopSize = numC;
  for (i=1; i<= maxLoopSize; i++) {
    glp_set_obj_coef(lp, i, 1.0);
  }

  const int numRows = numC + numR;
  glp_add_rows(lp, numRows);

  // Set bounds of the rows to 0
  maxLoopSize = numRows;
  for (i = 1; i<=maxLoopSize; ++i){
    glp_set_row_bnds(lp, i, GLP_FX, 0.0, 0.0);
  }

  // Set C_m constraints
  maxLoopSize = numC;
  for(int ii = 1, i = 1; i<=maxLoopSize; ++i, ++ii){
    index[1] = i, row[1] = 1.0;
    for (int c = 0, k = 2 ,j = 1 + numC; j<=numC + size; ++j, ++k, c+=size){
      const int indexValue = numC + numR + c + ii;
      index[k] = indexValue; row[k] = -1.0;
    }
    glp_set_mat_row(lp, i, 1 + size, index, row);
  }

  // Set R_n constraints
  maxLoopSize = numC + numR;
  for (int c = 0, i = numC + 1; i<= maxLoopSize; i++, c += size) {
    index[1] = i, row[1] = 1.0;
    for (int k = 2, j = 1; j<= size; ++j, ++k){
      const int indexValue = numC + numR + c + j;
      index[k] = indexValue, row[k] = -1.0;
    }
    glp_set_mat_row(lp, i, 1 + size, index, row);
  }

  glp_term_out(0); // switch off debug output from GLPK 

  // A return value of 0 means success!
  int returnValue = glp_simplex(lp, NULL);
  int success = 0 == returnValue;

  // Fill in solution
  double perfectValue = glp_get_obj_val(lp);
  for (j = 0, i = numC + numR + 1; i<=numC + numR + numX; ++i, ++j) {
      solution[j] = glp_get_col_prim(lp, i);
  }

  double rowVals[numRows+1];
  memset(rowVals, 0.0, sizeof(rowVals));
  for (i = 1; i<=numRows; ++i){
    rowVals[i] = glp_get_col_prim(lp, i);
  }

  glp_delete_prob(lp); // release memory used for LP

  return success; /* This is not always correct, of course, and needs to be changed. */
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
