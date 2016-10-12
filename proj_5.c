#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include "cpgplot.h"

#ifndef EXIT_FAILURE
#define EXIT_FAILURE 1
#endif
#ifndef EXIT_SUCCESS
#define EXIT_SUCCESS 0
#endif


typedef struct _matrix matrix;

struct _matrix
{
    int m_nrows;
    int m_ncols;
    double ** array;
};
//
// Create an array of dimensions (rows,cols) of double in memory
//
double ** create_2d_arrayd(int rows, int cols)
{
    double ** p2dArray;
    p2dArray = (double **) malloc(rows*sizeof(double *));
    int i = 0;
    for(i=0; i< rows; i++)
    {
        p2dArray[i] = (double *) malloc(cols*sizeof(double));
    }
    return p2dArray;
}
//
// Destroy a 2D array of dimensions (rows,[]) of double in memory
//
void destroy_2d_arrayd( double *array[], int rows)
{
    int i=0;
    for(i=0;i<rows;i++)
    {
        free(array[i]);
    }
    free(array);
}

//
// return the value of the matrix me at (row,col)
//
double val(matrix * me, int row, int col)
{
    return me->array[row][col];
}

//
// Create a matrix. caller must detroy it after use.
//
matrix * create_matrix(int nrows, int ncols)
{
    matrix * my_matrix = malloc(sizeof(matrix));
    my_matrix->m_nrows = nrows;
    my_matrix->m_ncols = ncols;
    my_matrix->array = create_2d_arrayd(nrows,ncols);
    return my_matrix;
}

//
// destory a matrix
//
void destroy_matrix(matrix * me)
{
    destroy_2d_arrayd( me->array, me->m_nrows);
    free(me);
}

//
// Transpose a matrix. Caller must destroy returned matrix;
//
matrix * transpose(matrix * me)
{
    matrix * trans = create_matrix(me->m_ncols,me->m_nrows);
    int i, j;
    for(i=0; i< me->m_nrows;i++)
    {
        for(j=0; j < me->m_ncols; j++)
        {
            trans->array[j][i] = me->array[i][j];
        }
    }
    return trans;
}

//
// Add a matrix to a matrix.
// Return 0 if dimensions don't match.
// Returns 1 on success
//
int AddMatrix(matrix * me, matrix * other)
{
    if(me->m_nrows != other->m_nrows)
        return 0;
    if(me->m_ncols != other->m_ncols)
        return 0;
    int i,j;
    for(i=0; i< me->m_nrows; i++)
    {
        for(j=0;j< me->m_ncols; j++)
        {
            me->array[i][j] += other->array[i][j];
        }
    }
    return 1;
}

//
// Multiply a matrix by a scalar
// Matrix passed in is scaled.
//
void ScaleMatrix(matrix * me, double scale)
{
    int i,j;
    for(i=0; i < me->m_nrows;i++)
    {
        for(j=0; j< me->m_ncols; j++)
        {
            me->array[i][j] = me->array[i][j]*scale;
        }
    }
}
//
// Multiply two matricies.
// Number of columns of left matrix must equal number of rows of right
// Result has dimensions of (left->m_nrows,right->n_cols)
// If rows and cols don;t match return a NULL pointer.
// User must destory returned matrix
//
matrix * MultMatrix(matrix * left, matrix * right)
{
    if(left->m_nrows != right->m_ncols)
        return NULL;
    matrix *newm = create_matrix(left->m_nrows,right->m_ncols);
    int i,j,k;
    for(i=0; i< left->m_nrows; i++)
    {
        for(j=0; j< right->m_ncols; j++)
        {
            newm->array[i][j] = 0.0;
            for(k=0;k< left->m_ncols;k++)
            {
                newm->array[i][j] += left->array[i][k]*right->array[k][j];
            }
        }
    }
    return newm;
}

//
// Make a copy of a matrix. Caller must destroy the matrix
//
matrix * clone( matrix * me)
{
    matrix * cme = create_matrix(me->m_nrows,me->m_ncols);
    int i,j;
    for(i=0; i< me->m_nrows; i++)
    {
        for(j= 0; j< me->m_ncols; j++)
        {
            cme->array[i][j] = me->array[i][j];
        }
    }
    return cme;
}

void print_matrix(matrix * me)
{
    int i,j;
    for(i= 0; i< me->m_nrows; i++)
    {
        for(j=0; j < me->m_ncols; j++)
        {
            printf("%8.4lf ",me->array[i][j]);
        }
        printf("\n");
    }
}

//
// Invert a matrix
// Caller should destroy inverse matrix
// If the input matrix is singular, the function returns a NULL pointer
//
matrix * Invert(matrix * me)
{
    if(me->m_nrows != me->m_ncols)
        return NULL;
    matrix * inverse =  create_matrix(me->m_nrows,me->m_ncols);
    int i,j,k;
    int size = me->m_nrows;
    matrix * cme = clone(me);
    double p = 0.0;
    double q = 0.0;

    for(i=0; i< size; i++)
    {
        for(j= 0; j< size; j++)
        {
            inverse->array[i][j] = 0.0;
            if(i==j)
            {
                inverse->array[i][j] = 1.0;
            }
        }
    }
    if( size == 1)
    {
        inverse->array[0][0] = 1.0/me->array[0][0];
        return inverse;
    }
    //
    // Put into lower triangular form
    //
    for(i=0;i<size;i++)
    {
        p = cme->array[i][i];
        if(fabs(p) < 1.0e-100)
        {
            destroy_matrix(inverse);
            destroy_matrix(cme);
            return NULL;
        }
        for(j=i;j< size;j++)
        {
            if(j == i)
            {
                for(k= 0; k< size; k++)
                {
                    cme->array[i][k] = cme->array[i][k]/p;
                    inverse->array[i][k] = inverse->array[i][k]/p;
                }
            }
            if( j != i)
            {
                q = cme->array[j][i];
                for(k =0; k < size; k++)
                {
                    cme->array[j][k] = cme->array[j][k] - q*cme->array[i][k];
                    inverse->array[j][k] = inverse->array[j][k] - q*inverse->array[i][k];
                }
            }
        }
    }
    //
    // Upper diagonal
    //
    for( i = size-1; i > 0; i--)
    {
        for(j = i-1; j>= 0; j--)
        {
            q = cme->array[j][i];
            for(k = 0; k < size; k++)
            {
                cme->array[j][k] = cme->array[j][k] - q*cme->array[i][k];
                inverse->array[j][k] = inverse->array[j][k] - q*inverse->array[i][k];
            }
        }
    }
    destroy_matrix(cme);
    return inverse;
}

//
// Reads data from the input files
// Creates the arrays xval[] and yval[] and loads these with data
// returns the number of elements in xval and yval in the pointers
// sizeX and sizeY
//
void readData( double xval[], double yval[], int * sizeX, int * sizeY)
{
    FILE* fx_data = fopen("x_data_p06.dat", "rb");
    FILE* fy_data = fopen("y_data_p06.dat", "rb");
    fseek(fx_data, 0, SEEK_END);
    fseek(fy_data, 0, SEEK_END);
    int isizeX = (int) ftell(fx_data);
    int isizeY = (int) ftell(fy_data);
    xval = malloc(isizeX*sizeof(double));
    yval = malloc(isizeY*sizeof(double));
    rewind(fx_data);
    rewind(fy_data);
    fread(xval, 8, isizeX, fx_data);
    fread(yval, 8, isizeY, fy_data);
    (*sizeX) = isizeX;
    (*sizeY) = isizeY;
}

int main(void)
{
    int rows, cols;
    matrix * testm = NULL;
    matrix * testi = NULL;
    matrix * multmi = NULL;
    int i,j;
    char inp[20];
    inp[0] = 'y';
    srand48(time(NULL));
    while(inp[0] == 'y')
    {
        printf("input number rows, cols :");
        scanf("%d,%d",&rows,&cols);
        //    printf("Number rows %d no cols %d \n",rows,cols);
        testm = create_matrix(rows,cols);
        printf("Created array \n");
        for(i=0;i<rows;i++)
        {
            for(j=0;j<cols;j++)
            {
                testm->array[i][j] = 1.0 - 2.0*drand48();
            }
        }
        printf("Initial is \n");
        print_matrix(testm);
        if(rows == cols)
        {
            testi = Invert(testm);
            if(testi != NULL)
            {
                printf("Inverted matrix is \n");
                print_matrix(testi);
                printf("Matrix times invert is \n");
                multmi =  MultMatrix(testm, testi);
                print_matrix(multmi);
                destroy_matrix(testi);
                destroy_matrix(multmi);
            }
            else
            {
                printf("Matrix is singular \n");
            }
        }
        else
        {
            testi = transpose(testm);
            printf("transposed matrix is \n");
            print_matrix(testi);
            multmi =  MultMatrix(testi,testm);
            printf("transposed matrix times original is \n");
            print_matrix(multmi);
            destroy_matrix(testi);
            destroy_matrix(multmi);
        }
        destroy_matrix(testm);
        printf("Another try? (y or n):");
        scanf("%s",inp);
    }
}