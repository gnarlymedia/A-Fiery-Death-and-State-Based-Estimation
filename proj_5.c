#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
//#include <memory.h>
//#include <errno.h>
#include "cpgplot.h"

#define pi 3.1415926
#define w 0.2625

#define ARRAY_SIZE 50272

#define initial_x 2
#define initial_y 2
#define initial_A 2
#define initial_B 2
#define initial_P_diagonal_val 0.1
#define R_diagonal_std_dev 0.06

#define initial_t_multiplier 2.0
#define final_t_multiplier 3.0

#define dist_multiplier 10000000

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
static void destroy_2d_arrayd( double *array[], int rows)
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
// Create a matrix. caller must destroy it after use.
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
static void destroy_matrix(matrix * me)
{
    destroy_2d_arrayd( me->array, me->m_nrows);
    free(me);
}

//
// Transpose a matrix. Caller must destroy returned matrix;
//
matrix * Transpose(matrix *me)
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
// Returns sum matrix on success
// Dimensions of me matrix must equal number of dimensions of other
// Result has dimensions of (me->m_nrows, me->n_cols)
// If dimensions don't match return a NULL pointer.
// User must destroy returned matrix
//
matrix * AddOrSubMatrix(matrix *me, matrix *other, char add_or_sub)
{
    matrix * sum = create_matrix(me->m_nrows, me->m_ncols);

    if (me->m_nrows != other->m_nrows)
        return NULL;
    if (me->m_ncols != other->m_ncols)
        return NULL;

    int i,j;
    for (i=0; i< me->m_nrows; i++)
    {
        for (j=0;j< me->m_ncols; j++)
        {
            if (add_or_sub == 'a') {
                sum->array[i][j] = me->array[i][j] + other->array[i][j];
            } else {
                sum->array[i][j] = me->array[i][j] - other->array[i][j];
                double test1 = me->array[i][j];
                double test2 = other->array[i][j];
                double test3 = me->array[i][j];
            }
        }
    }
    return sum;
}

//
// Multiply a matrix by a scalar
// Matrix passed in is scaled.
//
static void ScaleMatrix(matrix * me, double scale)
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
// User must destroy returned matrix
//
matrix * MultMatrix(matrix * left, matrix * right)
{
    if (left->m_ncols != right->m_nrows)
        return NULL;
    matrix *newm = create_matrix(left->m_nrows, right->m_ncols);
    int i,j,k;
    for (i=0; i< left->m_nrows; i++)
    {
        for (j=0; j< right->m_ncols; j++)
        {
            newm->array[i][j] = 0.0;
            for (k=0;k< left->m_ncols;k++)
            {
                double temp = left->array[i][k] * right->array[k][j];
                newm->array[i][j] += left->array[i][k] * right->array[k][j];
            }
        }
    }
    return newm;
}

//
// Make a copy of a matrix. Caller must destroy the matrix
//
matrix * Clone(matrix *me)
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

static void print_matrix(matrix * me)
{
    int i,j;
    for(i= 0; i< me->m_nrows; i++)
    {
        for(j=0; j < me->m_ncols; j++)
        {
            printf("%.8lf ",me->array[i][j]);
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
    matrix * cme = Clone(me);
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
    if(size == 1)
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

        // if any values are very small, return NULL
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
static void readData(double ** xval, double ** yval, int * sizeX, int * sizeY)
{
    // use on compphys server
//    FILE* fx_data = fopen("/home/shudson1/proj_5/x_data_p06.dat", "rb");
//    FILE* fy_data = fopen("/home/shudson1/proj_5/y_data_p06.dat", "rb");

    // use on virtual machine on mac
    FILE* fx_data = fopen("/home/compphys/__uni_current/x_data_p06.dat", "rb");
    FILE* fy_data = fopen("/home/compphys/__uni_current/y_data_p06.dat", "rb");

//    printf("Results of read()! %s\n", strerror(errno));
    fseek(fx_data, 0, SEEK_END);
    fseek(fy_data, 0, SEEK_END);
    int isizeX = (int) ftell(fx_data);
    int isizeY = (int) ftell(fy_data);
    (*xval) = (double *) malloc(isizeX*sizeof(double));
    (*yval) = (double *) malloc(isizeY*sizeof(double));
    rewind(fx_data);
    rewind(fy_data);
    size_t unused_1 = fread((*xval), 8, isizeX, fx_data);
    size_t unused_2 = fread((*yval), 8, isizeY, fy_data);
    (*sizeX) = isizeX;
    (*sizeY) = isizeY;
}

//
// create 1 D array of double in memory
//
static double * create_1d_array_d(int cols)
{
    double * p1dArray = malloc(cols*sizeof(double));
    return p1dArray;
}

//
// Destroy a created 1D array of double
//
static void destroy_1d_array_d(double p1dArray[])
{
    free(p1dArray);
}

static void fill_diagonal(matrix *m, double val) {
    int j, k;

    for (j = 0; j < (m->m_nrows); j = j + 1) {
        for (k = 0; k < (m->m_ncols); k = k + 1) {
            m->array[j][k] = 0.0;

            if (j == k) {
                m->array[j][k] = val;
            }
        }
    }
}

static void fill_with_rand(matrix * ma, int nrows, int ncols) {
    int j, k;

    for (j = 0; j < (ma->m_nrows); j = j + 1) {
        for (k = 0; k < (ma->m_ncols); k = k + 1) {
            ma->array[j][k] = drand48();
        }
    }
}

static void fill_with_zeros(matrix * m) {
    int j, k;

    for (j = 0; j < (m->m_nrows); j = j + 1) {
        for (k = 0; k < (m->m_ncols); k = k + 1) {
            m->array[j][k] = 0;
        }
    }
}

static void initialise_propagation_H(matrix *m) {
    fill_with_zeros(m);

    m->array[0][0] = 1;
    m->array[1][1] = 1;
}

static void update_propagation_H(matrix *m, double t) {
    m->array[0][2] = cos(w * t);
    m->array[1][3] = sin(w * t);
}

static void fill_estimate_covariance_P(matrix *m) {
    fill_diagonal(m, initial_P_diagonal_val);
}

double calculate_time(double initial_t, int counter) {
    return initial_t + (counter * (1.0/2000.0));
}

double add_fraction(double val, double add, char add_sub) {
    if (add_sub == 'a') {
        return val + (val * add);
    } else {
        return val - (val * add);
    }
}

static void disp_plot_traj(int num, double *x_vals, double *y_vals, double x_min, double x_max, double y_min,
                           double y_max, matrix * m_theta, matrix * m_P, double initial_t, double time_incr, char *heading, char *x_label, char *y_label)
{
    static float f_x_vals[ARRAY_SIZE];
    static float f_y_vals[ARRAY_SIZE];
    static float f_err_vals[ARRAY_SIZE];

    int i;

    for (i = 0; i < num; i++)
    {
        f_x_vals[i] = (float) x_vals[i];
        f_y_vals[i] = (float) y_vals[i];
        f_err_vals[i] = (float) R_diagonal_std_dev;
    }

    float fxmin, fxmax, fymin, fymax;
    fxmin = (float) x_min;
    fxmax = (float) x_max;
    fymin = (float) y_min;
    fymax = (float) y_max;

    // set up ellipse
    static float f_ell_x_vals[ARRAY_SIZE];
    static float f_ell_y_vals[ARRAY_SIZE];
    double x_0 = m_theta->array[0][0];
    double y_0 = m_theta->array[1][0];
    double A = m_theta->array[2][0];
    double B = m_theta->array[3][0];
    double w_t;

    static float f_ell_x_vals_smallest[ARRAY_SIZE];
    static float f_ell_y_vals_smallest[ARRAY_SIZE];
    double x_0_smallest = x_0 - sqrt(m_P->array[0][0]);
    double y_0_smallest = x_0 - sqrt(m_P->array[1][1]);
    double A_smallest = x_0 - sqrt(m_P->array[2][2]);
    double B_smallest = x_0 - sqrt(m_P->array[3][3]);

    static float f_ell_x_vals_biggest[ARRAY_SIZE];
    static float f_ell_y_vals_biggest[ARRAY_SIZE];
    double x_0_biggest = x_0 + sqrt(m_P->array[0][0]);
    double y_0_biggest = x_0 + sqrt(m_P->array[1][1]);
    double A_biggest = x_0 + sqrt(m_P->array[2][2]);
    double B_biggest = x_0 + sqrt(m_P->array[3][3]);

    static float f_ell_centre_x_val[1];
    f_ell_centre_x_val[0] = x_0;
    static float f_ell_centre_y_val[1];
    f_ell_centre_y_val[0] = y_0;

    double planet_centre_x = -0.1;
    double planet_centre_y = 1.3;
    double planet_radius = 0.02;
    double planet_smallest_x = planet_centre_x - (planet_radius / 2.0);
    double planet_largest_x = planet_centre_x + (planet_radius / 2.0);
    double planet_smallest_y = planet_centre_y - (planet_radius / 2.0);
    double planet_largest_y = planet_centre_y + (planet_radius / 2.0);

    int ellipse_counter = 0;
    double t;
    int impact = 0;
    int impact_sm = 0;
    int impact_lg = 0;
    double impact_x, impact_y, impact_t;
    double impact_x_sm, impact_y_sm, impact_t_sm;
    double impact_x_lg, impact_y_lg, impact_t_lg;
    double ell_x_val, ell_y_val;

    for (t = initial_t; t < (pi * 6.0); t = t + time_incr) {
        w_t = w * t;
        ell_x_val = x_0 + A * cos(w_t);
        ell_y_val = y_0 + B * sin(w_t);

        if (ell_x_val > planet_largest_x && ell_y_val > planet_largest_y) {
            // only plot if to the right and above planet
            f_ell_x_vals[ellipse_counter] = (float) ell_x_val;
            f_ell_y_vals[ellipse_counter] = (float) ell_y_val;
            ellipse_counter++;
        }

        if (ell_x_val >= planet_smallest_x && ell_x_val <= planet_largest_x && ell_y_val >= planet_smallest_y && ell_y_val <= planet_largest_y) {
            // impact
            impact = 1;
            impact_x = ell_x_val;
            impact_y = ell_y_val;
            impact_t = t;
        }
    }

    int ellipse_counter_sm = 0;
    for (t = initial_t; t < (pi * 6.0); t = t + time_incr) {
        w_t = w * t;
        ell_x_val = x_0_smallest + A_smallest * cos(w_t);
        ell_y_val = y_0_smallest + B_smallest * sin(w_t);

        if (ell_x_val > planet_largest_x && ell_y_val > planet_largest_y) {
            // only plot if to the right and above planet
            f_ell_x_vals_smallest[ellipse_counter_sm] = (float) ell_x_val;
            f_ell_y_vals_smallest[ellipse_counter_sm] = (float) ell_y_val;
            ellipse_counter_sm++;
        }

        if (ell_x_val >= planet_smallest_x && ell_x_val <= planet_largest_x && ell_y_val >= planet_smallest_y && ell_y_val <= planet_largest_y) {
            // impact
            impact_sm = 1;
            impact_x_sm = ell_x_val;
            impact_y_sm = ell_y_val;
            impact_t_sm = t;
        }
    }

    int ellipse_counter_lg = 0;
    for (t = initial_t; t < (pi * 6.0); t = t + time_incr) {
        w_t = w * t;
        ell_x_val = x_0_biggest + A_biggest * cos(w_t);
        ell_y_val = y_0_biggest + B_biggest * sin(w_t);

        if (ell_x_val > planet_largest_x && ell_y_val > planet_largest_y) {
            // only plot if to the right and above planet
            f_ell_x_vals_biggest[ellipse_counter_lg] = (float) ell_x_val;
            f_ell_y_vals_biggest[ellipse_counter_lg] = (float) ell_y_val;
            ellipse_counter_lg++;
        }

        if (ell_x_val >= planet_smallest_x && ell_x_val <= planet_largest_x && ell_y_val >= planet_smallest_y && ell_y_val <= planet_largest_y) {
            // impact
            impact_lg = 1;
            impact_x_lg = ell_x_val;
            impact_y_lg = ell_y_val;
            impact_t_lg = t;
        }
    }

    cpgbbuf();

    cpgsci(14);
    cpgenv(fxmin, fxmax, fymin, fymax, 0, 1);

    cpgsci(14);
    cpglab(x_label, y_label, heading);

    // target
    cpgsci(1);
    cpgcirc((float) planet_centre_x, planet_centre_y, planet_radius);

    // data points scatterplot
    cpgsci(1);
    cpgpt(num, f_x_vals, f_y_vals, -1);

    // error bars on data points
    cpgsci(15);
    cpgerrb(5, num, f_x_vals, f_y_vals, f_err_vals, 1.0);
    cpgerrb(6, num, f_x_vals, f_y_vals, f_err_vals, 1.0);

    // ellipse normal
    cpgsci(3);
    cpgline(ellipse_counter, f_ell_x_vals, f_ell_y_vals);

    // ellipse centre
    cpgsci(3);
    cpgpt(1, f_ell_centre_x_val, f_ell_centre_y_val, 2);

    // ellipse smallest
    cpgsci(4);
    cpgline(ellipse_counter_sm, f_ell_x_vals_smallest, f_ell_y_vals_smallest);

    // ellipse largest
    cpgsci(5);
    cpgline(ellipse_counter_lg, f_ell_x_vals_biggest, f_ell_y_vals_biggest);

    char centre_pt_label[256];
    snprintf(centre_pt_label, sizeof(centre_pt_label), "Ellipse centre at x = %.3e km, y = %.3e km", x_0, y_0);
    cpgsci(6);
    cpgtext((float) add_fraction(x_0, 0.9, 's'), (float) add_fraction(y_0, 0.06, 'a'), centre_pt_label);

    if (impact == 1) {
        // print label for impact point
        char pt_label[256];
        snprintf(pt_label, sizeof(pt_label), "Impact! at t = %.3lf yrs, x = %.2e km, y = %.2e km", impact_t, impact_x * dist_multiplier, impact_y * dist_multiplier);
        cpgsci(6);
        cpgtext((float) 0.0, (float) add_fraction(impact_y, 0.1, 'a'), pt_label);
    }

    if (impact_sm == 1) {
        // print label for impact point
        char pt_label[256];
        snprintf(pt_label, sizeof(pt_label), "Impact! at t = %.3lf yrs, x = %.2e km, y = %.2e km", impact_t_sm, impact_x_sm * dist_multiplier, impact_y_sm * dist_multiplier);
        cpgsci(7);
        cpgtext((float) 0.0, (float) add_fraction(impact_y_sm, 0.1, 's'), pt_label);
    }

    if (impact_lg == 1) {
        // print label for impact point
        char pt_label[256];
        snprintf(pt_label, sizeof(pt_label), "Impact! at t = %.3lf yrs, x = %.2e km, y = %.2e km", impact_t_lg, impact_x_lg * dist_multiplier, impact_y_lg * dist_multiplier);
        cpgsci(8);
        cpgtext((float) 0.0, (float) add_fraction(impact_y_lg, 0.2, 'a'), pt_label);
    }

    cpgebuf();
}

static void disp_plot_tr_P(int num, float *x_vals, float *y_vals, double x_min, double x_max, double y_min,
                           double y_max, char *heading, char *x_label, char *y_label)
{
    float fxmin = (float) x_min;
    float fxmax = (float) x_max;
    float fymin = (float) y_min;
    float fymax = (float) y_max;

    cpgbbuf();

    cpgsci(14);
    cpgenv(fxmin, fxmax, fymin, fymax, 0, 1);

    cpgsci(14);
    cpglab(x_label, y_label, heading);

    // ellipse
    cpgsci(3);
    cpgline(num, x_vals, y_vals);

    cpgebuf();
}

static double trace(matrix * m) {
    int j, k;
    double ret_val = 0.0;

    for (j = 0; j < (m->m_nrows); j = j + 1) {
        for (k = 0; k < (m->m_ncols); k = k + 1) {
            if (j == k) {
                ret_val += m->array[j][k];
            }
        }
    }

    return ret_val;
}

static void findAndSetMax(double val, double * max_val) {
    if (val > * max_val) {
        * max_val = val;
    }
}

int main(void)
{
    srand48((int) time(NULL));

    double initial_t = pi * initial_t_multiplier;
    double final_t = pi * final_t_multiplier;

    double R_diagonal_variance_val  = pow(R_diagonal_std_dev, 2.0);

    double * data_meas_arr_x = NULL;
    double * data_meas_arr_y = NULL;
    int data_size_x = 0;
    int data_size_y = 0;

    readData(&data_meas_arr_x, &data_meas_arr_y, &data_size_x, &data_size_y);

//    printf("First two measurements - x: %lf, y: %lf\n\n", data_meas_arr_x[0], data_meas_arr_y[0]);

    // theta
    matrix * state_vector_theta_k = create_matrix(4, 1);
    matrix * state_vector_theta_k_minus_1 = create_matrix(4, 1);
    state_vector_theta_k_minus_1->array[0][0] = initial_x;
    state_vector_theta_k_minus_1->array[1][0] = initial_y;
    state_vector_theta_k_minus_1->array[2][0] = initial_A;
    state_vector_theta_k_minus_1->array[3][0] = initial_B;

    printf("Initial value of state_vector_theta_k_minus_1\n");
    print_matrix(state_vector_theta_k_minus_1);
    printf("\n\n");

    matrix * measurement_m = create_matrix(2, 1);
    fill_with_zeros(measurement_m);

    matrix * measurement_covariance_R = create_matrix(2, 2);
    fill_diagonal(measurement_covariance_R, R_diagonal_variance_val);

    matrix * propagation_H = create_matrix(2, 4);

    matrix * residuals_e = create_matrix(2, 1);

    matrix * estimate_covariance_P_k = create_matrix(4, 4);

    matrix * estimate_covariance_P_k_minus_1 = create_matrix(4, 4);
    fill_estimate_covariance_P(estimate_covariance_P_k_minus_1);

    matrix * covariance_in_residuals_S = create_matrix(2, 2);

    matrix * kalman_gain_K = create_matrix(4, 2);

    matrix * identity_I = create_matrix(4, 4);
    fill_diagonal(identity_I, 1.0);

    double t_k, A, B;
    double time_incr = 1.0 / 2000.0;

    static double plot_x_vals[ARRAY_SIZE];
    static double plot_y_vals[ARRAY_SIZE];

    static float plot_tr_P_vals[ARRAY_SIZE];
    static float plot_tr_P_t_vals[ARRAY_SIZE];

    double biggest_meas_val_x = 0.0;
    double biggest_meas_val_y = 0.0;

    double biggest_meas_tr_P_val = 0.0;
    double biggest_meas_tr_P_time_val = 0.0;

    double state_vector_theta_x, state_vector_theta_y;
    double measurement_x, measurement_y;
    double trace_P_val;

    initialise_propagation_H(propagation_H);

    int counter = 0;
    for (t_k = initial_t; t_k <= final_t; t_k = t_k + time_incr) {
        measurement_x = data_meas_arr_x[counter];
        measurement_y = data_meas_arr_y[counter];

        measurement_m->array[0][0] = measurement_x;
        measurement_m->array[1][0] = measurement_y;

        update_propagation_H(propagation_H, t_k);

        residuals_e = AddOrSubMatrix(measurement_m, MultMatrix(propagation_H, state_vector_theta_k_minus_1), 's');

        covariance_in_residuals_S = AddOrSubMatrix(MultMatrix(MultMatrix(propagation_H, estimate_covariance_P_k_minus_1),
                                                          Transpose(propagation_H)), measurement_covariance_R, 'a');

        kalman_gain_K = MultMatrix(MultMatrix(estimate_covariance_P_k_minus_1, Transpose(propagation_H)), Invert(covariance_in_residuals_S));

//        if (t_k == initial_t) {
//            printf("residuals_e\n");
//            print_matrix(residuals_e);
//            printf("\n");
//        }

        state_vector_theta_k = AddOrSubMatrix(state_vector_theta_k_minus_1, MultMatrix(kalman_gain_K, residuals_e), 'a');

        estimate_covariance_P_k = MultMatrix(AddOrSubMatrix(identity_I, MultMatrix(kalman_gain_K, propagation_H), 's'), estimate_covariance_P_k_minus_1);

        state_vector_theta_x = state_vector_theta_k->array[0][0];
        state_vector_theta_y = state_vector_theta_k->array[1][0];

        measurement_x = data_meas_arr_x[counter];
        measurement_y = data_meas_arr_y[counter];

        plot_x_vals[counter] = measurement_x;
        plot_y_vals[counter] = measurement_y;

        findAndSetMax(measurement_x, &biggest_meas_val_x);
        findAndSetMax(measurement_y, &biggest_meas_val_y);

        trace_P_val = trace(estimate_covariance_P_k);

        plot_tr_P_vals[counter] = trace_P_val;
        plot_tr_P_t_vals[counter] = t_k;

        findAndSetMax(trace_P_val, &biggest_meas_tr_P_val);
        findAndSetMax(t_k, &biggest_meas_tr_P_time_val);

        state_vector_theta_k_minus_1 = state_vector_theta_k;
        estimate_covariance_P_k_minus_1 = estimate_covariance_P_k;

        counter++;
    }

    printf("Initial value of estimate_covariance_P_k diagonal: %lf\n", initial_P_diagonal_val);
    printf("Initial value of measurement_covariance_R diagonal: %lf\n\n", R_diagonal_variance_val);

    printf("Final value of estimate_covariance_P_k\n");
    print_matrix(estimate_covariance_P_k);
    printf("\n");

//    printf("Final value of t: %lf\n", t);
//    printf("Value of data_size_x: %d\n", data_size_x);

    printf("Final trace value of P_k_minus_1: %.8lf\n\n", trace_P_val);

//    printf("Final value of state_vector_theta_k\n");
//    print_matrix(state_vector_theta_k);
//    printf("\n");

    char plot_file_name[256];
    snprintf(plot_file_name, sizeof(plot_file_name), "proj5plot-trP-%.0lf-%.0lf.ps/CPS", initial_t_multiplier, final_t_multiplier);

    //    if (cpgbeg(0, "?", 1, 1) != 1) {
    if (cpgbeg(0, "/XWINDOW", 1, 1) != 1) {
//    if (cpgbeg(0, "proj5plot.ps/CPS", 1, 1) != 1) {
        exit(EXIT_FAILURE);
    }
    cpgask(1);

    disp_plot_traj(counter, plot_x_vals, plot_y_vals, -0.5, add_fraction(biggest_meas_val_x, 0.1, 'a'), 1.0, add_fraction(biggest_meas_val_y, 0.1, 'a'), state_vector_theta_k, estimate_covariance_P_k, initial_t, time_incr, "Plot of the trajectory of the asteroid", "Distance (units of 10,000,000 km)", "Distance (units of 10,000,000 km)");

    char plot_label[256];
    snprintf(plot_label, sizeof(plot_label), "Trace of P versus t, for t: pi * %.1lf - pi * %.1lf (%.3lf - %.3lf to 3 dec. plc.)", initial_t_multiplier, final_t_multiplier, initial_t, final_t);

    disp_plot_tr_P(counter, plot_tr_P_t_vals, plot_tr_P_vals, initial_t, add_fraction(biggest_meas_tr_P_time_val, 0.1, 'a'), 0.0, add_fraction(biggest_meas_tr_P_val, 0.1, 'a'), plot_label, "Time (years)", "Trace of propagation matrix P (x10E14 km squared)");

    cpgend();
}