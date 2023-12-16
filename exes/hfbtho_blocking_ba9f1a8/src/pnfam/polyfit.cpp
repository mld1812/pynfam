/* ****************************************************************************
   polyfit.cpp
   
   Subroutine for fitting a polynomial Sum_{n=1,N} (x-1)^n to a set of points;
   this version is designed to be easily wrappable to Fortran.  The purpose of
   the routine is to help finding a good polynomial approximation to the phase
   space integrals.
   
   M.T. Mustonen and T. Shafer, UNC Chapel Hill, 2013-15
**************************************************************************** */

#include <gsl/gsl_vector.h>
#include <gsl/gsl_multifit.h>

extern "C" {

// Input parameters:
// nr_fit = the order of the polynomial (= the number of terms, because the
//          constant term is omitted to force f(x=1) = 0).
// nr_data_points = number of points to fit to
// xs = x values of the data points
// ys = y values of the data points

// Output parameters:
// cs = fitted constants
int fit_poly(const int nr_fit, const int nr_data_points, const double xs[],
   const double ys[], double cs[]) {
   
   // Move data into GSL vector types (required by the fitting routine)
   gsl_vector *x = gsl_vector_alloc(nr_data_points);
   gsl_vector *y = gsl_vector_alloc(nr_data_points);
   for(int i = 0; i < nr_data_points; ++i) {
      gsl_vector_set(x, i, xs[i]); gsl_vector_set(y, i, ys[i]);
   }
   
   // Utilize GSL to fit a polynomial
   gsl_matrix* X = gsl_matrix_alloc(nr_data_points, nr_fit);
   gsl_vector* c = gsl_vector_alloc(nr_fit);
   gsl_matrix* cov = gsl_matrix_alloc(nr_fit, nr_fit);
   double chisq;
   // This is where we define the curve to fit:
   for (int i = 0; i < nr_data_points; ++i) {
      double aux = 1;
      double xval = gsl_vector_get(x, i);
      for (int j = 0; j < nr_fit; ++j) {
         aux = aux*(xval - 1);
         gsl_matrix_set(X, i, j, aux);
      }
   }
   gsl_multifit_linear_workspace *work = gsl_multifit_linear_alloc(nr_data_points, nr_fit);
   gsl_multifit_linear(X, y, c, cov, &chisq, work);
   gsl_multifit_linear_free(work);
   
   // Store the result on the return array
   for (int i = 0; i < nr_fit; ++i) {
      cs[i] = gsl_vector_get(c, i);
   }
   
   gsl_vector_free(x);
   gsl_vector_free(y);
   gsl_matrix_free(X);
   gsl_matrix_free(cov);
   gsl_vector_free(c);
   
   return 0;
}

}
