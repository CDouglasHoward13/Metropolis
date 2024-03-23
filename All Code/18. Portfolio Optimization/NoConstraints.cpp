/*
MIT License

Copyright (c) 2024 CDouglasHoward13

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

////////////////////////////////////////////////////////////////////////////////
// This code finds the minimum variance portfolio with no constraints, as
// explained in Section 18.
////////////////////////////////////////////////////////////////////////////////

// These functions are found below.
void     GetData ();
double   Energy (double *);
void     Optimal ();
int      Stable ();
void     Metropolis ();
void     Report ();

// These functions, also below, involve matrix manipulation.
double **Array (int, int);
double **Multiply (double **, double **);
double **Invert (double **);
double **Transpose (double **);
double **Copy (double **);
double **Identity (int);
int      Rows (double **);
int      Columns (double **);
void     Free (double **);

// Global variables:
char **ticker;
double epsilon = 0.001, *x, *xstar, **V;

#include "MetropolisFunctions.h"

///////////////////////////////////////////////////////////////////////////////
// Main program.
///////////////////////////////////////////////////////////////////////////////
int main () {

   // Calculate the covariance matrix.
   GetData ();

   // Best portfolio via Metropolis.
   Metropolis ();

   // Calculate the true optimal using quadratic optimization.
   Optimal ();

   // Report the results.
   Report ();

}

///////////////////////////////////////////////////////////////////////////////
// Compute the minimum variance portfolio via Metropolis. Here we use the
//   zero temperature dynamics.
///////////////////////////////////////////////////////////////////////////////
void Metropolis () {

   int i, j;
   double t, t1, E, E_new;

   // Seed the RNG and get the temperature.
   printf ("I'm looking for the minimum variance unconstrained portfolio.\n");
   MTUniform ();
   printf ("\nI'll be done when I find a stable state. ");

   // Start with everything in the portfolio at $2.  This is arbitrary. The
   //    starting portfolio should not affect the outcome.
   for (i = 1; i <= 50; i++) {
      x[i] = 2.0;
   }

   // Compute the initial variance.
   E = Energy (x);

   // Start the timer.
   t1 = Time ();

   // Run the Markov chain until a stable state is found.
   while (1) {

      // Every five seconds indicate that it's still thinking. Break if the state is stable.
      t = Time ();
      if (t > t1 + 5.0) {
         if (Stable ()) break;
         printf (". ");
         t1 = t;
      }

      // Select a stock at random to decrease.
      i = RandomInteger (1, 50);

      // Select a different stock at random to increase.
      j = i;
      while (j == i) {
         j = RandomInteger (1, 50);
      }

      //  Alter the portfolio.
      x[i] -= epsilon;
      x[j] += epsilon;

      // Compute the new energy.
      E_new = Energy (x);

      // If not worse, accept the change.
      // Use zero temperature dynamics here.
      // The energy decreases monotonically in this application.
      if (E_new <= E) {
         E = E_new;
      }

      // Otherwise restore previous configuration.
      else {
         x[i] += epsilon;
         x[j] -= epsilon;
      }

   }

}

///////////////////////////////////////////////////////////////////////////////
// Report the results.
///////////////////////////////////////////////////////////////////////////////
void Report () {

   int i;

   // Report the best found portfolio and the true optimal.
   printf ("\n\n");
   printf ("                  True\n");
   printf ("Metropolis       Optimal\n");
   printf ("==========    ==========\n");
   for (i = 1; i <= 50; i++) {
         printf ("%8.2f      %8.2f  ", x[i], xstar[i]);
         printf (ticker[i]);
    }
   printf ("\n");

   printf ("The smallest variance found via Metropolis is %.5f\n", Energy (x));
   printf ("True optimal portfolio variance is %.5f\n", Energy (xstar));

   Pause ();

}

///////////////////////////////////////////////////////////////////////////////
// Using quadratic optimization compute the true optimal portfolio.
///////////////////////////////////////////////////////////////////////////////
void Optimal () {

   int i;
   double c;

   double **e, **VInv, **opt;

   e    = Array (50, 1);
   for (i = 1; i <= 50; i++) {
      e[i][1] = 1;
   }

   VInv = Invert (V);

   opt = Multiply (VInv, e);

   c = Multiply (Transpose(e), opt)[1][1];

   for (i = 1; i <= 50; i++) {
      xstar[i] = 100.0 * opt[i][1] / c;
   }

   return;

}

////////////////////////////////////////////////////////////////////////////////
// Compute the return variance of the p portfolio.
////////////////////////////////////////////////////////////////////////////////
double Energy (double *p) {

   int i, j;
   double variance=0;

   for (i = 1; i <= 50; i++)  {
      for (j = 1; j <= 50; j++) {
         variance += p[i] * V[i][j] * p[j];
      }
   }

   return variance;

}

////////////////////////////////////////////////////////////////////////////////
// Check to see if the portfolio x is stable.
////////////////////////////////////////////////////////////////////////////////
int Stable () {

   int i, j;
   double E0, E;

   // Current energy.
   E0 = Energy (x);

   // Does x have a neighbor with lower energy?
   for (i = 1; i <= 50; i++) {
      for (j = 1; j <= 50; j++) if (i != j) {

         x[i] += epsilon;
         x[j] -= epsilon;
         E = Energy (x);
         x[i] -= epsilon;
         x[j] += epsilon;
         if (E < E0) return 0;

      }
   }

   return 1;

}

////////////////////////////////////////////////////////////////////////////////
// Allocate space for and get stock price return covariance data.
////////////////////////////////////////////////////////////////////////////////
void GetData () {

   int i, j;
   double V0;
   char input[100];
   FILE *fp;

   // Allocate array space. These are global variables.
   V     = Array (50, 50);
   x     = (double *) calloc (51, sizeof (double));
   xstar = (double *) calloc (51, sizeof (double));

   // Allocate space to hold ticker names; "ticker" is a global variable.
   ticker = (char **) calloc (50+1, sizeof (char *));
   for (i = 0; i <= 50; i++) {
      ticker[i] = (char *) calloc (10, sizeof (char));
   }

   // Read in the data.
   fp = fopen ("V.txt", "r");
   for (i = 1; i <= 50; i++) {

      // Read in stock i's covariance data.

      // Name of the stock ticker.
      fgets (ticker[i], 8, fp);

      // The covariances.
      for (j = 1; j <= 50; j++) {

         // Read in V[i][j].
         fgets (input, 99, fp);
         sscanf (input, "%lf", &V0);

         // Put data into the V array.
         V[i][j] = V0;

      }

   }

   fclose (fp);

   // Send the data back to the main program.
   return;

}

////////////////////////////////////////////////////////////////////////////////
// Allocate array space for an m x n array.
////////////////////////////////////////////////////////////////////////////////
double **Array (int m, int n) {

   int i;
   double **A;

   A = (double **) calloc (m+1, sizeof (double *));
   for (i = 0; i <= m; i++) {
      A[i] = (double *) calloc (n+1, sizeof (double));
   }

   // Record dimensions in the 0^th row, which is unused in matrix operations.
   A[0][0] = m;
   A[0][1] = n;

   return (A);

}






//**** The following functions are used only in the function Optimal. ****//

////////////////////////////////////////////////////////////////////////////////
// Compute the transpose of the matrix A.
////////////////////////////////////////////////////////////////////////////////
double **Transpose (double **A) {

   int i, j, n, m;
   double **At;

   m = Rows (A);
   n = Columns (A);

   At = Array (n,m);

   for (i = 1; i <= n; i++) {
      for (j = 1; j <= m; j++) {
         At[i][j] = A[j][i];
      }
   }

   return At;

}

////////////////////////////////////////////////////////////////////////////////
// Use Gaussian elimination to invert the square matrix A0.
////////////////////////////////////////////////////////////////////////////////
double **Invert (double **A0) {

   int i, j, k, n, rmax;
   double **A, **Ainv, c;

   // Make a copy so that the original matrix is not changed during G.E.
   A = Copy (A0);

   n = Rows (A);

   if (n != Columns (A)) {
      printf ("Trying to invert a non-square matrix.\n");
      Pause ();
      exit(1);
   }

   // Start with the n x n identity matrix.  This matrix will eventually hold
   //    the inverse of A.
   Ainv = Identity (n);

   // Work on column j of A.
   for (j = 1; j <= n; j++) {

      // Find the largest non-zero entry in the column on or below the diagonal.
      c = 0;
      rmax = 0;
      for (i = j; i <= n; i++) {
         if (fabs(A[i][j]) > c) {
            c = fabs (A[i][j]);
            rmax = i;
         }
      }

      // If they are all 0 the matrix is singular.
      if (c < 1e-10) {
         printf ("Trying to invert a singular matrix.\n");
         Exit();
      }

      // Swap rows j and rmax in both A and Ainv.
      i = j;
      for (k = 1; k <= n; k++) {
         c = A[i][k];
         A[i][k] = A[rmax][k];
         A[rmax][k] = c;
         c = Ainv[i][k];
         Ainv[i][k] = Ainv[rmax][k];
         Ainv[rmax][k] = c;
      }

      // Scale so the pivot is 1.
      c = A[i][i];
      for (k = 1; k <= n; k++) {
         A[i][k] /= c;
         Ainv[i][k] /= c;
      }

      // Make rest of column j equal to 0. Apply same row operations to Ainv.
      for (i = 1; i <= n; i++) if (i != j) {
         c = A[i][j];
         for (k = 1; k <= n; k++) {
            A[i][k] += -c * A[j][k];
            Ainv[i][k] += -c * Ainv[j][k];
         }
      }

   }

   // We're done with A.
   Free (A);

   return Ainv;

}

////////////////////////////////////////////////////////////////////////////////
// Compute the matrix product AB.
////////////////////////////////////////////////////////////////////////////////
double **Multiply (double **A, double **B) {

   int i, j, k, nA, mA, nB, mB;
   double **AB;

   mA = Rows (A);
   nA = Columns (A);

   mB = Rows (B);
   nB = Columns (B);

   if (nA != mB) {
      printf ("Dimensions don't match in matrix multiplication.\n");
      Pause ();
      exit(1);
   }

   // Allocate space for the product.
   AB = Array (mA, nB);

   // Compute (AB)_ij.
   for (i = 1; i <= mA; i++) {
      for (j = 1; j <= nB; j++) {
         for (k = 1; k <= nA; k++) {
            AB[i][j] += A[i][k] * B[k][j];
         }
      }
   }

   return AB;

}

////////////////////////////////////////////////////////////////////////////////
// Report the number of rows in the matrix A.  This is held in A[0][0].  Note
//   that row 0 of A is not used in any matrix computations.
////////////////////////////////////////////////////////////////////////////////
int Rows (double **A) {

   return (int) A[0][0];

}

////////////////////////////////////////////////////////////////////////////////
// Report the number of columns in the matrix A. Held in A[0][1].
////////////////////////////////////////////////////////////////////////////////
int Columns (double **A) {

   return (int) A[0][1];

}

////////////////////////////////////////////////////////////////////////////////
// Copy the matrix A into another matrix (C).
////////////////////////////////////////////////////////////////////////////////
double **Copy (double **A) {

   double **C;
   int i, j, m, n;

   // Allocate space for C.
   m = Rows (A);
   n = Columns (A);
   C = Array (m, n);

   // Copy A_ij into C_ij.
   for (i = 1; i <= m; i++) {
      for (j = 1; j <= n; j++) {
         C[i][j] = A[i][j];
      }
   }

   return C;

}

////////////////////////////////////////////////////////////////////////////////
// Generate the n x n identity matrix.
////////////////////////////////////////////////////////////////////////////////
double **Identity (int n) {

   double **I;
   int i, j;

   I = Array (n, n);

   for (i = 1; i <= n; i++) {
      for (j = 1; j <= n; j++) {
         I[i][j] = (i == j);
      }
   }

   return I;

}

////////////////////////////////////////////////////////////////////////////////
// Release allocated space for the matrix A.
////////////////////////////////////////////////////////////////////////////////
void Free (double **A) {

   int i, m;

   m = Rows (A);
   for (i = 0; i <= m; i++) {
      free (A[i]);
   }
   free (A);

   return;

}








