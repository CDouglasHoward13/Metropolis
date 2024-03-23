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
// This code finds the minimum variance portfolio with no short positions
// permitted, as explained in Section 18.
////////////////////////////////////////////////////////////////////////////////

// These functions are found below.
void     GetData ();
double   Energy ();
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
double epsilon = 0.001, *x, **V;

#include "MetropolisFunctions.h"

///////////////////////////////////////////////////////////////////////////////
// Main program.
///////////////////////////////////////////////////////////////////////////////
int main () {

   // Calculate the covariance matrix.
   GetData ();

   // Best portfolio via Metropolis.
   Metropolis ();

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
   printf ("I'm looking for the minimum variance no-shorts portfolio.\n");
   MTUniform ();
   printf ("\nI'll be done when I find a stable state. ");

   // Start with everything in the portfolio at $2.  This is arbitrary. The
   //    starting portfolio should not affect the outcome.
   for (i = 1; i <= 50; i++) {
      x[i] = 2.0;
   }

   // Compute the initial variance.
   E = Energy ();

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
      E_new = Energy ();

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
   for (i = 1; i <= 50; i++) if (x[i] > epsilon/2.0) {
         printf ("%8.2f  ", x[i]);
         printf (ticker[i]);
    }
   printf ("\n");

   printf ("The smallest variance found via Metropolis is %.5f\n", Energy ());

   Pause ();

}

////////////////////////////////////////////////////////////////////////////////
// Compute the energy of the x portfolio.
////////////////////////////////////////////////////////////////////////////////
double Energy () {

   int i, j;
   double variance=0;

   for (i = 1; i <= 50; i++)  {
      if (x[i] < -epsilon/2.0) return 1000.0;
      for (j = 1; j <= 50; j++) {
         variance += x[i] * V[i][j] * x[j];
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
   E0 = Energy ();

   // Does x have a neighbor with lower energy?
   for (i = 1; i <= 50; i++) {
      for (j = 1; j <= 50; j++) if (i != j) {

         x[i] += epsilon;
         x[j] -= epsilon;
         E = Energy ();
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

   // Allocate array space.
   V     = Array (50, 50);
   x     = (double *) calloc (51, sizeof (double));

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


