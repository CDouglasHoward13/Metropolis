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
// This code finds the minimum variance "simple" portfolio, as explained
// in Section 18.
////////////////////////////////////////////////////////////////////////////////

// These functions are found below.
void     GetData ();
double   Energy ();
void     Optimal ();
int      Stable ();
void     Metropolis ();
void     Report ();
double **Array (int, int);

// Global variables:
char **ticker;
double **V;
int *x, *x_min;

#include "MetropolisFunctions.h"

////////////////////////////////////////////////////////////////////////////////
// Main program.
////////////////////////////////////////////////////////////////////////////////
int main () {

   // Calculate the covariance matrix.
   GetData ();

   // Best portfolio via Metropolis.
   Metropolis ();

   // Report the results.
   Report ();

}

///////////////////////////////////////////////////////////////////////////////
// Compute the minimum variance portfolio via Metropolis. Here we cannot use
//   the zero temperature dynamics.
///////////////////////////////////////////////////////////////////////////////
void Metropolis () {

   int i, j;
   double T, t, t1, E, E_new, E_min, U, p, accept;


   // Start with everything in the portfolio at $2.  This is arbitrary. The
   //    starting portfolio should not affect the outcome.
   for (i = 1; i <= 50; i++) {
      x[i] = 1;
   }

   // Compute the initial variance.
   E_min = E = Energy ();

   // Seed the RNG and get the temperature.
   printf ("I'm looking for the minimum variance simple portfolio.\n");
   MTUniform ();
   T = GetDouble ("\nWhat is the temperature (0.1 is good)?... ");

   printf ("\nI'll be done in 60 seconds. ");
   t = t1 = Time ();

   // Run the Markov chain for 60 seconds.
   while (t < 60.0) {

      // Every five seconds indicate that it's still thinking.
      t = Time ();
      if (t > t1 + 5.0) {
         printf (". ");
         t1 = t;
      }

      // Select a stock at random to "flip".
      i = RandomInteger (1, 50);

      // "Flip" that stock.
      x[i] = 1 - x[i];

      // Compute the new energy.
      E_new = Energy ();

      accept = 0;

      // If not worse, accept the change.
      if (E_new <= E) {
         accept = 1;
         // See if energy is a new minimum.  If so record data.
         if (E_new < E_min) {
            E_min = E_new;
            for (j = 1; j <= 50; j++) {
               x_min[j] = x[j];
            }
         }
      }

      // If worse and T>0, accept with some positive probability.
      else if (T > 0.0) {
         U = MTUniform ();
         p = exp (-(E_new-E)/T);
         accept = (U <= p);
      }

      // If accepted, keep new configuration and update energy.
      if (accept) {
         E = E_new;
      }

      // Otherwise restore previous configuration.
      else {
         x[i] = 1 - x[i];
      }

   }

   // Copy x_min into x.
   for (i = 1; i <= 50; i++) {
      x[i] = x_min[i];
   }

}

///////////////////////////////////////////////////////////////////////////////
// Report the results.
///////////////////////////////////////////////////////////////////////////////
void Report () {

   int i, n;

   // How many stocks are in the portfolio?
   n = 0;
   for (i = 1; i <= 50; i++) if (x[i]) {
      n++;
   }

   // Report the best found portfolio and the true optimal.
   printf ("\n\n");
   for (i = 1; i <= 50; i++) if (x[i]) {
         printf ("%8.2f  ", (100.0/n));
         printf (ticker[i]);
    }
   printf ("\n");

   printf ("The smallest variance found via Metropolis is %.5f\n", Energy ());

   // See if this portfolio is stable.
   printf ("\n");
   if (Stable ()) {
      printf ("This portfolio is a stable state.\n");
   } else {
      printf ("This portfolio is not a stable state.\n");
   }   


   Pause ();

}

////////////////////////////////////////////////////////////////////////////////
// Compute the energy of the x portfolio.
////////////////////////////////////////////////////////////////////////////////
double Energy () {

   int i, j, n;
   double variance;

   n = variance = 0;
   for (i = 1; i <= 50; i++)  if (x[i]) {
      n ++;
      for (j = 1; j <= 50; j++) if (x[j]) {
         variance += V[i][j]; // x[i] = x[j] = 1 here.
      }
   }

   if (!n) return 1000.0;
   else    return pow (100.0/n, 2.0)*variance;

}

////////////////////////////////////////////////////////////////////////////////
// Check to see if the portfolio x is stable.
////////////////////////////////////////////////////////////////////////////////
int Stable () {

   int i;
   double E0, E;

   // Current energy.
   E0 = Energy ();

   // Does x have a neighbor with lower energy?
   for (i = 1; i <= 50; i++) {

      x[i] = 1 - x[i];
      E = Energy ();
      x[i] = 1 - x[i];
      if (E < E0) return 0;

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

   // Allocate array space; these are global variables.
   V     = Array (50, 50);
   x     = (int *) calloc (51, sizeof (int));
   x_min = (int *) calloc (51, sizeof (int));

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


