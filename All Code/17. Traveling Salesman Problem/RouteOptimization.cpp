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
// This code looks for near-optimal TSP routes using the Metropolis
// algorithm as described in Section 17. The sites toured are drill sites
// on a circuit board.
// Look at www.math.uwaterloo.ca/tsp for images and information on TSP.
////////////////////////////////////////////////////////////////////////////////

// Global variables.
int K, n_min, *c, *best, *temp, i0, j0;
double *X, *Y, **d, E, E_min;

// These functions are found below.
void    InitializeArrays ();
void    RandomRoute ();
void    ReportRoute ();
void    Proposal ();
void    Reverse ();
void    Metropolis ();

// These functions are in common to all applications.
#include "MetropolisFunctions.h"

////////////////////////////////////////////////////////////////////////////////
// Main program.
////////////////////////////////////////////////////////////////////////////////
int main() {

   // Generate a random initial route through the sites.
   RandomRoute ();

   // Find a low energy (short) route via Metropolis.
   Metropolis ();

   // Pause, then exit program.
   Exit ();

}

////////////////////////////////////////////////////////////////////////////////
// Minimize route length via Matropolis.
////////////////////////////////////////////////////////////////////////////////
void Metropolis () {

   double T, DeltaE, p, U, t, t1;
   int k, n, AcceptTransition, NextReport;



   // Report next route at Markov chain period number "NextReport".
   // Reports are generated at 1000, 10000, 100000, 1000000, 10000000, and
   //    100000000 steps of the Markov chain.
   NextReport = 1000;

   // Get the temperature parameter.
   T = GetDouble ("\nWhat is the temperature (best is .07)?... ");

   printf ("\nI'll be done in 60 seconds. ");
   t = t1 = Time ();
   n = 0;

   // Run the Markov chain for 60 seconds.
   while (t < 60.0) {

      // Every five seconds indicate that it's still thinking.
      t = Time ();
      if (t > t1 + 5.0) {
         printf (". ");
         t1 = t;
      }

      // Update the Markov chain step counter.
      n ++;

      // Get a proposed random change.
      Proposal ();

      // Compute the change in energy associated with reversing that portion
      //   of the route.
      DeltaE =  d [c[i0-1]] [c[j0]] + d [c[i0]] [c[j0+1]]
              - d [c[i0-1]] [c[i0]] - d [c[j0]] [c[j0+1]];

      // See if the proposed transition is accepted.
      AcceptTransition = 0;

      if (DeltaE <= 0) {
         AcceptTransition = 1;
      }

      else if (T > 0) {
         p = exp (-DeltaE / T);
         U = MTUniform ();
         if (U <= p) {
            AcceptTransition = 1;
         }
      }

      // Effect the transition. If the new energy is a best-so-far, record the
      //    route (in best[*]) and the new minimal energy (Emin).
      if (AcceptTransition) {

         // Reverse the part of route "c" from i0 to j0.
         Reverse ();

         // Update the length of the current route (the route's "energy").
         E += DeltaE;

         // Record data for the best-route-so-far, if appropriate.
         if (E < E_min) {
            E_min = E;
            n_min = n;
            for (k = 1; k <= K+1; k++) {
               best[k] = c[k];
            }
         }

      } // End of "if" statement.

      // Periodically report the route.
      if (n == NextReport) {

         // Report the current route for viewing with TeX software.
         ReportRoute ();

         // Update the next period to be reported.
         NextReport *= 10;

      }

   } // This ends the Markov chain simulation loop.

   // Report the best route found throughout the Markov chain.
   E = E_min;
   ReportRoute ();

   // Finish up; report best-found route length to the screen.
   printf ("\n\n");
   printf ("%.1f million Markov chain steps completed in 60 seconds.\n\n", n/1000000.0);
   printf ("Shortest route was number %d with length %.3f\n\n", n_min, E_min);
   printf ("View the solution with ShowRoutes.tex using Plain TeX.\n");

}

////////////////////////////////////////////////////////////////////////////////
// Determine the proposed transition.
////////////////////////////////////////////////////////////////////////////////
void Proposal () {

   int k;

   // Randomly choose a "neighbor" of the current route; use accept/reject.
   // First pick a pair (i0,j0), uniformly from {2,...,K} x {2,...,K} with
   //   (i) i0 < j0 and (ii) (i0,j0) != (2,K).  Condition (ii) prevents
   //   a simple reversal of direction.
   while (1) {

      // Pick i and j independently and uniformly from {2,...,K}.
      i0 = RandomInteger (2, K);
      j0 = RandomInteger (2, K);

      // Now switch i0 and j0, if necessary, so that i0 <= j0.
      if (j0 < i0) {
         k = j0;
         j0 = i0;
         i0 = k;
      }

      // See if they are acceptable, i.e, if they satisfy (i) and (ii) above.
      k = j0 - i0;
      if (0 < k && k < K-2) {
         break;
      }

   }


}

////////////////////////////////////////////////////////////////////////////////
// This function reverses the part of route c from i0 to j0.
////////////////////////////////////////////////////////////////////////////////
void Reverse () {

   int k;

   // Make a copy of the part of the route to be reversed.
   for (k = i0; k <= j0; k++) {
      temp[k] = c[k];
   }

   // Now reverse it. Observe that when k = i0 we get c[i0] = temp[j0],
   //    and when k = j0 we get c[j0] = temp[i0].
   for (k = i0; k <= j0; k++) {
      c[k] = temp[i0+j0-k];
   }

   return;

}

////////////////////////////////////////////////////////////////////////////////
// Allocate array space for "K" sites and specify X and Y coordinates of the
//   drill holes.
////////////////////////////////////////////////////////////////////////////////
void InitializeArrays () {

   int i, j;
   double dx, dy;

   // X and Y coordinates in millimeters of the 183 drill sites.
   double X0[] =
   {0, 5,   5,   7,   7,   8,   8,   8,  12,  12,  12,  12,  12,  11,  16,  16,
      16,  18,  18,  18,  18,  20,  20,  20,  20,  23,  23,  23,  23,  23,  22,
      28,  28,  28,  27,  26,  26,  25,  25,  33,  32,  30,  28,  31,  30,  30,
      34,  34,  34,  34,  36,  36,  36,  36,  38,  40,  40,  40,  40,  40,  40,
      41,  41,  43,  43,  43,  44,  44,  44,  46,  46,  46,  49,  48,  48,  48,
      48,  52,  52,  52,  52,  52,  52,  52,  52,  52,  52,  52,  54,  54,  54,
      56,  56,  56,  56,  56,  56,  57,  57,  61,  61,  61,  68,  68,  72,  72,
      72,  72,  71,  73,  73,  76,  76,  76,  78,  78,  77,  78,  80,  80,  79,
      79,  82,  82,  81,  83,  86,  86,  86,  89,  89,  89,  92,  92,  92,  84,
      94,  94,  94,  94,  90,  93, 100, 100, 102, 102, 105, 105, 105, 105, 105,
     105, 110, 110, 113, 113, 113, 117, 117, 117, 117, 121, 121, 121, 121, 121,
     122, 124, 124, 125, 129, 129, 130, 130, 135, 137, 137, 137, 136, 142, 142,
     146, 148, 148};

   double Y0[] =
   {0, 20,  35,  32,  42,  13,  19,  25,  12,  18,  26,  30,  33,  40,  30,  35,
       39,  15,  22,  38,  42,   9,  12,  26,  35,  25,  28,  31,  33,  37,  42,
       10,  17,  20,  25,  32,  36,  39,  43,   9,  33,  33,  33,  35,  38,  43,
       23,  35,  38,  43,  15,  19,  28,  33,  22,   7,  11,  36,  39,  42,  46,
       20,  33,   5,  37,  42,  13,  23,  33,  22,  26,  44,   5,   8,  32,  36,
       40,   5,  11,  18,  23,  26,  30,  33,  36,  38,  41,  45,  35,  40,  43,
        7,  12,  17,  22,  27,  32,  37,  42,  36,  39,  45,  36,  39,  13,  18,
       23,  28,  32,  35,  38,  10,  22,  33,  24,  29,  35,  37,  33,  37,  10,
       21,  10,  21,  24,  41,  25,  36,  40,  25,  32,  40,  25,  32,  40,  35,
       10,  18,  32,  40,  35,  35,  30,  39,  21,  35,   2,  16,  24,  28,  36,
       41,  14,  22,  28,  36,  41,   7,  32,  35,  40,   3,  15,  24,  27,  36,
       39,  17,  31,  12,  13,  18,  23,  38,  13,  13,  24,  30,  40,  24,  30,
       43,  20,  30};

   // The number of sites, here K = 183.
   K = 183;

   // Allocate additional necessary array space.
   X = (double *) calloc (K+1, sizeof (double));
   Y = (double *) calloc (K+1, sizeof (double));
   d = (double **) calloc (K+1, sizeof (double *));
   for (i = 1; i <= K; i++) {
      d[i] = (double *) calloc (K+1, sizeof (double));
   }
   c    = (int *) calloc (K+3, sizeof (int));
   best = (int *) calloc (K+3, sizeof (int));
   temp = (int *) calloc (K+3, sizeof (int));

   // Copy the above coordinates to the global variables. Perturb them slightly
   //   to avoid distance ties. (The coordinates are currently integer-valued.)
   for (i = 1; i <= K; i++) {
      X[i] = X0[i] + 0.001 * MTUniform();
      Y[i] = Y0[i] + 0.001 * MTUniform();
   }

   // Compute the distance between each pair of sites in centimeters (the data
   //    is in millimeters, so divide by 10).
   for (i = 1; i <= K; i++) {
      for (j = 1; j <= K; j++) {
         dx = (X[i] - X[j]) / 10.0;
         dy = (Y[i] - Y[j]) / 10.0;
         d[i][j] = sqrt(dx*dx + dy*dy);
      }
   }


   return;

}

////////////////////////////////////////////////////////////////////////////////
// Randomly select the initial route, starting and ending at site 1.
////////////////////////////////////////////////////////////////////////////////
void RandomRoute () {

   int i, j, k;

   printf ("I'm looking for the minimal route through a circuit board.\n");

   // Seed the RNG.
   MTUniform();
   
   // Allocate array space for 183 sites and specify site coordinates.
   InitializeArrays ();

   // Initially tour them in numerical order.
   for (i = 1; i <= K; i++) {
      c[i] = i;
   }

   // End the tour at site 1, where it begain.
   c[K+1] = 1;

   // Now randomly permute sites 2,3,...,K.
   for (i = 2; i < K; i++) {
      j = RandomInteger (i,K);
      k = c[i];
      c[i] = c[j];
      c[j] = k;
   }

   // Compute the initial route distance.
   E = 0;
   for (i = 1; i <= K; i++) {
      E += d [c[i]] [c[i+1]];
   }

   // Initialize the minimal energy and where it occurs in the Markov chain.
   E_min = E;
   n_min = 0;

   // Record the initial random route as the best-so-far.
   for (i = 1; i <= K+1; i++) {
      best[i] = c[i];
   }

   // Report the initial route.
   ReportRoute ();

   return;

}

////////////////////////////////////////////////////////////////////////////////
// Report data for viewing with TeX software. 
////////////////////////////////////////////////////////////////////////////////
void ReportRoute () {

   static int n = 0;
   int i, *r;
   char filename[10][100] = {"RandomRoute.txt", "1000.txt", "10000.txt",
                             "100000.txt", "1000000.txt", "10000000.txt",
                             "100000000.txt", "BestRoute.txt"};
   int N[10] = {0, 1000, 10000, 100000, 1000000, 10000000, 100000000};
   FILE *fp;

   // Open the appropriate output file.
   fp = fopen (filename[n], "w");

   // The last root displayed is the best route; earlier routes are the current route.
   r = (n < 7 ? c : best);

   // Header for optimization in progress -- current route.
   if (n < 7) {
      fprintf (fp, "\\put {\\sl Route distance = %8.3f cm} at 80 58\n", E);
      fprintf (fp, "\\put {\\sl $n$ = %9d} at 80 53\n", N[n]);
   }

   // Header for best-found route.
   else {
      fprintf (fp, "\\put {\\sl Best-found route distance = %8.3f cm} at 80 58\n", E);
   }

   // Display the route.
   for (i = 1; i <= K; i++) {
      fprintf (fp, "\\plot %8.4f %8.4f  %8.4f %8.4f /\n",
                 X[r[i]], Y[r[i]], X[r[i+1]], Y[r[i+1]]);
   }
   fclose (fp);

   // Report the site coordinates to another output file the first time through.
   if (n == 0) {
      fp = fopen ("Sites.txt", "w");
      for (i = 1; i <= K; i++) {
         fprintf (fp, "%8.3f %8.3f\n", X[i], Y[i]);
      }
      fclose (fp);
   }

   n ++;

   return;

}


