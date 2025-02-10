////////////////////////////////////////////////////////////////////////////////
// This code solves "Traffic Signals" puzzles by Prasanna Seshadri via
//    Metropolis. The puzzle lives on an 8x8 grid. (S=8, below.)
////////////////////////////////////////////////////////////////////////////////

// Global variables.
int S=8, K, i0, j0, R, G;
int *X, *Y, **d, *c, *temp, *R0, *R1, *G0, *G1;
double x0, x1, y0, y1;

// These functions are found below.
void    InitializeArrays ();
void    GetPuzzle ();
void    Metropolis ();
void    ReportRoute ();
void    Proposal ();
void    Reverse ();
int     Energy ();
void    DualSegment ();

// These functions are in common to all applications.
#include "MetropolisFunctions.h"

////////////////////////////////////////////////////////////////////////////////
// Main program.
////////////////////////////////////////////////////////////////////////////////
int main() {

   // Getting started.
   InitializeArrays ();

   // Get the puzzle data.
   GetPuzzle ();

   // Find the optimal route via Metropolis.
   Metropolis ();

   // Report the solution.
   ReportRoute ();

   // Pause, then exit program.
   Exit ();

}

////////////////////////////////////////////////////////////////////////////////
// Find puzzle solution via Metropolis.
////////////////////////////////////////////////////////////////////////////////
void Metropolis () {

   double T, p, U;
   int n, DeltaE, E, AcceptTransition;

   // Starting energy.
   E = Energy ();

   // Get the temperature parameter.
   T = GetDouble ("What is the temperature (.2 seems good)?... ");

   // Markov chain counter.
   n = 0;

   // Time the calculations.
   Time ();

   // Run until energy is minimal (route length=K with no violations).
   while (E > K) {

      // Update the Markov chain step counter.
      n ++;

      // Periodically report the energy to the screen. The puzzle should
      //    be solved before this happens.
      if (n % 5000000 == 0) printf ("%9d %3d\n", n, E);

      // Get a proposed random change to the route.
      Proposal ();

      // Compute the change in energy.
      DeltaE = Energy() - E;

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

      // Effect the transition -- keep the new route and update the energy.
      if (AcceptTransition) {
         E += DeltaE;
      }

      // Reinstate previous path if proposal is not accepted.
      else {
         Reverse ();
      }   

   } // This ends the Markov chain simulation loop.

   printf ("Solved in %.1f seconds.  View the solution with TS.tex.", Time());

}

////////////////////////////////////////////////////////////////////////////////
// Determine the proposed transition.
////////////////////////////////////////////////////////////////////////////////
void Proposal () {

   int k;

   // Randomly choose a "neighbor" of the current route; use accept/reject.
   // First pick a pair (i0,j0), uniformly from {2,...,K} x {2,...,K} with
   //   i0 < j0.

   i0 = j0 = 0;
   while (i0 == j0) {
      i0 = RandomInteger (2, K);
      j0 = RandomInteger (2, K);
   }

   // Now switch i0 and j0, if necessary, so that i0 <= j0.
   if (j0 < i0) {
      k = j0;
      j0 = i0;
      i0 = k;
   }

   // Reverse that part of the path.
   Reverse ();

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
// Determine current energy.
////////////////////////////////////////////////////////////////////////////////
int Energy () {

   int E, V, r, g, i;

   // Compute the route distance.
   E = 0;
   for (i = 1; i <= K; i++) {
      E += d [c[i]] [c[i+1]];
   }

   // Add traffic signal violations. Initially assume all green lights
   //    are violated and no red lights are violated.
   V = G;
   for (i = 1; i <= K; i++) if (d [c[i]] [c[i+1]] == 1) {

      // Red light violations augment V.
      for (r = 1; r <= R; r++) {
         if (   (c[i] == R0[r] && c[i+1] == R1[r])
             || (c[i] == R1[r] && c[i+1] == R0[r])  ) V++;
      }

      // Green light non-violations decrement V.
      for (g = 1; g <= G; g++) {
         if (   (c[i] == G0[g] && c[i+1] == G1[g])
             || (c[i] == G1[g] && c[i+1] == G0[g])  ) V--;
      }

   }

   // Energy is route length plus traffic signal violations.
   E += V;

   return E;

}

////////////////////////////////////////////////////////////////////////////////
// Allocate array space for "K" cells and specify X and Y coordinates of the
//   puzzle cells.
////////////////////////////////////////////////////////////////////////////////
void InitializeArrays () {

   int i, j;

   // The grid is SxS, so the number of cells is...
   K = S*S;

   // Allocate additional necessary array space.
   X = (int *) calloc (K+1, sizeof (int));
   Y = (int *) calloc (K+1, sizeof (int));
   d = (int **) calloc (K+1, sizeof (int *));
   for (i = 1; i <= K; i++) {
      d[i] = (int *) calloc (K+1, sizeof (int));
   }
   c    = (int *) calloc (K+3, sizeof (int));
   temp = (int *) calloc (K+3, sizeof (int));
   R0   = (int *) calloc (K+3, sizeof (int));
   R1   = (int *) calloc (K+3, sizeof (int));
   G0   = (int *) calloc (K+3, sizeof (int));
   G1   = (int *) calloc (K+3, sizeof (int));


   // Compute the coodinates of each cell. This is all integer arithmetic.
   for (i = 1; i <= K; i++) {
      Y[i] = (i-1)/S + 1;
      X[i] = i - (Y[i]-1) * S;
      // Re-number from bottom to top.
      Y[i] = S+1 - Y[i];
   }

   // Compute the Manhattan distance between each pair of cells.
   for (i = 1; i <= K; i++) {
      for (j = 1; j <= K; j++) {
         d[i][j] = abs(X[i] - X[j]) + abs(Y[i] - Y[j]);
      }
   }

   // Initially tour them in numerical order.
   for (i = 1; i <= K; i++) {
      c[i] = i;
   }

   // End the tour at cell 1, where it begain.
   c[K+1] = 1;

   // Seed the RNG.
   MTUniform ();

}

////////////////////////////////////////////////////////////////////////////////
// Get the puzzle data.
////////////////////////////////////////////////////////////////////////////////
void GetPuzzle () {

   FILE *fp, *fpr, *fpg;
   int i;
   char input[100];

   // First get the file name.
   printf ("Please input the name of the puzzle input file... ");
   fgets (input, 99, stdin);

   // Now make sure to terminate the file name with ".txt".
   for (i = 1; i <= 20; i++) {
      if (input[i] == '\n' || input[i] == '.') {
         input[i] = '.';
         input[i+1] = 't';
         input[i+2] = 'x';
         input[i+3] = 't';
         input[i+4] = '\0';
         break;
      }
   }

   // Now open the file...
   fp = fopen (input, "r");

   // Open some output files.
   fpr = fopen ("Red.txt", "w");
   fpg = fopen ("Green.txt", "w");

   fgets (input, 99, fp);
   sscanf (input, "%d", &G);
   for (i = 1; i <= G; i++) {
      fgets (input, 99, fp);
      sscanf (input, "%d %d", G0+i, G1+i);
      // Make sure these cells are adjacent.
      if (d [G0[i]] [G1[i]] != 1) {
         printf ("Green error %d %d\n", G0[i], G1[i]);
         Exit ();
      }
      // Compute green line segment for TeX viewing...
      // Endpoints of must-used segment.
      x0 = X[G0[i]];
      x1 = X[G1[i]];
      y0 = Y[G0[i]];
      y1 = Y[G1[i]];
      // Dual graph segment.
      DualSegment ();
      // Report green segment to output file.
      fprintf (fpg, "\\plot %.3f %.3f  %.3f %.3f /\n", x0,y0, x1,y1);
   }

   fgets (input, 99, fp);
   sscanf (input, "%d", &R);
   for (i = 1; i <= R; i++) {
      fgets (input, 99, fp);
      sscanf (input, "%d %d", R0+i, R1+i);
      // Make sure these cells are adjacent.
      if (d [R0[i]] [R1[i]] != 1) {
         printf ("Red error %d %d\n", R0[i], R1[i]);
         Exit ();
      }
      // Compute red line segment for TeX viewing...
      // Endpoints of cannot-use segment.
      x0 = X[R0[i]];
      x1 = X[R1[i]];
      y0 = Y[R0[i]];
      y1 = Y[R1[i]];
      // Dual graph segment.
      DualSegment ();
      // Report red segment to output file.
      fprintf (fpr, "\\plot %.3f %.3f  %.3f %.3f /\n", x0,y0, x1,y1);
   }

   fclose (fp);
   fclose (fpg);
   fclose (fpr);

}

////////////////////////////////////////////////////////////////////////////////
// Rotate the segment (x0,y0),(x1,y1) 90 degrees counter-clockwise about
//    its midpoint.
////////////////////////////////////////////////////////////////////////////////
void DualSegment () {

   double a, b, t;

   // Midpoint of segment is (a,b).
   a = (x0+x1)/2.0;
   b = (y0+y1)/2.0;
   // Segment shifted so it's midpoint is (0,0).
   x0 -= a;
   x1 -= a;
   y0 -= b;
   y1 -= b;
   // Segment rotated 90 degrees counter-clockwise.
   t = x0;
   x0 = -y0;
   y0 = t;
   t = x1;
   x1 = -y1;
   y1 = t;
   // Segment shifted back so midpoint is (a,b)
   x0 += a;
   x1 += a;
   y0 += b;
   y1 += b;

}

////////////////////////////////////////////////////////////////////////////////
// Report solution for viewing with TeX software. 
////////////////////////////////////////////////////////////////////////////////
void ReportRoute () {

   int i;
   FILE *fp;

   // Open the appropriate output file.
   fp = fopen ("Solution.txt", "w");

   // Display the route.
   for (i = 1; i <= K; i++) {
      fprintf (fp, "\\plot %d %d  %d %d /\n",
                 X[c[i]], Y[c[i]], X[c[i+1]], Y[c[i+1]]);
   }
   fclose (fp);

   fp = fopen ("CellCenters.txt", "w");
   for (i = 1; i <= K; i++) {
      fprintf (fp, "%d %d\n", X[c[i]], Y[c[i]]);
   }
   fclose (fp);

   return;

}





