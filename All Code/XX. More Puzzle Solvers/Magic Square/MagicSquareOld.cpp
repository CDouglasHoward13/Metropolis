/*
MIT License

Copyright (c) 2025 CDouglasHoward13

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
// This code generates magic squares. Works well for orders up to 30.  I've
//    gotten it to produce an order 60 square.
////////////////////////////////////////////////////////////////////////////////

// Global variables.
int N, **S;
double *P, T;

// These functions are found below.
void Initialize ();
void Metropolis ();
int  Energy ();
void Report ();

// These functions are in common to all applications.
#include "MetropolisFunctions.h"

////////////////////////////////////////////////////////////////////////////////
// Main program.
////////////////////////////////////////////////////////////////////////////////
int main() {

   // Allocate array space and initialize certain quantities.
   Initialize ();

   // Find a magic square via Metropolis.
   Metropolis ();

   // Report the magic square to the screen and to a LaTeX output file.
   Report ();

   // Pause, then exit program.
   Exit ();

}

////////////////////////////////////////////////////////////////////////////////
// Find an N x N magic square via Metropolis.
////////////////////////////////////////////////////////////////////////////////
void Metropolis () {

   double U, p;
   int DeltaE, E, AcceptTransition, i0, j0, i1, j1, temp, n, n0;

   // Markov chain counter.
   n0 = n = 0;

   // Starting energy.
   E = Energy ();

   // Run until the energy is 0.
   while (E > 0) {

      // Update the Markov chain step counter.
      n0 ++;

      // Periodically report the energy to the screen. 
      if (n0 == 1000000) {
         n++;
         n0 = 0;
         printf ("n = %4d million, E = %3d\n", n, E);
      }   

      // Get a proposed random change to the configuration...

      // First choose a square at random.
      i1 = i0 = RandomInteger (1, N);
      j1 = j0 = RandomInteger (1, N);

      // Now choose a (different) second square.
      while (i1 == i0 && j1 == j0) {
         i1 = RandomInteger (1, N);
         j1 = RandomInteger (1, N);
      }

      // Then swap the numbers in the two squares.
      temp = S[i0][j0];
      S[i0][j0] = S[i1][j1];
      S[i1][j1] = temp;

      // Compute the change in energy.
      DeltaE = Energy() - E;

      // See if the proposed transition is accepted.
      AcceptTransition = 0;

      if (DeltaE <= 0) {
         AcceptTransition = 1;
      }

      else if (T > 0) {
         U = MTUniform ();
         p = (DeltaE < 100 ? P[DeltaE] : 0.0);
         if (U <= p) {
            AcceptTransition = 1;
         }
      }

      // Effect the transition --- keep the proposal and update the energy.
      if (AcceptTransition) {
         E += DeltaE;
      }

      // Reinstate previous configuration if proposal is not accepted.
      else {
         temp = S[i0][j0];
         S[i0][j0] = S[i1][j1];
         S[i1][j1] = temp;
      }

   } // This ends the Markov chain simulation loop.

   printf ("\nFound an order %d magic square.\n", N);
   printf ("View it by processing MS.tex with LaTeX or printing\n");
   printf ("the file MSout.txt.\n");
   printf ("Verify it by running VerifyMS.cpp.\n\n");

   return;

}

////////////////////////////////////////////////////////////////////////////////
// Determine energy of configuration S.
////////////////////////////////////////////////////////////////////////////////
int Energy () {

   int E=0, sum, i, j, rowsum, colsum, d1sum=0, d2sum=0;

   // All the rows and columns and the two diagonals should sum to this.
   sum = (N*(N*N+1))/2;

   // Compute and sum all deviations from this.
   for (i = 1; i <= N; i++) {
      rowsum = colsum = 0;
      for (j = 1; j <= N; j++) {
         rowsum += S[i][j];
         colsum += S[j][i];
      }
      E += abs(rowsum-sum);
      E += abs(colsum-sum);
      d1sum += S[i][i];
      d2sum += S[i][N+1-i];
   }
   E += abs(d1sum-sum);
   E += abs(d2sum-sum);

   return E;

}

////////////////////////////////////////////////////////////////////////////////
// Initialize certain quantities.
////////////////////////////////////////////////////////////////////////////////
void Initialize () {

   int i, j, k, temp, *n;

   // Seed the RNG.
   MTUniform ();

   // Get the square size.
   N = GetInteger ("What is the square's order (N, for an N x N square)?... ");

   // Array space for the square.
   S = (int **) calloc (N+1, sizeof (int *));
   for (i = 1; i <= N; i++) {
      S[i] = (int *) calloc (N+1, sizeof (int));
   }

   // This list holds the numbers 1 through N^2 in numerical order...
   n = (int *) calloc (N*N+1, sizeof (int));
   for (i = 1; i <= N*N; i++) {
      n[i] = i;
   }

   // ... now randomly scramble the list.
   for (i = 1; i < N*N; i++) {
      j = RandomInteger (i, N*N);
      temp = n[i];
      n[i] = n[j];
      n[j] = temp;
   }

   // Place the randomly ordered numbers 1 through N^2 on the N x N grid.
   k = 0;
   for (i = 1; i <= N; i++) {
      for (j = 1; j <= N; j++) {
         k++;
         S[i][j] = n[k];
      }
   }

   // Acceptance probabilities.
   P = (double *) calloc (101, sizeof (double));

   // Pre-calculate acceptance probabilities.
   T = GetDouble ("What is the temperature?... ");
   if (T > 0) {
      for (i = 1; i <= 100; i++) {
         P[i] = exp(-i/T);
      }
   }

   return;

}

////////////////////////////////////////////////////////////////////////////////
// Report the magic square to the screen and to MSout.txt and also to the
// file MS.tex for viewing with LaTeX.
////////////////////////////////////////////////////////////////////////////////
void Report () {

   int i, j;
   FILE *fp;

   // Print the square to the screen and to MSout.txt.
   fp = fopen ("MSout.txt", "w");
   fprintf (fp, "%d x %d magic square found by Metropolis\n", N, N);
   for (i = N; i >= 1; i--) {
      for (j = 1; j <= N; j++) {
         printf (" %5d", S[i][j]);
         fprintf (fp, " %5d", S[i][j]);
      }
      printf ("\n");
      fprintf (fp, "\n");
   }
   fclose (fp);
   printf ("\n\n");

   // Now create a TeX file.
   fp = fopen ("MS.tex", "w");

   fprintf (fp, "\\documentclass[12pt]{article}\\usepackage{pictex}\n");
   fprintf (fp, "\\pagestyle{empty}\\begin{document}\\beginpicture\n");

   fprintf (fp, "\\setcoordinatesystem units <0.3 truein, 0.3 truein>\n");
   fprintf (fp, "\\setplotarea x from 0 to %d, y from  0 to %d\n", N, N);
   fprintf (fp, "\\grid %d %d\n", N, N);
   fprintf (fp, "\\put {\\bf Order %d Magic Square Generated by Metropolis} at %.2f %.2f\n",
                N,  0.5*N, N+0.5);

   for (i = 1; i <= N; i++) {
      for (j = 1; j <= N; j++) {
         // j is the row, i the column.
         fprintf (fp, "\\put {%d} at %.1f %.1f\n", S[i][j], j-0.5, i-0.5);
      }
   }

   fprintf (fp, "\\endpicture\\vfill\\eject\\end{document}\n");

   fclose (fp);

   return;

}







