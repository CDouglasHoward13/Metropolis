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
// Metropolis "Not Alone" puzzle solver. Not Alone appears in the New
// York Times Magazine under Prasanna Seshadri's byline.
////////////////////////////////////////////////////////////////////////////////

// Global variables.
int **x, **y, **clue, N=8, halfN=4, R, C1, C2; // The puzzle is now 8x8.
double *P, T;

// Functions found below.
void GetPuzzle ();
int  Energy ();
void Proposal ();
void Metropolis ();
void Report ();
void AllocateMemory ();

#include "MetropolisFunctions.h"

////////////////////////////////////////////////////////////////////////////////
// Main program. ///////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
int main () {

   // Get the puzzle, allocate memory, compute acceptance probabilities,
   // generate initial configuration.
   GetPuzzle ();

   // Report the puzzle.
   Report ();

   // Solve the puzzle via Metropolis.
   Metropolis ();

   // Generate output files and exit program.
   Report ();

}

////////////////////////////////////////////////////////////////////////////////
// Implement the Metropolis algorithm. /////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
void Metropolis () {

   int E, deltaE, AcceptTransition, n = 0;

   // Calculate the energy for the initial configuration.
   E = Energy ();

   // Time the computations.
   Time ();

   // Keep going until the solution is found.
   while (E > 0) {

      // Count Markov chain steps.
      n ++;

      // Propose a random change to the current configuration.
      Proposal ();

      // Compute the resulting change in energy.
      deltaE = Energy () - E;

      // Start with the zero-temperature dynamics.
      AcceptTransition = 0;

      // If deltaE <= 0, accept the transition.
      if (deltaE <= 0) {
         AcceptTransition = 1;
      }

      // If T > 0, accept with small probability P[deltaE] as previously calculated.
      else if (T > 0) {
         if (MTUniform() <= P[deltaE]) {
            AcceptTransition = 1;
         }
      }

      // If accepted, keep the change and update the energy.
      if (AcceptTransition) {
         E += deltaE;
      }

      // Otherwise change back to previous configuration.
      else {
         x[R][C1] = 1 - x[R][C1];
         x[R][C2] = 1 - x[R][C2];
      }

      // Periodically report progress.
      if (n % 1000000 == 0) printf ("%10d %2d\n", n, E);

   }

   printf ("Solved after %d steps of the Markov chain in %.3f seconds.\n", n, Time ());

}

////////////////////////////////////////////////////////////////////////////////
// Propose a random change to the configuration. ///////////////////////////////
////////////////////////////////////////////////////////////////////////////////
void Proposal () {

   // R, C1, and C2 are global variables.

   // Choose a row to change.
   R  = RandomInteger (1,N);

   // Choose a random circle in that row.
   C1 = RandomInteger (1,N);

   // Choose another random circle of a different color in that row.
   while (1) {
      C2 = RandomInteger (1,N);
      if (x[R][C1] != x[R][C2]) break;
   }

   // Switch the colors. This takes 0 to 1, and 1 to 0.
   x[R][C1] = 1 - x[R][C1];
   x[R][C2] = 1 - x[R][C2];

   return;

}

////////////////////////////////////////////////////////////////////////////////
// Calculate the energy of the current configuration col[][]. //////////////////
////////////////////////////////////////////////////////////////////////////////
int Energy () {

   int r, c, count, E=0;

   // Deviation from clues.
   for (r = 1; r <= N; r++) {
      for (c = 1; c <= N; c++) {
         if (clue[r][c] && (x[r][c] != y[r][c])) E += 5;
      }
   }

   // Columns deviation from halfN black circle count.
   for (c = 1; c <= N; c++) {
      count = 0;
      for (r = 1; r <= N; r++) {
         count += x[r][c];
      }
      E += abs(count-halfN);
   }

   // The row black circle counts are always halfN.

   // Now count the circles that are alone horizontally.
   for (r = 1; r <= N; r++) {
      for (c = 2; c <= N-1; c++) {
         E += ((x[r][c-1] != x[r][c]) && (x[r][c+1] != x[r][c]));
      }
   }

   // Now count the circles that are alone vertically.
   for (c = 1; c <= N; c++) {
      for (r = 2; r <= N-1; r++) {
         E += ((x[r-1][c] != x[r][c]) && (x[r+1][c] != x[r][c]));
      }
   }

   return E;

}

////////////////////////////////////////////////////////////////////////////////
// This function reads in the NxN puzzle's cluess.
// Array space is allocated.
// An initial configuration is generated.
// Proposal acceptance probabilities are pre-calculated.
////////////////////////////////////////////////////////////////////////////////
void GetPuzzle () {

   int i, r, c, dE, nclues, y0;
   char input[100];
   FILE *fp;

   // First, seed the RNG.
   MTUniform ();

   printf ("Please input the name of the puzzle input file... ");
   fgets (input, 99, stdin);

   // Now terminate the file name with ".txt".
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

   // Read in the number of clues.
   fgets (input, 99, fp);
   sscanf (input, "%d", &nclues);

   // Allocate array space.
   AllocateMemory ();

   // Read in the clue data.
   for (i = 1; i <= nclues; i++) {
      fgets (input, 99, fp);
      sscanf (input, "%d %d %d", &r, &c, &y0); // Row, column, color.
      y[r][c] = y0;
      clue[r][c] = 1;
   }
   // Now clue[r][c] if 1 (0) if there is (isn't) a clue at site (r,c).
   // If there is a clue there, y[r][c] is 1 (0) if the circle there is black (white).
   fclose (fp);

   // Now generate an initial configuration.
   for (r = 1; r <= N; r++) {
      for (c = 1; c <= halfN; c++) {
         x[r][c] = 1;
      }
   }

   // Now get the temperature and calculate acceptance probabilities.
   T = GetDouble ("What is the temperature parameter (0.25 seems good)?... ");
   if (T > 0.0) {
      for (dE = 1; dE < 25; dE++) {
         P[dE] = exp (-dE/T);
      }
   }

}

////////////////////////////////////////////////////////////////////////////////
// Create output files. ////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
void Report () {

   int r, c;
   double u, v;
   FILE *fp;
   static int a=0;

   // PUZZLE.
   fp = fopen ("PuzzleClues.txt", "w");
   // Loop through all the squares.
   for (r = 1; r <= N; r++) for (c = 1; c <= N; c++) {
      // Compute the coordinates of square (r,c)'s center.
      u = c - 0.5;
      v = N+0.5 - r;
      if (clue[r][c]) {
         if (y[r][c]) {
            fprintf (fp, "\\put {$\\bullet$} at %f %f\n", u, v);
         } else {
            fprintf (fp, "\\put {$\\circ$} at %f %f\n", u, v);
         }
      }
   }
   fclose (fp);

   if (a == 0) {
      // Empty out solution file.
      fp = fopen ("PuzzleSolution.txt", "w");
      fclose (fp);
      // Next call to Report() is for the solution.
      a = 1;
      return;
   }

   // SOLUTION.
   fp = fopen ("PuzzleSolution.txt", "w");
   // Loop through all the squares.
   for (r = 1; r <= N; r++) for (c = 1; c <= N; c++) {
      // Compute the coordinates of square (r,c)'s center.
      u = c - 0.5;
      v = N+0.5 - r;
      if (x[r][c]) {
         fprintf (fp, "\\put {$\\bullet$} at %f %f\n", u, v);
      } else {
         fprintf (fp, "\\put {$\\circ$} at %f %f\n", u, v);
      }
   }
   fclose (fp);

   printf ("View the puzzle and solution using Plain TeX with NA.tex.\n");
   Exit();

}

////////////////////////////////////////////////////////////////////////////////
// Allocate array space. ///////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
void AllocateMemory () {

   int r;

   x    = (int **) calloc (N+1, sizeof (int *));
   y    = (int **) calloc (N+1, sizeof (int *));
   clue = (int **) calloc (N+1, sizeof (int *));

   for (r = 1; r <= N; r++) {
      x[r]    = (int *) calloc (N+1, sizeof (int));
      y[r]    = (int *) calloc (N+1, sizeof (int));
      clue[r] = (int *) calloc (N+1, sizeof (int));
   }

   P = (double *) calloc (25, sizeof (double));

}


