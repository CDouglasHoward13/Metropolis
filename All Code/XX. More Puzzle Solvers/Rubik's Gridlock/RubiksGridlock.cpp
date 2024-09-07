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
// Metropolis "Rubik's Gridlock" puzzle solver.
////////////////////////////////////////////////////////////////////////////////

// Global variables.
int *x0, *y0, *x1, *y1, *or, **grid;
int I, x0_cur, y0_cur, or_cur;
double *P, T;
// The puzzle's three white "puzzle rectangles" are 1x1, 1x2, and 1x3 (numbers 9, 10, 11).
// The "solution rectangles" are 1x4, 1x5, 2x2, 2x3, 2x4, 2x5, 3x3, and 3x4 (numbers 1 - 8).
// Width-1 and length-1 of the eleven rectangles:
int dim[2][12] = { {0,    0, 0, 1, 1, 1, 1, 2, 2,    0, 0, 0},   // Width-1.
                   {0,    3, 4, 1, 2, 3, 4, 2, 3,    0, 1, 2} }; // Length-1.

// Functions found below.
void GetPuzzle ();
int  Energy ();
void Proposal ();
void Restore ();
int  Defective (int, int, int, int);
void Counts ();
void Metropolis ();
void Report ();
void AllocateMemory ();

#include "MetropolisFunctions.h"

////////////////////////////////////////////////////////////////////////////////
// Main program. ///////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
int main () {

   // Initialize certain quantities, get the puzzle, allocate memory,
   //   compute acceptance probabilities, generate initial configuration,
   //   Report the size of the state space and the number of neighbors.
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
         Restore ();
      }

      // Periodically report progress.
      if (n % 10000000 == 0) printf ("%10d %2d\n", n, E);

   }

   printf ("Solved after %d steps of the Markov chain in %.3f seconds.\n", n, Time ());

}

////////////////////////////////////////////////////////////////////////////////
// Propose a random change to the configuration. ///////////////////////////////
////////////////////////////////////////////////////////////////////////////////
void Proposal () {

   int x0_new, y0_new, x1_new, y1_new, or_new;

   // Work until a suitable proposal is found.
   while (1) {

      while (1) {
         // Choose a random puzzle piece to move.
         I = RandomInteger (1,8);
         // Choose its new location on the grid.
         x0_new = RandomInteger(1,8);
         y0_new = RandomInteger(1,8);
         or_new = (MTUniform() <= 0.5);
         // Make sure it's different from its current location.
         if (x0_new != x0[I] || y0_new != y0[I] || or_new != or[I]) break;
      }

      // Calculate upper-right corner of new rectangle location.
      x1_new = x0_new + dim[or_new][I];
      y1_new = y0_new + dim[1-or_new][I];

      // A proposal that covers a puzzle rectangle site or protrudes beyond the
      //    8x8 grid is defective.
      if (!Defective (x0_new, y0_new, x1_new, y1_new)) break;

   }

   // Record current lower-left coordinates and orientation in case they need to be restored.
   // These (and I) are global variables.
   x0_cur  = x0[I];
   y0_cur  = y0[I];
   or_cur  = or[I];

   // Switch to proposed new coordinates and orientation.
   x0[I]  = x0_new;
   y0[I]  = y0_new;
   x1[I]  = x1_new;
   y1[I]  = y1_new;
   or[I]  = or_new;

}

////////////////////////////////////////////////////////////////////////////////
// Restore current configuration.                               ////////////////
////////////////////////////////////////////////////////////////////////////////
void Restore () {

   // Restore.
   x0[I]  = x0_cur;
   y0[I]  = y0_cur;
   or[I]  = or_cur;

   // Recompute upper-right corner.
   x1[I]  = x0[I] + dim[or[I]][I];
   y1[I]  = y0[I] + dim[1-or[I]][I];

}

////////////////////////////////////////////////////////////////////////////////
// Calculate the energy of the current configuration. //////////////////////////
////////////////////////////////////////////////////////////////////////////////
int Energy () {

   int i, x, y, E;

   // Initialize grid coverage counts.
   for (x = 1; x <= 8; x++) {
      for (y = 1; y <= 8; y++) {
         grid[x][y] = 0;
      }
   }

   // Loop through all puzzle pieces.
   for (i = 1; i <= 11; i++) {
      // See which grid sites rectangle number i covers.
      for (x = x0[i]; x <= x1[i]; x++) {
         for (y = y0[i]; y <= y1[i]; y++) {
            grid[x][y] ++;
         }
      }
   }

   // Penalize deviations from one.
   E = 0;
   for (x = 1; x <= 8; x++) {
      for (y = 1; y <= 8; y++) {
         E += abs (grid[x][y] - 1);
      }
   }

   return E;

}

////////////////////////////////////////////////////////////////////////////////
// This function reads in the locations of the three white puzzle rectangles.
// The file format is:
/*
2 1 h
5 1 v
4 6 h
*/
// This means: the 1x1 piece is at site (2,1); the 1x2 piece's lower-left part
//   is at site (5,1) and is oriented vertically; and the 1x3 piece's lower-left part
//   is at site (4,6) and is oriented horizontally. This is puzzle 01.txt.

// An initial configuration is generated.

// The number of states and neighbors is counted (just out of curiosity).

// Proposal acceptance probabilities are pre-calculated for efficiency.
////////////////////////////////////////////////////////////////////////////////
void GetPuzzle () {

   int i, or0, dE;
   char input[100], c;
   FILE *fp;

   // First, seed the RNG.
   MTUniform ();

   // Now allocate array space.
   AllocateMemory ();

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

   // Read in the clue data.
   for (i = 9; i <= 11; i++) {
      fgets (input, 99, fp);
      sscanf (input, "%d %d %c", x0+i, y0+i, &c);
      // Determine the orientation --- horizontal (1) or vertical (0).
      or0 = (c == 'h' || c == 'H');
      x1[i] = x0[i] + dim[or0][i];
      y1[i] = y0[i] + dim[1-or0][i];
   }
   fclose (fp);

   // Now generate a random initial non-defective configuration.
   for (i = 1; i <= 8; i++) {
      // Work until a "non-defective" position is found for solution rectangle number i.
      while (1) {
         x0[i] = RandomInteger(1,8); // Lower-left corner.
         y0[i] = RandomInteger(1,8);
         or[i] = (MTUniform() <= 0.5); // Horizontal/vertical orientation.
         x1[i] = x0[i] + dim[or[i]][i]; // Upper-right corner.
         y1[i] = y0[i] + dim[1-or[i]][i];
         if (!Defective (x0[i], y0[i], x1[i], y1[i])) break;
      }
   }

   // Report how many states and neighbors there are for this puzzle.
   Counts ();

   // Now get the temperature and calculate acceptance probabilities.
   T = GetDouble ("What is the temperature parameter (best = 1.04)?... ");
   if (T > 0.0) {
      for (dE = 1; dE <= 25; dE++) {
         P[dE] = exp (-dE/T);
      }
   }
   // P[dE] = 0 for 25 < dE < 200, per AllocateMemory().


}

////////////////////////////////////////////////////////////////////////////////
// Determine if a rectangle with corners at (a,b) and (c,d) covers a
//   puzzle rectangle square or lies partly off the 8x8 grid.
////////////////////////////////////////////////////////////////////////////////
int Defective (int a, int b, int c, int d) {

   int j, x, y;

   // Does proposed rectangle location extend off the grid?
   if (c > 8 || d > 8) return 1;

   // Does proposed rectangle location cover a puzzle rectangle site?
   for (j = 9; j <= 11; j++) {
      for (x = x0[j]; x <= x1[j]; x++) {
         for (y = y0[j]; y <= y1[j]; y++) {
            // (x,y) is a puzzle rectangle site; does it lie in the proposed rectangle?
            if (a <= x && x <= c && b <= y && y <= d) return 1;
         }
      }
   }

   // If all is good, it's not defective.
   return 0;

}

////////////////////////////////////////////////////////////////////////////////
// Count how many states there are and the number of neighbors of each state.
////////////////////////////////////////////////////////////////////////////////
void Counts () {

   int i, x, y, z, x1, y1, c, neighbors, a;
   double states, b;

   neighbors = 0;
   states = 1;

   // Loop through the eight solution rectangles.
   for (i = 1; i <= 8; i++) {
      // Count how many non-defective positions rectangle number i has.
      c = 0;
      // Loop through all possible values of (x,y,z).
      for (x = 1; x <= 8; x++) {
         for (y = 1; y <= 8; y++) {
            for (z = 0; z <= 1; z++) {
               // (x1,y1) is the upper-right corner of rectangle number i in position (x,y,z).
               x1 = x + dim[z][i];
               y1 = y + dim[1-z][i];
               c += !Defective (x, y, x1, y1);
            }
         }
      }

      neighbors += c-1; // A state is not a neighbor of itself.
      states *= c;
   }

   // Convert the number of states to scientific notation.
   a = log(states)/log(10);
   b = states/pow(10,a);

   printf ("This puzzle has %.2f x 10^%d states, each of which has %d neighbors.\n",
            b, a, neighbors);

   return;

}

////////////////////////////////////////////////////////////////////////////////
// Create output files for viewing with RG.tex. ////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
void Report () {

   int i, m, n;
   double x, y, u, v;
   FILE *fp;
   static int a=0;

   // PUZZLE.
   fp = fopen ("PuzzleClues.txt", "w");
   // Loop through the three puzzle rectangles.
   for (i = 9; i <= 11; i++) {
      // Compute the coordinates of rectangle i's corners.
      x = x0[i] - 1.0;
      y = y0[i] - 1.0;
      u = x1[i];
      v = y1[i];
      fprintf (fp, "\\plot %.0f %.0f  %.0f %.0f  %.0f %.0f  %.0f %.0f  %.0f %.0f /\n",
                   x,y, u,y, u,v, x,v, x,y);
      for (m = x0[i]; m <= x1[i]; m++) for (n = y0[i]; n <= y1[i]; n++) {
         fprintf (fp, "\\put {$\\bullet$} at %.1f %.1f\n", m-0.5, n-0.5);
      }   
   }
   fclose (fp);

   if (a == 0) {
      // Empty out solution file.
      fp = fopen ("PuzzleSolution.txt", "w");
      fclose (fp);
      printf ("\n");
      printf ("View the puzzle without solution if you wish by Plain Texing RG.tex, or...");
      Pause ();
      // The next call to Report() is for the solution.
      a = 1;
      return;
   }

   // SOLUTION.
   fp = fopen ("PuzzleSolution.txt", "w");
   // Loop through all eleven rectangles.
   for (i = 1; i <= 11; i++) {
      // Compute the coordinates of rectangle i's corners.
      x = x0[i] - 1.0;
      y = y0[i] - 1.0;
      u = x1[i];
      v = y1[i];
      fprintf (fp, "\\plot %.0f %.0f  %.0f %.0f  %.0f %.0f  %.0f %.0f  %.0f %.0f /\n",
                   x,y, u,y, u,v, x,v, x,y);
   }
   // Highlight the puzzle rectangles.
   for (i = 9; i <= 11; i++) {
      for (m = x0[i]; m <= x1[i]; m++) for (n = y0[i]; n <= y1[i]; n++) {
         fprintf (fp, "\\put {$\\bullet$} at %.1f %.1f\n", m-0.5, n-0.5);
      }
   }      

   fclose (fp);

   printf ("View the puzzle and solution using Plain TeX with RG.tex.\n");
   Exit();

}

////////////////////////////////////////////////////////////////////////////////
// Allocate more than enough array space. //////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
void AllocateMemory () {

   int i;

   x0  = (int *) calloc (20, sizeof (int));
   y0  = (int *) calloc (20, sizeof (int));
   x1  = (int *) calloc (20, sizeof (int));
   y1  = (int *) calloc (20, sizeof (int));
   or  = (int *) calloc (20, sizeof (int));

   grid = (int **) calloc (20, sizeof (int *));
   for (i = 1; i < 20; i++) {
      grid[i] = (int *) calloc (20, sizeof (int));
   }

   P = (double *) calloc (200, sizeof (double));

}


