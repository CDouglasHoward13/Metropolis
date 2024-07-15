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
// Metropolis "One Up" puzzle solver. These are puzzles by Rodolfo
// Kurchan that have appeared in the Sunday New York Times Magazine.
// Rules:
// Between each pair of walls, place the numbers 1 through
// n, where n is the number of squares between the walls. Do
// this both vertically and horizontally. You determine
// the order of the numbers.

////////////////////////////////////////////////////////////////////////////////

// Global variables.
int **rowreg, **colreg, **x;
int nclues, *cluerow, *cluecol, *clue;
int R, C, V;
double T, *P;

// These functions are found below.
void GetPuzzle ();
void Metropolis ();
int  Energy ();
void Proposal ();
void Probabilities ();
void Report ();
void AllocateMemory ();

// These functions are in common among all puzzle codes.
#include "MetropolisFunctions.h"


////////////////////////////////////////////////////////////////////////////////////
// Main program.
////////////////////////////////////////////////////////////////////////////////////
int main () {

   // Allocate array space.
   AllocateMemory ();

   // Get the puzzle and generate random initial configuration.
   GetPuzzle ();

   // Report the puzzle.
   Report ();

   // Solve the puzzle via Metropolis.
   Metropolis ();

   // Report the solution.
   Report ();

}

////////////////////////////////////////////////////////////////////////////////////
// Implement the Metropolis algorithm.
////////////////////////////////////////////////////////////////////////////////////
void Metropolis () {

   int E, deltaE, AcceptTransition, n=0;

   // Get the temperature parameter.
   T = GetDouble ("What is the temperature parameter (.4 seems ok)?... ");

   // Compute acceptance probabilities.
   if (T > 0) Probabilities ();

   // Time the calculations.
   Time ();
   
   // Calculate the energy for the initial random configuration.
   E = Energy ();

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

      // If T > 0, accept with small probability prob[deltaE] as previously calculated.
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
         x[R][C] = V;
      }

      // Periodically report the energy.
      if (n % 1000000 == 0) {
         if (n == 1000000) {
            printf ("        n    E\n");
            printf ("=========  ===\n");
         }   
         printf ("%9d  %3d\n", n, E);
      }

   }

   printf ("Solved after %.1f million steps of the Markov chain in %.3f seconds.\n",
            n/1000000.0, Time());

   return;

}

////////////////////////////////////////////////////////////////////////////////////
// Propose a random change to the configuration.
// Pick a random grid site and randomly change its allocated value.
////////////////////////////////////////////////////////////////////////////////////
void Proposal () {

   // R, C, and V are global variables.

   // Choose a grid site.
   R = RandomInteger (1,8);
   C = RandomInteger (1,8);

   // Record the current value for possible reinstatement.
   V = x[R][C];

   // Choose a new value.
   while (1) {
      x[R][C] = RandomInteger (1,8);
      if (x[R][C] != V) break;
   }

   return;

}

////////////////////////////////////////////////////////////////////////////////////
// Pre-calculate the probabilities of accepting a change in configuration as a
//    function of the change in energy. This improves run-time.
// If T = 0 this is not necessary.
////////////////////////////////////////////////////////////////////////////////////
void Probabilities () {

   int dE;

   for (dE = 1; dE <= 20; dE++) {
      P[dE] = exp (-dE/T);
   }

   return;

}

////////////////////////////////////////////////////////////////////////////////////
// Calculate the energy of the current configuration x[][].
////////////////////////////////////////////////////////////////////////////////////
int Energy () {

   int i, r, c, region, size, count[9], E;

   // Count clue violations.
   E = 0;
   for (i = 1; i <= nclues; i++) {
      r = cluerow[i];
      c = cluecol[i];
      E += 5 * (x[r][c] != clue[i]);  // The penalty of 5 is arbitrary.
   }

   // Row count discrepencies.
   for (r = 1; r <= 8; r++) {

      // Initialize region data.
      region = rowreg[r][1];                 // What region are we in?
      size = 0;                              // How big is the region?
      for (i = 1; i <= 8; i++) count[i] = 0; // How often does the number "i"
                                             //   appear in the region?
      // Move across row r.
      for (c = 1; c <= 8; c++) {
         size ++;
         i = x[r][c];
         count[i] ++;
         if (c == 8 || rowreg[r][c+1] != region) {
            // Update E to reflect current region's data.
            for (i = 1; i <= 8; i++) {
               if (i <= size) E += abs(count[i]-1);  // Should appear exactly once.
               else           E += 5*count[i];       // Should not appear at all.
            }
            // Re-initialize region data if not at row's end.
            if (c != 8) {
               region = rowreg[r][c+1];
               size = 0;
               for (i = 1; i <= 8; i++) count[i] = 0;
            }
         }
      }

   }

   // Column count discrepencies.
   for (c = 1; c <= 8; c++) {

      // Initialize region data.
      region = colreg[1][c];
      size = 0;
      for (i = 1; i <= 8; i++) count[i] = 0;

      // Move down column c.
      for (r = 1; r <= 8; r++) {
         size ++;
         i = x[r][c];
         count[i] ++;
         if (r == 8 || colreg[r+1][c] != region) {
            // Update E to reflect current region's data.
            for (i = 1; i <= 8; i++) {
               if (i <= size) E += abs(count[i]-1);  // Should appear exactly once.
               else           E += 5*count[i];       // Should not appear at all.
            }
            // Re-initialize region data if not at column's end.
            if (r != 8) {
               region = colreg[r+1][c];
               size = 0;
               for (i = 1; i <= 8; i++) count[i] = 0;
            }
         }
      }

   }

   return E;

}

////////////////////////////////////////////////////////////////////////////////////
// This function reads in the 8x8 puzzle.
////////////////////////////////////////////////////////////////////////////////////
void GetPuzzle () {

   int r, c, i;
   char input[100];
   FILE *fp=NULL;

   // Seed the RNG.
   MTUniform();

   while (fp == NULL) {

      printf ("Please input the name of the puzzle input file... ");
      fgets (input, 99, stdin);

      // Now terminate the file name with ".txt".
      for (i = 0; i <= 20; i++) {
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

   }

   // Read in row data. In the ASCII code, the value of '0' is 48 and the values
   //    of '1' through '9' are 49 through 57.
   for (r = 1; r <= 8; r++) {
      fgets (input, 99, fp);
      for (c = 1; c <= 8; c++) {
         rowreg[r][c] = input[c-1] - '0';
      }
   }

   // Blank line
   fgets (input, 99, fp);

   // Read in column data.
   for (r = 1; r <= 8; r++) {
      fgets (input, 99, fp);
      for (c = 1; c <= 8; c++) {
         colreg[r][c] = input[c-1] - '0';
      }
   }

   // Blank line
   fgets (input, 99, fp);

   // Read in clues.
   fgets (input, 99, fp);
   sscanf (input, "%d", &nclues);
   for (i = 1; i <= nclues; i++) {
      fgets (input, 99, fp);
      sscanf (input, "%d %d %d", cluerow+i, cluecol+i, clue+i);
   }

   // Close the data file.
   fclose (fp);

   // Generate a random initial configuration.
   for (r = 1; r <= 8; r++) {
      for (c = 1; c <= 8; c++) {
         x[r][c] = RandomInteger (1,8);
      }
   }

   return;

}

////////////////////////////////////////////////////////////////////////////////////
// Generate the output file OneUp.tex for viewing.
////////////////////////////////////////////////////////////////////////////////////
void Report () {

   int r, c;
   static int a=0;
   FILE *fp;

   // PUZZLE.
   fp = fopen ("PuzzleData.txt", "w");

   // Plot the vertical walls in each row.
   for (r = 1; r <= 8; r++) {
      for (c = 2; c <= 8; c++) if (rowreg[r][c-1] != rowreg[r][c]) {
         fprintf (fp, "\\plot %d %d %d %d /\n", c-1, 9-r, c-1, 8-r);
      }
   }

   // Plot the horizontal walls in each column.
   for (c = 1; c <= 8; c++) {
      for (r = 2; r <= 8; r++) if (colreg[r][c] != colreg[r-1][c]) {
         fprintf (fp, "\\plot %d %d %d %d /\n", c-1, 9-r, c, 9-r);
      }
   }

   // Report the puzzle clues.
   for (c = 1; c <= nclues; c++) {
      fprintf (fp, "\\put %d at %.1f %.1f\n", clue[c], cluecol[c]-0.5, 8.5-cluerow[c]);
   }

   fclose (fp);

   if (a == 0) {
      // Clear out solution data.
      fp = fopen ("Solution.txt", "w");
      fclose (fp);
      // Next call to Report() is for the solution.
      a = 1;
      return;
   }

   // SOLUTION.
   fp = fopen ("Solution.txt", "w");
   // Report the solution.
   for (r = 1; r <= 8; r++) {
      for (c = 1; c <= 8; c++) {
         fprintf (fp, "\\put {%d} at %.1f %.1f\n", x[r][c], c-0.5, 8.5-r);
      }
   }
   fclose (fp);

   printf ("View the puzzle and solution using Plain TeX with OneUp.tex.\n");
   Exit();

}

////////////////////////////////////////////////////////////////////////////////
// Allocate array space. All array entries are initialized to 0 by calloc().
////////////////////////////////////////////////////////////////////////////////
void AllocateMemory () {

   int i;

   P = (double *) calloc (21, sizeof (double));

   rowreg = (int **) calloc (9, sizeof (int *));
   for (i = 0; i < 9; i++) {
      rowreg[i] = (int *) calloc (9, sizeof (int));
   }

   colreg = (int **) calloc (9, sizeof (int *));
   for (i = 0; i < 9; i++) {
      colreg[i] = (int *) calloc (9, sizeof (int));
   }

   x = (int **) calloc (9, sizeof (int *));
   for (i = 0; i < 9; i++) {
      x[i] = (int *) calloc (9, sizeof (int));
   }

   cluerow = (int *) calloc (9, sizeof (int));
   cluecol = (int *) calloc (9, sizeof (int));
   clue    = (int *) calloc (50, sizeof (int));

   return;

}


