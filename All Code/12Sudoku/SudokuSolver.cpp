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
// Metropolis Sudoku solver as described in Section 12. The data file
// Puzzle1.txt is the example in the text. You may make additional puzzle
// data files using SudokuPuzzleMaker.cpp.
////////////////////////////////////////////////////////////////////////////////

// Global variables.
int  **x, **clue, i0, j0, x0;
double T, *prob;

// Algorithm functions found below.
void    Report ();                                              
void    GetPuzzle ();
void    Metropolis ();
int     Energy ();
int     Conflicts (int, int, int);
void    Proposal ();
void    Probabilities ();
void    AllocateMemory ();

// These functions are in common among all applications.
#include "MetropolisFunctions.h"

////////////////////////////////////////////////////////////////////////////////
// Main program. ///////////////////////////////////////////////////////////////
int main () {

   // Allocate array space.
   AllocateMemory ();

   // Read in the puzzle's clues and set up an initial configuration.
   GetPuzzle ();

   // Report the puzzle to the screen.
   Report ();

   // Solve the puzzle.
   Metropolis ();

   // Report the puzzle solution.
   Report ();

   // Pause so the window doesn't close.
   Pause ();

}

////////////////////////////////////////////////////////////////////////////////
// Solve the puzzle via Metropolis. ////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
void Metropolis () {

   int E, AcceptTransition, deltaE, n=0;

   // Get the temperature parameter.
   T = GetDouble ("What is the temperature (best = 0.39)?... ");

   // Unless T = 0, compute the acceptance probabilities once-and-for-all.
   if (T > 0) Probabilities ();

   // Time the computations.
   Time ();

   // Record initial energy.
   E = Energy ();

   // Continue as long as there is improvement to be had.
   while (E > 0) {

      // Count the number of Markov chain steps.
      n ++;

      // Periodically report the energy.
      if (n % 1000000 == 0) {
         if (n == 1000000) {
            printf ("\n");
            printf ("       n    E\n");
            printf ("========  ===\n");
         }   
         printf ("%8d  %3d\n", n, E);
      }

      // Determine the proposed transition. This function generates the cell (i0,j0)
      //    to be changed and the digit x0 that it is changed to. i0, j0, and x0 are
      //    global variables.
      Proposal ();

      // Compute the change in energy. Conflicts (c, r, d) is the number of row, column,
      //    and 3x3 block conflicts with cell (c,r) if that cell has digit d.
      deltaE = Conflicts (i0, j0, x0) - Conflicts (i0, j0, x[i0][j0]);

      // Determine if proposed transition is accepted.
      AcceptTransition = 0;

      if (deltaE <= 0) {
         AcceptTransition = 1;
      }

      // So far it's just the zero temperature dynamics.
      // If T > 0, accept a transition to higher energy with probability prob[deltaE].
      else if (T > 0) {
         if (MTUniform () <= prob[deltaE]) {
            AcceptTransition = 1;
         }
      }

      // Change the state and update the energy if the proposed transition is accepted.
      if (AcceptTransition) {
         x[i0][j0] = x0;
         E += deltaE;
      }

   }

   // Show n when E = 0.
   printf ("\nSolved after %d steps.  Computations took %.2f seconds", n, Time ());

}

////////////////////////////////////////////////////////////////////////////////
// Pick a cell (i0,j0) to randomly change the digit from x[i0][j0] to x0. //////
void Proposal () {

   // Determine the proposed transition...

   // Choose a random neighbor of the configuration.
   // - First select the cell (i0,j0).
   while (1) {
      i0 = RandomInteger (1, 9);
      j0 = RandomInteger (1, 9);
      // Accept if (i0,j0) is not a clue cell.
      if (!clue[i0][j0]) {
         break;
      }
   }

   // - Now select the new digit for that cell.
   while (1) {
      x0 = RandomInteger (1, 9);
      // Accept if x0 is different from the current value at (i0,j0).
      if (x0 != x[i0][j0]) {
         break;
      }
   }

   return;

}

////////////////////////////////////////////////////////////////////////////////
// Compute the acceptance probabilities at temperature T. //////////////////////
////////////////////////////////////////////////////////////////////////////////
void Probabilities () {

   int deltaE;

   for (deltaE = 1; deltaE <= 20; deltaE++) {
      prob[deltaE] = exp (-deltaE/T);
   }

   return;

}

////////////////////////////////////////////////////////////////////////////////
// Get the puzzle from a file and initialize certain information. //////////////
////////////////////////////////////////////////////////////////////////////////
void GetPuzzle () {

   int i, j, *d;
   char input[100];
   FILE *fp=NULL;


   printf ("I will solve any Sudoku puzzle for you!\n\n");

   // First get the file name.
   while (!fp) {

      printf ("Please input the name of the puzzle's data file... ");
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
      fp = fopen (input, "r");

   }

   // ...and read in the data.
   for (i = 1; i <= 9; i++) {

      // Read in row number i.
      d = x[i];
      fgets (input, 99, fp);
      sscanf (input, "%d %d %d %d %d %d %d %d %d",
                   d+1, d+2, d+3, d+4, d+5, d+6, d+7, d+8, d+9);

   }
   fclose (fp);

   // Identify the locations that are puzzle clues and randomize the
   //   other locations. This generates a random initial configuration.
   for (i = 1; i <= 9; i++) {
      for (j = 1; j <= 9; j++) {

         // Flag the clue locations.
         clue[i][j] = (x[i][j] != 0);

         // Randomize the non-clue locations.
         if (!clue[i][j]) {
            x[i][j] = RandomInteger (1, 9);
         }

      }
   }

   return;

}

////////////////////////////////////////////////////////////////////////////////
// Report the puzzle or the solution to the screen. ////////////////////////////
////////////////////////////////////////////////////////////////////////////////
void Report () {

   int i, j;
   static int a=0;

   printf ("\n");

   // If it's the solution...
   if (a) {
      printf ("\n");
      printf ("       Sudoku Solution\n");
   }

   // If it's the puzzle...
   else {
      printf ("        Sudoku Puzzle\n");
   }

   // Now print the clues/solution in a nice format.
   printf ("  =========================\n");
   for (i = 1; i <= 9; i++) {
      printf ("  |");
      for (j = 1; j <= 9; j++) {
         if ((!a && clue[i][j]) || a) {
            printf ("%2d", x[i][j]);
         } else {
            printf ("  ");
         }
         if (j == 3 || j == 6  || j == 9) {
            printf (" |");
         }
      }
      printf ("\n");
      if (i == 3 || i == 6) {
         printf ("  |=======|=======|=======|\n");
      }
   }
   printf ("  =========================\n");
   printf ("\n\n");

   // The second call to Report is for the solution.
   a = 1;

   return;

}

////////////////////////////////////////////////////////////////////////////////
// This function totals up the number conflicts in the initial configuration. //
////////////////////////////////////////////////////////////////////////////////
int Energy () {

   int E, i, j;
   E = 0;
   for (i = 1; i <= 9; i++) {
      for (j = 1; j <= 9; j++) {
         // How many conflicts are there with cell (i,j)?
         E += Conflicts (i, j, x[i][j]);
      }
   }

   // Conflicts are double counted: if x[a][b] conflicts with x[c][d] then
   //                                  x[c][d] conflicts with x[a][b].
   // Divide the energy by 2 to account for this.
   E /= 2;

   return (E);

}   

////////////////////////////////////////////////////////////////////////////////
// This function computes the number of conflicts with cell (r0,c0) in     /////
//   the configuration x when the digit x[r0][c0] is changed to d.         /////
// Conflicts involving a clue are given an energy of 5, other conflicts 1. /////
////////////////////////////////////////////////////////////////////////////////
int Conflicts (int r0, int c0, int d) {

   int c, r, k, left, lower, conflicts = 0;

   for (k = 1; k <= 9; k++) {

      // Row conflicts.
      if (k != c0) {
         conflicts += (x[r0][k] == d) * (clue[r0][c0] || clue[r0][k] ? 5 : 1);
      }

      // Column conflicts.
      if (k != r0) {
         conflicts += (x[k][c0] == d) * (clue[r0][c0] || clue[k][c0] ? 5 : 1);
      }
   }
   
   // Block conflicts: find the coordinates of the lower-left cell in (r0,c0)'s 3x3 block.
   lower = 1 + 3*((c0-1)/3); // Integer arithmetic!
   left  = 1 + 3*((r0-1)/3);

   // Now loop through the 9 cells in that block to find block conflicts that
   //    have not already been counted as row or column conflicts.
   for (r = left; r <= left+2; r++) {
      for (c = lower; c <= lower+2; c++) {

         // Insist that (r,c) is in both a different row and column then (r0,c0).
         if (r != r0  && c != c0) {
            conflicts += (x[r][c] == d) * (clue[r0][c0] || clue[r][c] ? 5 : 1);
         }

      }
   }

   return (conflicts);

}

///////////////////////////////////////////////////////////////////////////////////////////////
// Allocate array space. //////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////
void AllocateMemory () {

   int i;

   // Array space for transition-acceptance probabilities.
   prob = (double *) calloc (21, sizeof (double));

   // Allocate array space to hold the puzzle data.
   // x[i][j] is the digit in row i, column j of the puzzle.
   // Here 1 <= i <= 9 and 1 <= j <= 9.
   x = (int **) calloc (10, sizeof (int *));
   for (i = 1; i <= 9; i++) {
      x[i] = (int *) calloc (10, sizeof (int));
   }

   // clue[i][j] indicates whether (clue[i][j] = 1) or not (clue[i][j] = 0) the
   //  puzzle location row i, column j, is a clue.
   clue = (int **) calloc (10, sizeof (int *));
   for (i = 1; i <= 9; i++) {
      clue[i] = (int *) calloc (10, sizeof (int));
   }


}





