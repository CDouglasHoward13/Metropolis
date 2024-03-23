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
// Metropolis Sudoku puzzle maker, as described in Section 13. This program
// creates a random Sudoku puzzle and reports the clues to the file Clues.txt.
// Solve the puzzle using SudokuSolver.cpp.
////////////////////////////////////////////////////////////////////////////////

// Global variables -- metropolis1 stage.
int  **s, I0, J0, s0;
double *prob;

// Global variables -- metropolis2 stage.
int *r, *c, **Gtilde, **n, *m, **C, istar, Csize;

// Algorithm functions found below.
void    Metropolis1 ();
int     Energy ();
int     Conflicts (int, int, int);
void    Proposal ();
void    Probabilities ();
void    AllocateMemory ();
void    Randomize ();

// Reduction functions.
void    Metropolis2 ();
void    Report (int);
int     Consistent (int, int, int, int);
int     Flip (int, int);
int     Unique ();
void    GettingStarted ();
int     Advance ();
int     Retreat ();

// These functions are in common among all puzzle codes.
#include "MetropolisFunctions.h"

////////////////////////////////////////////////////////////////////////////////
// Main program. ///////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
int main () {

   // Allocate array space.
   AllocateMemory ();

   // Random initial configuration.
   Randomize ();

   // Solve the clueless puzzle. The solution is held in the s[i][j] array.
   Metropolis1 ();

   // Reduce puzzle solution to essential clues.
   Metropolis2 ();

   Pause ();


}

/***************************************************************************

             THESE FUNCTIONS ARE FOR SOLUTION GENERATION PHASE

***************************************************************************/


////////////////////////////////////////////////////////////////////////////////
// Solve the cluless puzzle via Metropolis.
////////////////////////////////////////////////////////////////////////////////
void Metropolis1 () {

   int E, AcceptTransition, deltaE;

   // Compute the acceptance probabilities once-and-for-all.
   Probabilities ();

   // Record initial energy.
   E = Energy ();

   // Continue as long as there is improvement to be had.
   while (E > 0) {

      // Determine the proposed transition. This function generates the cell (I0,J0)
      //    to be changed and the digit s0 that it is changed to. I0, J0, and s0 are
      //    global variables.
      Proposal ();

      // Compute the change in energy. Conflicts (r, c, z) is the number of row, column,
      //    and 3x3 block conflicts with cell (r,c) if that cell has digit z.
      deltaE = Conflicts (I0, J0, s0) - Conflicts (I0, J0, s[I0][J0]);

      // Determine if proposed transition is accepted.
      AcceptTransition = 0;

      if (deltaE <= 0) {
         AcceptTransition = 1;
      }

      // So far it's just the zero temperature dynamics.
      // If deltaE > 0, accept a transition to higher energy with probability prob[deltaE].
      else if (MTUniform () <= prob[deltaE]) {
         AcceptTransition = 1;
      }

      // Change the state and update the energy if the proposed transition is accepted.
      if (AcceptTransition) {
         s[I0][J0] = s0;
         E += deltaE;
      }

   }

}

/////////////////////////////////////////////////////////////////////////////////
// Generate a random initial configuration. /////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
void Randomize () {

   int x, y;

   // Seed the RNG.
   printf ("I will generate a random Sudoku puzzle for you!\n");
   MTUniform ();

   for (x = 1; x <= 9; x++) {
      for (y = 1; y <= 9; y++) {
         s[x][y] = RandomInteger (1, 9);
      }
   }


}

////////////////////////////////////////////////////////////////////////////////
// Pick a cell (I0,J0) to randomly change the digit from s[I0][J0] to s0.
////////////////////////////////////////////////////////////////////////////////
void Proposal () {

   // Determine the proposed transition...

   // Choose a random neighbor of the configuration.
   // - First select the cell (I0,J0).
   I0 = RandomInteger (1, 9);
   J0 = RandomInteger (1, 9);

   // - Now select the new digit for that cell.
   s0 = s[I0][J0];
   while (s0 == s[I0][J0]) {
      s0 = RandomInteger (1, 9);
   }

   return;

}

////////////////////////////////////////////////////////////////////////////////
// Compute the acceptance probabilities at temperature T.
////////////////////////////////////////////////////////////////////////////////
void Probabilities () {

   int deltaE;
   double T = 0.2;

   for (deltaE = 1; deltaE <= 20; deltaE++) {
      prob[deltaE] = exp (-deltaE/T);
   }

   return;

}

////////////////////////////////////////////////////////////////////////////////
// This function totals up the number conflicts in the initial configuration.
////////////////////////////////////////////////////////////////////////////////
int Energy () {

   int E, x, y;

   E = 0;
   for (x = 1; x <= 9; x++) {
      for (y = 1; y <= 9; y++) {
         // How many conflicts are there with cell (x,y)?
         E += Conflicts (x, y, s[x][y]);
      }
   }

   // Conflicts are double counted: if s[a][b] conflicts with s[c][d] then
   //                                  s[c][d] conflicts with s[a][b].
   // Divide the energy by 2 to account for this.
   E /= 2;

   return (E);

}   

////////////////////////////////////////////////////////////////////////////////
// This function computes the number of conflicts with cell (x0,y0) in
//   the configuration x when the digit s[x0][y0] is set to d.
// Conflicts involving a clue are given an energy of 5, other conflicts 1.
////////////////////////////////////////////////////////////////////////////////
int Conflicts (int x0, int y0, int d) {

   int x, y, k, left, lower, conflicts = 0;

   for (k = 1; k <= 9; k++) {

      // Row conflicts.
      if (k != y0) {
         conflicts += (s[x0][k] == d);
      }

      // Column conflicts.
      if (k != x0) {
         conflicts += (s[k][y0] == d);
      }
   }
   
   // Block conflicts: find the coordinates of the lower-left cell in (c0,r0)'s 3x3 block.
   left  = 1 + 3*((y0-1)/3); // Integer arithmetic!
   lower = 1 + 3*((x0-1)/3);

   // Now loop through the 9 cells in that block to find block conflicts that
   //    have not already been counted as row or column conflicts.
   for (x = lower; x <= lower+2; x++) {
      for (y = left; y <= left+2; y++) {

         // Insist that (x,y) is in both a different row and column then (r0,c0).
         if (x != x0  && y != y0) {
            conflicts += (s[x][y] == d);
         }

      }
   }

   return (conflicts);

}

/***************************************************************************

            THESE FUNCTIONS ARE FOR THE REDUCTION-TO-CLUES PHASE

***************************************************************************/


////////////////////////////////////////////////////////////////////////////////
// Remove clues from the configuration until a specified difficulty is reached
//    also via Metropolis. The resulting clue sites have symmetry about
//    the center site (5,5).
////////////////////////////////////////////////////////////////////////////////
void Metropolis2 () {

   int x, y, d, E, DeltaClues;
   double p[3];

   // Acceptance probabilities.
   p[1] = 0.1581; p[2] = .0250;

   // In the initial configuration every site is a clue.
   for (x = 1; x <= 9; x++) {
      for (y = 1; y <= 9; y++) {
         C[x][y] = 1;
      }
   }

   // Report the puzzle solution.
   Report (0);

   // Get the desired level of puzzle difficulty.
   d = -1;
   while (d < 0 || d > 10) {
      d =  GetInteger ("What difficulty do you want (integer from 0 = Easy to 10 = Evil)?... ");
   }
   if (d == 10) {
      printf ("Be patient, this could take a few seconds...\n");
   }   

   // The number of clues for this level of difficulty is...
   d = 35 - d;

   // The energy of initial configuration is 81.
   E = 81;

   // Work the number of clues down to d.
   while (E > d) {

      /* Pick a random site from rows 1 - 4, and row 5 through column 5.  There
         are 41 such sites:
      =========================
      | * * * | * * * | * * * |
      | * * * | * * * | * * * |
      | * * * | * * * | * * * |
      |=======|=======|=======|
      | * * * | * * * | * * * |
      | * * * | * *   |       |
      |       |       |       |
      |=======|=======|=======|
      |       |       |       |
      |       |       |       |
      |       |       |       |
      =========================    */
      x = RandomInteger (1,5);
      y = RandomInteger (1,(x == 5 ? 5 : 9));

      // Flip it and its reflection about (5,5).  Compute the change in the number of clues.
      DeltaClues = Flip (x, y);

      // Accept if clues are reduced and solution is still unique. Otherwise, flip back.
      if (DeltaClues < 0) {
         if (Unique ()) {
            E += DeltaClues;
         } else {
            Flip (x, y);
         }
      }

      // If clues are increased accept with small probability. Otherwise, flip back.
      else {
         if (MTUniform() <= p[DeltaClues]) {
            E += DeltaClues;
         } else {
            Flip (x, y);
         }
      }

   }

   // Report the essential clues.
   Report (1);

   return;

}

////////////////////////////////////////////////////////////////////////////////
// Implement the flips and report the number of clue changes.
////////////////////////////////////////////////////////////////////////////////
int Flip (int x0, int y0) {

   int x1, y1, both55, DeltaClues;

   // Flip location (x0,y0).
   C[x0][y0] = 1 - C[x0][y0];

   // Sites (x0,y0) and (10-x0,10-y0) may both be (5,5).
   both55 = (x0 == 5 && y0 == 5);

   // If not, flip the reflection about (5,5).
   if (!both55) {
      // It's reflection about (5,5) is...
      x1 = 10 - x0;
      y1 = 10 - y0;
      C[x1][y1] = 1 - C[x1][y1];
   }

   // Change in the number of clues.
   DeltaClues = (C[x0][y0] == 1 ? +1 : -1) * (both55 ? 1 : 2);

   return DeltaClues;

}

////////////////////////////////////////////////////////////////////////////////
// Is the solution to P(C) unique?   1 = yes, 0 = no.
////////////////////////////////////////////////////////////////////////////////
int Unique () {

   int status;

   // Set up the clue data.
   GettingStarted ();

   // Status values:
   // - 0 means Advance/Retreat is done. istar is either |C| or 81.
   // - 1 means we're in the Advance state.
   // - 2 means we're in Retreat.

   // Find the first solution. Begin by advancing.
   status = 1;
   while (status != 0) {

      if (status == 1) {
         status = Advance ();
      } else {
         status = Retreat ();
      }

   }

   // This should never happen:
   if (istar != 81) {
      printf ("Problem! No solution was found.\n");
      Exit ();
   }

   // Look for a second solution. Begin with a retreat.
   status = 2;
   while (status != 0) {

      if (status == 1) {
         status = Advance ();
      } else {
         status = Retreat ();
      }

   }

   // If istar is 81, a second solution was found.
   return (istar == 81 ? 0 : 1);

}

////////////////////////////////////////////////////////////////////////////////
// Setup for the advance/retreat algorithm.
////////////////////////////////////////////////////////////////////////////////
void GettingStarted () {

   int i, x, y;

   // Getting started.
   i = 0;
   for (x = 1; x <= 9; x++) {
      for (y = 1; y <= 9; y++) {

         if (C[x][y]) {
            i++;
            r[i] = x;
            c[i] = y;
            m[i] = 1;
            n[i][1] = s[x][y];
            Gtilde[x][y] = 1; // 1 indicates that (x,y) is in Gtilde.
         }

         else {
            Gtilde[x][y] = 0; // 0 indicates that (x,y) is not in Gtilde.
         }

      }
   }
   Csize = istar = i;   // Csize is |C| in the text.

   return;

}

////////////////////////////////////////////////////////////////////////////////
// Continue advancing toward a solution.
////////////////////////////////////////////////////////////////////////////////
int Advance () {

   int j, mprime, rprime, cprime, x, y, m0, v;

   // If istar is 81, we have found a solution and are done.
   if (istar == 81) return 0;

   // Compute m', r', and c'.
   mprime = 10; // m' must be less than this.

   // Loop through grid sites not in Gtilde.
   for (x = 1; x <= 9; x++) {
      for (y = 1; y <= 9; y++) if (!Gtilde[x][y]) {

         // Grid site (x,y) is not in Gtilde.

         // How many digit values v at (x,y) are consistent with the trial
         //    solution through i*?  Put that number in m0.
         m0 = 0;
         for (v = 1; v <= 9; v++) {
            if (Consistent (x, y, v, istar)) {
               m0 ++;
            }
         }

         // Is that the smallest number so far? If so, update m', r', and c'.
         if (m0 < mprime) {
            mprime = m0;
            rprime = x;
            cprime = y;
         }

        // Case 1 (m' is 0), something must be wrong here. Retreat.
        if (mprime == 0) {
           return 2;
        }

      }
   }

   // Case 2, record the digit values at (r',c') that are
   //    consistent with the trial solution so far (through i*-1,
   //    after augmenting i*).
   istar ++;
   j = 0;
   for (v = 1; v <= 9; v++) {
      if (Consistent (rprime, cprime, v, istar-1)) {
         j ++;
         n[istar][j] = v;
      }
   }
   // Update r, c, Gtilde, and m data at i*.
   r[istar] = rprime;
   c[istar] = cprime;
   m[istar] = mprime;
   Gtilde[rprime][cprime] = 1; // Put (r',c') in Gtilde.

   // Continue advancing.
   return 1;

}

////////////////////////////////////////////////////////////////////////////////
// Retreat from current proposed solution.
////////////////////////////////////////////////////////////////////////////////
int Retreat () {

   int x, y;

   // Stop if we've retreated back to the clues.
   if (istar == Csize) return 0;

   // Case 1 (m[i*] > 1), try next consistent value at (r[i*],c[i*]) in the trial
   //    solution through i*. Switch to advancing.
   if (m[istar] > 1) {
      m[istar] --;
      return 1;
   }

   // Case 2 (m[i*] is 1), remove i* data from the record.
   x = r[istar];
   y = c[istar];
   Gtilde[x][y] = 0; // Remove (r[i*],c[i*]) from Gtilde.
   istar --; // Reduce i* by one.

   // Continue retreating.
   return 2;

}

///////////////////////////////////////////////////////////////////////////////
// Is the digit v at (x, y) consistent with the trial solution through k?
// Yes = 1, no = 0.
///////////////////////////////////////////////////////////////////////////////
int Consistent (int x, int y, int v, int k) {

   int i, b0, b;

   // Compute (x,y)'s 3x3 block number. Integer arithmetic.
   b0 = (x-1)/3 + 3*((y-1)/3);

   for (i = 1; i <= k; i++) {

      // n[i][m[i]] is the digit value at (r[i],c[i]) in the trial solution
      //    through k.
      if (v == n[i][m[i]]) {
         if (r[i] == x) return 0; // Same row.
         if (c[i] == y) return 0; // Same column.
         b = (r[i]-1)/3 + 3*((c[i]-1)/3);
         if (b == b0)   return 0; // Same 3x3 block.
      }
   }

   // If no conflicts found...
   return 1;

}

////////////////////////////////////////////////////////////////////////////////
// Report the puzzle or the solution to the screen.
// Create the Clues.txt data file for SudokuSolver.cpp input.
////////////////////////////////////////////////////////////////////////////////
void Report (int a) {

   int i, j;
   FILE *fp;

   printf ("\n");

   // If it's the puzzle...
   if (a == 0) {
      printf ("\n");
      printf ("       Sudoku Solution\n");
   }

   // If it's the solution...
   else {
      printf ("        Sudoku Puzzle\n");
   }

   if (a) fp = fopen ("Clues.txt", "w");

   // Now print the clues/solution in a nice format.
   printf ("  =========================\n");
   for (i = 1; i <= 9; i++) {
      printf ("  |");
      for (j = 1; j <= 9; j++) {
         if (a) {
            if (C[i][j]) fprintf (fp, "%2d", s[i][j]);
            else         fprintf (fp, "%2d", 0);
         }   
         if (C[i][j]) {
            printf ("%2d", s[i][j]);
         } else {
            printf ("  ");
         }
         if (j == 3 || j == 6  || j == 9) {
            printf (" |");
         }
      }
      printf ("\n");
      if (a) fprintf (fp, "\n");
      if (i == 3 || i == 6) {
         printf ("  |=======|=======|=======|\n");
      }
   }
   printf ("  =========================\n");
   printf ("\n\n");

   // Now make the file for Puzzle.tex.
   if (a) {
      fclose (fp);
      fp = fopen ("PuzzleClues.txt", "w");
      for (i = 1; i <= 9; i++) {
         for (j = 1; j <= 9; j++) if (C[i][j]) {
            fprintf (fp, "\\put {%d} at %3.1f %3.1f\n", s[i][j], j-0.5, 9.5-i);
         }
      }
      fclose (fp);
      printf ("For a better rendition use Puzzle.tex with Plain TeX.\n");
   }

   return;

}

////////////////////////////////////////////////////////////////////////////////
// Allocate array space.
////////////////////////////////////////////////////////////////////////////////
void AllocateMemory () {

   int i;

   // Metropolis 2 phase, reducing to a set of clues:
   r = (int *) calloc (82, sizeof (int));
   c = (int *) calloc (82, sizeof (int));
   m = (int *) calloc (82, sizeof (int));
   n = (int **) calloc (82, sizeof (int *));

   for (i = 1; i < 82; i++) {
      n[i] = (int *) calloc (11, sizeof (int));
   }

   Gtilde = (int **) calloc (10, sizeof (int *));
   for (i = 0; i < 10; i++) {
      Gtilde[i] = (int *) calloc (10, sizeof (int));
   }
   
   C = (int **) calloc (10, sizeof (int *));
   for (i = 0; i < 10; i++) {
      C[i] = (int *) calloc (100, sizeof (int));
   }

   // Metropolis 1 phase, generating the solution:
   prob = (double *) calloc (21, sizeof (double));

   s = (int **) calloc (10, sizeof (int *));
   for (i = 1; i <= 9; i++) {
      s[i] = (int *) calloc (10, sizeof (int));
   }

}








