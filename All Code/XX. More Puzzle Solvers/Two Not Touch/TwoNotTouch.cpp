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
// Metropolis "Two Not Touch" puzzle solver. Two Not Touch appears in the Arts
// section of the New York Times under a www.krazydad.com copyright.
////////////////////////////////////////////////////////////////////////////////

// Global variables.
int **col, **reg;
double T, *P;
int R, S, C;

// Functions found below.
void GetPuzzle ();
int  Energy ();
int  NearestNeighbor (int, int, int, int);
void Proposal ();
void Metropolis ();
void Report ();
void AllocateMemory ();

#include "MetropolisFunctions.h"

////////////////////////////////////////////////////////////////////////////////
// Main program. ///////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
int main () {

   // Allocate array space.
   AllocateMemory ();

   // Get the puzzle, compute acceptance probabilities,
   //    generate random initial configuration.
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
         col[R][S] = C;
      }

      // Periodically report progress.
      if (n % 1000000 == 0) printf ("%10d %2d\n", n, E);

   }

   printf ("Solved after %.1f million steps of the Markov chain in %.3f seconds.\n",
            n/1000000.0, Time ());

}

////////////////////////////////////////////////////////////////////////////////
// Propose a random change to the configuration. ///////////////////////////////
////////////////////////////////////////////////////////////////////////////////
void Proposal () {

   int c;   // R, S, and C are global variables.

   // Choose a star to move.
   R = RandomInteger (1,10); // Choose a random row.
   S = RandomInteger (1,2);  // Choose a random star in that row.

   // Record the current column location of that star for possible reinstatement.
   C = col[R][S];

   // Column location of the other star (3-S sends 1 to 2, and 2 to 1).
   c = col[R][3-S];

   // Choose a random new column for star number S different from both C and c.
   while (1) {
      col[R][S] = RandomInteger (1,10);
      if (col[R][S] != C && col[R][S] != c) break;
   }

   return;

}

////////////////////////////////////////////////////////////////////////////////
// Calculate the energy of the current configuration col[][]. //////////////////
////////////////////////////////////////////////////////////////////////////////
int Energy () {

   int r, s, c, E, region, column, regcount[11], colcount[11];

   // Initialize region and column star counts.
   for (c = 1; c <= 10; c++) {
      regcount[c] = colcount[c] = 0;
   }

   // Count the number of stars in each column and region (the rows always have 2).
   for (r = 1; r <= 10; r++) {
      for (s = 1; s <= 2; s++) {
         column = col[r][s];
         region = reg[r][column];
         colcount[column] ++;
         regcount[region] ++;
      }
   }

   // Augment energy for regions and columns differing from a star count of 2.
   E = 0;
   for (c = 1; c <= 10; c++) {
      E += abs(regcount[c]-2) + abs(colcount[c]-2);
   }

   // Now augment energy for "touching" stars...
   // First check for touching in the same row.
   for (r = 1; r <= 10; r++) {
      E += (abs(col[r][1] - col[r][2]) == 1);
   }

   // Then check for touching in adjacent rows r and r+1.
   for (r = 1; r <= 9; r++) {
      E += (abs(col[r][1] - col[r+1][1]) <= 1);
      E += (abs(col[r][1] - col[r+1][2]) <= 1);
      E += (abs(col[r][2] - col[r+1][1]) <= 1);
      E += (abs(col[r][2] - col[r+1][2]) <= 1);
   }

   return E;

}

////////////////////////////////////////////////////////////////////////////////
// Return 1 if squares (a,b) and (c,d) are nearest neighbors; 0 if not.
////////////////////////////////////////////////////////////////////////////////
int NearestNeighbor (int a, int b, int c, int d) {

   return (abs(a-c) + abs(b-d) == 1);
}

////////////////////////////////////////////////////////////////////////////////
// This function reads in the 10x10 puzzle's regions.
// There are 10 regions each with at most 25 squares.
// The data is held in reg[r][c], which is the region number (1-10) of the square
//    at row r, column c.
// Input file format looks like:
/*
AAABBBBBBB
ACAAABDDDD
ACCCABBBBD
AEECAFFDDD
EEGCAAFDDD
EEGGGAFFDD
GGGHHAAFDD
GGGHAAIIJJ
GGHHHIIJJJ
GHHHHHHHHJ
*/
// Each letter corresponds to a particular region.  Letters must be A,B,C,D,E,F,G,H,I,J.
// An initial random configuration is generated.
// Proposal acceptance probabilities are pre-calculated.
////////////////////////////////////////////////////////////////////////////////
void GetPuzzle () {

   int i, r, c, dE;
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

   // Read in the data for each row.
   for (r = 1; r <= 10; r++) {
      fgets (input, 99, fp);
      // In the ASCII code, the value of @ is 64, and A's value is 65. So, for example,
      //    the numerical value of 'J' - '@' is 10.
      for (c = 1; c <= 10; c++) {
         reg[r][c] = input[c-1] - '@';
      }
   }
   fclose (fp);

   // Now generate a random initial configuration.
   // The configuration is held in col[r][s], where 1 <= r <= 10 and s = 1 or 2 (each row
   //    has two stars).
   // col[r][c] is the column location of star number c in row r.
   for (r = 1; r <= 10; r++) {
      col[r][1] = RandomInteger (1,10);  // Random column for star number 1 in row r.
      while (1) {              // Random different column for star number 2 in row r.
         col[r][2] = RandomInteger (1,10);
         if (col[r][2] != col[r][1]) break;
      }
      // Now col[r][1] and col[r][2] are two different randomly selected column numbers.
   }

   // Now get the temperature and calculate acceptance probabilities.
   T = GetDouble ("What is the temperature parameter (best = 0.448)?... ");
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

   int r, c, s, r0, c0;
   double x, y;
   FILE *fp;
   static int a=0;

   // PUZZLE.
   fp = fopen ("PuzzleRegions.txt", "w");
   // Loop through all the squares.
   for (r = 1; r <= 10; r++) for (c = 1; c <= 10; c++) {
      // Compute the coordinates of square (r,c)'s upper-left corner.
      x = c - 1;
      y = 11 - r;
      // Loop through all squares.
      for (r0 = 1; r0 <= 10; r0++) for (c0 = 1; c0 <= 10; c0++)
         if ( NearestNeighbor (r,c,r0,c0) && (reg[r][c] != reg[r0][c0]) ) {
            // Northern, southern, eastern or western neighbor?
            if (r0 == r-1) {
               fprintf (fp, "\\plot %f %f  %f %f /\n", x, y, x+1, y);
            }
            else if (r0 == r+1) {
               fprintf (fp, "\\plot %f %f  %f %f /\n", x, y-1, x+1, y-1);
            }
            else if (c0 == c+1) {
               fprintf (fp, "\\plot %f %f  %f %f /\n", x+1, y, x+1, y-1);
            }
            else if (c0 == c-1) {
               fprintf (fp, "\\plot %f %f  %f %f /\n", x, y, x, y-1);
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
   for (r = 1; r <= 10; r++) {
      for (s = 1; s <= 2; s++) {
         c = col[r][s];
         // Compute the coordinates of square (r,c)'s upper-left corner.
         x = c - 1;
         y = 11 - r;
         // Report the star.
         fprintf (fp, "\\put {$\\star$} at %f %f\n", x+0.5, y-0.5);
      }
   }
   fclose (fp);

   printf ("View the puzzle and solution using Plain TeX with TNT.tex.\n");
   Exit();

}

////////////////////////////////////////////////////////////////////////////////
// Allocate array space. ///////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
void AllocateMemory () {

   int r;

   col = (int **) calloc (11, sizeof (int *));
   reg = (int **) calloc (11, sizeof (int *));

   for (r = 1; r <= 10; r++) {
      col[r] = (int *) calloc ( 3, sizeof (int));
      reg[r] = (int *) calloc (11, sizeof (int));
   }

   P = (double *) calloc (25, sizeof (double));

}


