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
// Metropolis KenKen solver. See Section 15. The data file PuzzleA is
// for the puzzle shown in the text.  Look for 7x7 KenKen puzzles
// in the New York Time Magazine.
////////////////////////////////////////////////////////////////////////////////

// Global variables.
int **puzzle, *x1, I0, I1, *ans, **N;
char *op;
double T, *prob;

void GetPuzzle ();
void Initialize ();
void Neighbors ();
int  Region (int);
int  Energy ();
int  Operation (int);
void Probabilities ();
void Proposal ();
void ChangeBack ();
void Metropolis ();
void Report ();
void AllocateMemory ();

// These functions are in common among all metropolis applications.
#include "MetropolisFunctions.h"

////////////////////////////////////////////////////////////////////////////////////
// Main program.
////////////////////////////////////////////////////////////////////////////////////
int main () {

   // Allocate array space.
   AllocateMemory ();

   // Get the puzzle and initialize the configuration.
   GetPuzzle ();

   // Solve the puzzle via Metropolis.
   Metropolis ();

   // Set up the output files.
   Report ();

   // Pause so the window doesn't close.
   Pause ();

}

////////////////////////////////////////////////////////////////////////////////////
// Implement the Metropolis algorithm.
////////////////////////////////////////////////////////////////////////////////////
void Metropolis () {

   int E, deltaE, AcceptTransition;
   int n = 0;

   // Get the temperature parameter.
   T = GetDouble ("\nWhat is the temperature parameter (.5 seems ok)?... ");

   // Compute acceptance probabilities.
   if (T > 0) Probabilities ();

   // Time the calculations.
   Time ();
   
   // Calculate the energy for the initial configuration.
   E = Energy ();

   // Keep going until the solution is found.
   while (E > 0) {

      // Count Markov chain steps.
      n ++;

      // Propose a random change to the current configuration.
      Proposal ();

      // Compute the resulting change in energy.
      deltaE = Energy () - E;

      // Start with the zero temperature dynamics.
      AcceptTransition = 0;

      // If deltaE <= 0, accept the transition.
      if (deltaE <= 0) {
         AcceptTransition = 1;
      }

      // If T > 0, accept with small probability prob[deltaE] as previously calculated.
      else if (T > 0) {
         if (MTUniform() <= prob[deltaE]) {
            AcceptTransition = 1;
         }
      }

      // If accepted, keep the change and update the energy.
      if (AcceptTransition) {
         E += deltaE;
      }

      // Otherwise change back to previous configuration.
      else {
         ChangeBack ();
      }

      // Periodically report the energy.
      if (n % 10000000 == 0) {
         if (n == 10000000) {
            printf ("\n");
            printf ("        n    E\n");
            printf ("=========  ===\n");
         }   
         printf ("%9d  %3d\n", n, E);
      }

   }

   printf ("\nSolved after %.1f million steps of the Markov chain in %.1f seconds.\n\n",
   n/1000000.0, Time());

   return;

}
 
////////////////////////////////////////////////////////////////////////////////////
// Reverse the change made in Proposal().
////////////////////////////////////////////////////////////////////////////////////
void ChangeBack () {

   int temp;

   temp = x1[I0];
   x1[I0] = x1[I1];
   x1[I1] = temp;

}

////////////////////////////////////////////////////////////////////////////////////
// Propose a random change to the configuration.
// Pick two grid sites with different digits and swap their digits.
////////////////////////////////////////////////////////////////////////////////////
void Proposal () {

   int temp;

   I0 = RandomInteger (1, 49);

   I1 = I0;
   while (x1[I1] == x1[I0]) {
      I1 = RandomInteger (1, 49);
   }

   // Swap digits at sites I0 and I1.
   temp = x1[I0];
   x1[I0] = x1[I1];
   x1[I1] = temp;

}

////////////////////////////////////////////////////////////////////////////////////
// Pre-calculate the probabilities of accepting a change in configuration as a
//    function of the change in energy. This improves run-time.
// If T = 0 this is not necessary.
////////////////////////////////////////////////////////////////////////////////////
void Probabilities () {

   int deltaE;

   for (deltaE = 1; deltaE <= 10; deltaE++) {
      prob[deltaE] = exp (-deltaE/T);
   }

   return;

}

////////////////////////////////////////////////////////////////////////////////////
// Calculate the energy of the current configuration x1[*].
////////////////////////////////////////////////////////////////////////////////////
int Energy () {

   int count[8];
   int E, c, k, r, i;

  // Initialize.
   E = 0;

   // Column deviations from {1,2,3,4,5,6,7}.
   for (c = 1; c <= 7; c++) {
      for (i = 0; i <= 7; i++) {
         count[i] = 0;
      }
      k = c;
      count[x1[k]] ++;
      for (r = 2; r <= 7; r++) {
         k += 7;
         count[x1[k]] ++;
      }

      for (i = 1; i <= 7; i++) {
         E += abs(count[i] - 1);
      }
   }   

   // Row deviations from {1,2,3,4,5,6,7}.
   k = 0;
   for (r = 1; r <= 7; r++) {
      for (i = 0; i <= 7; i++) {
         count[i] = 0;
      }
      for (c = 1; c <= 7; c++) {
         k ++;
         count[x1[k]] ++;
      }

      for (i = 1; i <= 7; i++) {
         E += abs(count[i] - 1);
      }
   }

   // Regions where the target number is not hit.
   for (k = 1; k <= puzzle[0][0]; k++) {
      E += (Operation (k) != ans[k]);
   }

   return E;

}

////////////////////////////////////////////////////////////////////////////////////
// Perform the indicated operation on the digits in region k.
////////////////////////////////////////////////////////////////////////////////////
int Operation (int k) {

   int *p, i, answer;

   p = puzzle[k];

   if (op[k] == '+') {
      answer = 0;
      for (i = 1; i <= p[0]; i++) {
         answer += x1[p[i]];
      }
   }

   else if (op[k] == '-') {
      answer = x1[p[1]];
      for (i = 2; i <= p[0]; i++) {  // p[0] should be 2 here.
         answer -= x1[p[i]];
         if (answer < 0) answer *= -1;
      }
   }

   else if (op[k] == '*') {
      answer = 1;
      for (i = 1; i <= p[0]; i++) {
         answer *= x1[p[i]];
      }
   }

   else if (op[k] == '/') { // division, p[0] should be 2.
      answer = -1;
      if (x1[p[1]] % x1[p[2]] == 0) {
         answer = x1[p[1]] / x1[p[2]];
      } else if (x1[p[2]] % x1[p[1]] == 0) {
         answer = x1[p[2]] / x1[p[1]];
      } 
   }

   else {
      printf ("Cannot find operation.\n");
      Exit ();
   }   

   return answer;

}

////////////////////////////////////////////////////////////////////////////////////
// Construct initial random configuration.
////////////////////////////////////////////////////////////////////////////////////
void Initialize () {

   int i, k;

   // Allocate the digits {1,2,3,4,5,6,7} in equal number (7 each) to the grid sites.
   x1[0] = -1;
   for (k = 0; k < 49; k++) {
      i = 0;
      // Pick a random unoccupied grid site.
      while (x1[i] != 0) {
         i = RandomInteger (1, 49);
      }
      // Assign it a digit so that all seven digits are present in equal number.
      x1[i] = (k%7) + 1;
   }
}

////////////////////////////////////////////////////////////////////////////////////
// This function reads in the 7x7 puzzle's regions and numerical clues.
////////////////////////////////////////////////////////////////////////////////////
void GetPuzzle () {

   int i, j, k, n, r, m, *p;
   char input[100], *input0;
   FILE *fp=NULL;

   printf ("I'll solve any 7x7 KenKen puzzle for you.\n\n");

   while (fp == NULL) {

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

   }   

   // Read in the regions.
   // n holds the number of regions.
   n = k = 0;
   for (i = 1; i <= 7; i++) {
      fgets (input, 99, fp);
      for (j = 0; j < 7; j++) {
         k++; // Current grid site number.
         r = input[j] - 'A' + 1;
         if (r > n) n = r;
         p = puzzle[r];
         p[0] ++;
         m = p[0];
         p[m] = k;
      }
   }

   // Put the number of regions into puzzle[0][0].
   puzzle[0][0] = n;

   // Now get the numerical clues -- one for each of the n regions.
   for (i = 1; i <= n; i++) {
      fgets (input, 99, fp);
      op[i] = input[2];
      input0 = input+3;
      sscanf (input0, "%d", ans+i);
   }

   fclose (fp);

   // Construct a random initial configuration.
   Initialize ();

   return;

}

////////////////////////////////////////////////////////////////////////////////////
// List the neighboring squares of each square.
////////////////////////////////////////////////////////////////////////////////////
void Neighbors () {

   int i, n;

   for (i = 1; i <= 49; i++) {

      // Northern neighbor.
      if (i > 7) {
         // Northern.
         N[i][0] ++;
         n = N[i][0];
         N[i][n] = i - 7;
      }

      // Southern neighbor.
      if (i < 49 - 7 + 1) {
         N[i][0] ++;
         n = N[i][0];
         N[i][n] = i + 7;
      }

      // Eastern neighbor.
      if (i % 7 != 0) {
         N[i][0] ++;
         n = N[i][0];
         N[i][n] = i + 1;
      }

      // Western neighbor.
      if (i % 7 != 1) {
         N[i][0] ++;
         n = N[i][0];
         N[i][n] = i - 1;
      }

   }

   return;

}

////////////////////////////////////////////////////////////////////////////////////
// Which region is square k in?
////////////////////////////////////////////////////////////////////////////////////
int Region (int k) {

   int n, i, j, *p;

   n = puzzle[0][0];
   for (i = 1; i <= n; i++) {
      p = puzzle[i];
      for (j = 1; j <= p[0]; j++) {
         if (p[j] == k) return i;
      }
   }

   printf ("Could not find which region square %d is in.\n", k);
   Exit();

   return 0;

}

////////////////////////////////////////////////////////////////////////////////////
// Generate the output files for viewing.
////////////////////////////////////////////////////////////////////////////////////
void Report () {

   int k, n, m;
   double x, y;
   FILE *fp, *fps;

   Neighbors ();

   fp = fopen ("PuzzleRegions.txt", "w");
   fps = fopen ("PuzzleSolution.txt", "w");

   // Loop through all the squares.
   for (k = 1; k <= 49; k++) {

      // Compute the coordinates of square k's upper-left corner.
      x = (k-1) % 7;
      y = 7 - (k-1) / 7;

      fprintf (fps, "\\put {\\bf %d} at %f %f\n", x1[k], x+0.5, y-0.5);

      // Loop through all neighboring squares of k.
      for (n = 1; n <= N[k][0]; n++) {

         // Square m is a neighbor of k.
         m = N[k][n];

         // If k and m are in different regions draw a boundary segment.
         if (Region(k) != Region(m)) {

            // Is m k's northern, southern, eastern or western neighbor.
            if (m == k - 7) {
               fprintf (fp, "\\plot %f %f  %f %f /\n", x, y, x+1, y);
            }
            else if (m == k + 7) {
               fprintf (fp, "\\plot %f %f  %f %f /\n", x, y-1, x+1, y-1);
            }
            else if (m == k + 1) {
               fprintf (fp, "\\plot %f %f  %f %f /\n", x+1, y, x+1, y-1);
            }
            else if (m == k - 1) {
               fprintf (fp, "\\plot %f %f  %f %f /\n", x, y, x, y-1);
            }

         }

      }

   }

   // Loop through the regions.
   for (n = 1; n <= puzzle[0][0]; n++) {

      for (k = 1; k <= 49; k++) {

         if (Region (k) == n) {

            // Compute the coordinates of square k's upper-left corner.
            x = (k-1) % 7;
            y = 7 - (k-1) / 7;
            if (op[n] == '+' || op[n] == '-') {
               fprintf (fp, "\\put {$\\scriptscriptstyle %c$} [cr] at %f %f\n", op[n], x+.9, y-.17);
            }
            else if (op[n] == '/') {
               fprintf (fp, "\\put {$\\scriptscriptstyle \\div$} [cr] at %f %f\n", x+.9, y-.17);
            }
            else {
               fprintf (fp, "\\put {$\\scriptscriptstyle \\times$} [cr] at %f %f\n", x+.9, y-.17);
            }
            fprintf (fp, "\\put {$\\scriptscriptstyle %d$} [cl] at %f %f\n", ans[n], x+.1, y-.17);
            break;

         }

      }

   }   

      
   fclose (fp);
   fclose (fps);

   printf ("View the puzzle and solution using plain TeX with KK.tex.\n");
   Exit();


}

////////////////////////////////////////////////////////////////////////////////
// Allocate array space.
////////////////////////////////////////////////////////////////////////////////
void AllocateMemory () {

   int i;

   prob = (double *) calloc (11, sizeof (double));
   x1 = (int *) calloc (50, sizeof (int));
   N = (int **) calloc (50, sizeof (int *));
   puzzle = (int **) calloc (50, sizeof (int *));
   for (i = 0; i <= 49; i++) {
      puzzle[i] = (int *) calloc (50, sizeof (int));
      N[i] = (int *) calloc (9, sizeof (int));
   }
   ans = (int *) calloc (50, sizeof (int));
   op  = (char *) calloc (50, sizeof (char));

}


