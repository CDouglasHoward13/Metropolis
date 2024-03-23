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
// Metropolis 5x5 crossword puzzle constructor as explained in Section 14.
////////////////////////////////////////////////////////////////////////////////

// Some global variables known to all functions below.
int    N;      // Number of words in the dictionary.
int   *first;  // Where words beginning with each letter start in the dictionary.
int   *last;   // Where words beginning with each letter end in the dictionary.
char **words;  // The words in the dictionary.
char **X;      // The current configuration as the Markov chain evolves.
double AP[4];  // Proposed transition acceptance probabilities when deltaE > 0.
int    i0, j0; // Proposed site to change.
char   X0;     // Proposed new letter.

// Functions found below.
void   Initialize ();
void   Metropolis (int);
int    Energy ();
int    BlackSquares ();
int    NotSymmetric ();
int    IsAWord (char, char, char, char, char);
int    Twice ();
void   MakeOutputFiles ();
char   RandomLetter ();
void   Proposal ();

#include "MetropolisFunctions.h"

////////////////////////////////////////////////////////////////////////////////
// Main program.
////////////////////////////////////////////////////////////////////////////////
int main () {

   int n;

   // Read in the dictionary, etc.  Dictionary5.txt must be present
   //   in the same directory as the executable version of this software.
   Initialize ();

   // Generate 10 random 5x5 puzzle solutions.
   for (n = 1; n <= 10; n++) {
      Metropolis (n);
      MakeOutputFiles ();
      Pause ();
   }

}

////////////////////////////////////////////////////////////////////////////////
// This function generates a random 5x5 crossword puzzle solution via Metropolis.
////////////////////////////////////////////////////////////////////////////////
void Metropolis (int n) {

   int i, j, E, deltaE, AcceptTransition;
   double t;

   printf ("\n");
   printf ("I'm working on puzzle number %d . ", n);

   // Generate a random initial configuration.
   for (i = 1; i <= 5; i++) {
      for (j = 1; j <= 5; j++) {
         X[i][j] = RandomLetter ();
      }
   }

   // Calculate the initial energy.
   E = Energy ();

   // This is to report progress to the screen.
   t = Time() + 1.0;

   // Continue until a solution is found (E == 0).
   while (E > 0) {

      // Show that Metropolis is still working on it every 2 seconds.
      if (Time() > t) {
         printf (". ");
         t = Time() + 2.0;
      }

      // Make a random change to the current puzzle...
      Proposal ();

      // Compute the change in energy.
      deltaE = Energy () - E;

      // Keep change and update the energy if deltaE <= 0.
      AcceptTransition = 0;
      if (deltaE <= 0) {
         AcceptTransition = 1;
      }

      // Otherwise keep with probability AP[deltaE] = e^(-deltaE/T).
      // Here T hardwired at 0.137 which was empirically found to be best.
      else {

         // Keep the proposed transition if accepted.
         if (MTUniform () <= AP[deltaE]) {
            AcceptTransition = 1;
         }

      }

      // Update the energy if transition is accepted.
      if (AcceptTransition == 1) {
         E += deltaE;
      }

      // Otherwise restore the current configuration, so X_{n+1} = X_n, and
      //    leave the energy unchanged.
      else {
         X[i0][j0] = X0;
      }

   }

   return;

}

////////////////////////////////////////////////////////////////////////////////
// Generate proposed random change. 
////////////////////////////////////////////////////////////////////////////////
void Proposal () {

   // i0, j0, and X0 are global variables.

   // Choose a random grid site (i0,j0) at random using the Mersenne Twister.
   // Will have i0 and j0 uniformly distributed on {1,2,3,4,5} --- independently.
   i0 = RandomInteger (1, 5);
   j0 = RandomInteger (1, 5);

   // ...Save the letter at that site.
   X0 = X[i0][j0];

   // ...Now randomly change it.
   while (1) {
      X[i0][j0] = RandomLetter ();
      // Accept this random letter only if it is different from the current
      //    letter (usually is).
      if (X[i0][j0] != X0) {
         break;
      }
   }

   return;
}

////////////////////////////////////////////////////////////////////////////////
// Compute the energy of the current configuration X[][].
////////////////////////////////////////////////////////////////////////////////
int Energy () {

   int i, E;

   // Initialize.
   E = 0;

   // Add one for each non-word in the configuration.
   for (i = 1; i <= 5; i++) {

      // Check the i^th row of the configuration. Augment energy if it is not a
      //    word.
      E += !IsAWord (X[i][1], X[i][2], X[i][3], X[i][4], X[i][5]);

      // Check the i^th column of the configuration.
      E += !IsAWord (X[1][i], X[2][i], X[3][i], X[4][i], X[5][i]);

   }

   // Uncomment this if you want only 5-letter words (no 4-letter words).
   // E += BlackSquares ();

   // (1) Comment out if some word appearing twice in the configuration is OK.
   E += Twice ();

   // (2) Uncomment if you want only puzzles that are symmetric about the main diagonal.
   // E += NotSymmetric ();

   // If (2) is uncommented you must comment (1), because words appear twice in
   //   symmetric puzzles.

   return E;

}

////////////////////////////////////////////////////////////////////////////////
// This function:
// (1) Seeds the random number generator;
// (2) Reads in the dictionary of five-letter words;
// (3) Allocates space for the Markov chain configuration X[][];
// (4) Determines where in the dictionary words beginning with each letter
//      are located.  For example: HABIT is the first word beginning with
//      H and is word number 4607; HYPOS is the last word beginning with H
//      and is word number 4741. This data helps when looking up a potential
//      word beginning with H (see the function IsAWord, where this is used).
// (5) Calculates acceptance probabilies when deltaE > 0.  This speeds up
//      simulations.
////////////////////////////////////////////////////////////////////////////////
void Initialize () {

   int i, k, n, deltaE;
   char input[100];
   FILE *fp;

   printf ("I will generate 10 random solutions to 5x5 crossword puzzles.\n\n");
   printf ("For each puzzle I will generate a TeX file (Puzzle.tex) which, when\n");
   printf ("processed with Plain TeX, will generate a beautiful puzzle for you!\n");

   // (1)
   // Seed the Mersenne Twister.
   MTUniform ();

   // (2)
   // Open the dictionary. This dictionary of five letter words was compiled by
   //    Donald Knuth.  I have added some four letter words. --CDH
   fp = fopen ("Dictionary5.txt", "r");

   // Dictionary can hold up to 10000 words (Dictionary5.txt currently has 9408).
   words = (char **) calloc (11001, sizeof (char *));

   // Read through the dictionary.  "N" records how may words are in the
   //   dictionary.
   N = 0;
   while (1) {

      // Get the next word in the dictionary; break "while (1)" loop when done.
      if (fgets (input, 99, fp) == NULL || input[0] == ' ' || input[0] == '\n') {
         break;
      }

      // Words preceded by a "-" are deleted from the dictionary. Otherwise,
      //    add the current word to the list.
      if (input[0] != '-') {

         // Augment size the size of the dictionary.
         N++;

         // Allocate memory for the next word.
         words[N] = (char *) calloc (6, sizeof (char));

         // Record the word.
         for (k = 0; k <= 4; k++) {
            words[N][k] = input[k];
         }
         words[N][5] = '\0';

      }

   }

   // Close the dictionary.
   fclose (fp);

   // (3)
   // Allocate space for the configuration array.
   X = (char **) calloc (6, sizeof (char *));
   for (i = 0; i <= 5; i++) {
      X[i] = (char *) calloc (6, sizeof (char));
   }

   // (4)
   // See where in the dictionary words beginning with each letter are first and
   //  last found. This helps in looking up a word (see the function IsAWord).

   // For example, the first word beginning with 'H' is HABIT, which is word
   //   number 4607 in the dictionary.  The last word beginning with 'H' is HYPOS
   //   (word number 4941).  The code below will calculate first['H'] = 46076
   //   and last['H'] = 4941.  This facilitates searching for 5 letter strings
   //   beginning with 'H' in the function IsAWord ().

   // (Note: in the ASCII code system, '@' = 64, 'A' = 65, 'B' = 66, ... , 'Z' = 90.)

   // This data is not "hard-wired" because the user has the option to delete
   //   words in the dictionary by placing a '-' in front of the word and to
   //   add words to the dictionary.

   // First allocate array space.
   first = (int *) calloc ('Z' + 1, sizeof (int));
   last  = (int *) calloc ('Z' + 1, sizeof (int));

   // Now loop throught the alphabet.
   for (k = '@'; k <= 'Z'; k++) {

      // Initial values. If no word begins with the letter "k" these values
      //   are not changed below.
      first[k] = 0;
      last[k]  = -1;

      // Loop through words in the dictionary.
      for (n = 1; n <= N; n++) {

         // If have not yet found any words beginning with "k" and the current
         //   word does begin with "k"...
         if (first[k] == 0 && words[n][0] == k) {
            first[k] = n;
         }

         // If have found words beginning with "k" but have not yet found any
         //   subsequent words not beginning with "k" and the current word does
         //   not begin with "k"...
         if (first[k] != 0 && last[k] == -1 && words[n][0] != k) {
            last[k] = n - 1;
            break;
         }

      }

      // If found words beginning with "k" but no subsequent words not beginning
      //   with "k" (e.g., in full dictionary when k = 'Z')...
      if (first[k] != 0 && last[k] == -1) {
         last[k] = N;
      }

   }

   // (5)
   // Calculate the probability of accepting proposed transitions when the
   //   change in energy is positive.  This uses that, when positive, deltaE can
   //   only take on the values 1, 2, and 3. The temperature is hard-wired at 0.137.
   // Calculation of these values ahead of time speeds up simulation of the Markov
   //   chain.
   for (deltaE = 1; deltaE <= 3; deltaE++) {
      AP[deltaE] = exp (-deltaE / 0.137);
   }

   return;

}

////////////////////////////////////////////////////////////////////////////////
// See if there are any blackened squares in the puzzle.
////////////////////////////////////////////////////////////////////////////////
int BlackSquares () {

   int i, j;

   for (i = 1; i <= 5; i++) {
      for (j = 1; j <= 5; j++) {
         if (X[i][j] == '@') return 1;
      }
   }

   return 0;

}

////////////////////////////////////////////////////////////////////////////////
// See if the puzzle is symmetric about the main diagonal. Return 1 if no, 0
//    if yes.
////////////////////////////////////////////////////////////////////////////////
int NotSymmetric () {

   int i, j;

   for (i = 2; i <= 5; i++) {
      for (j = 1; j < i; j++) {
         if (X[i][j] != X[j][i]) return 1;
      }
   }

   return 0;

}

////////////////////////////////////////////////////////////////////////////////
// Returns 1 if word "abcde" (see function arguments) is found in the
//    dictionary, 0 if it is not.
////////////////////////////////////////////////////////////////////////////////
int IsAWord (char a, char b, char c, char d, char e) {

   int n, found=0;

   // Look through words beginning with "a" to see if "abcde" is found.
   for (n = first[a]; n <= last[a] && !found; n++) {

      // If the n^th word in the dictionary is "abcde"...
      if (words[n][0] == a && words[n][1] == b
                           && words[n][2] == c
                           && words[n][3] == d
                           && words[n][4] == e) {

         found = 1;

      }

   }

   return found;

}

////////////////////////////////////////////////////////////////////////////////
// This function returns 1 if a word appears twice in a configuration, 0 if not.
////////////////////////////////////////////////////////////////////////////////
int Twice () {

   int i, j, k, t;

   // t = 1 means some word appears twice in X[][].
   // t = 0 means no word appears twice. Initially t is 0. Then checking begins.
   t = 0;

   // First check rows "i" and "j".
   for (i = 1; i < 5 && !t; i++) {
      for (j = i + 1; j <= 5 && !t; j++) {
         t = 1;
         for (k = 1; k <= 5 && t; k++) {
            t &= (X[i][k] == X[j][k]);
         }
      }
   }

   // Now check columns "i" and "j" (if rows were OK).
   for (i = 1; i < 5 && !t; i++) {
      for (j = i + 1; j <= 5 && !t; j++) {
         t = 1;
         for (k = 1; k <= 5 && t; k++) {
            t &= (X[k][i] == X[k][j]);
         }
      }
   }

   // Now check row "i" against column "j" (if rows and columns were OK).
   for (i = 1; i <= 5 && !t; i++) {
      for (j = 1; j <= 5 && !t; j++) {
         t = 1;
         for (k = 1; k <= 5 && t; k++) {
            t &= (X[i][k] == X[k][j]);
         }
      }
   }

   return t;

}

////////////////////////////////////////////////////////////////////////////////
// This function generates an output file of the puzzle.
////////////////////////////////////////////////////////////////////////////////
void MakeOutputFiles () {

   FILE *fp;
   int i, j;
   double c[] = {0, 0.1, 0.3, 0.5, 0.7,.9}, x, y, d = .004;

   // On-screen output:
   printf ("\n\n");
   printf ("Done! Here is the puzzle:\n\n");
   for (i = 1; i <= 5; i++) {
      for (j = 1; j <= 5; j++) {
         if (X[i][j] == '@') printf ("  ");
         else                printf (" %c", X[i][j]);
      }
      printf ("\n");
   }
   printf ("\n\n");
   printf ("For a better rendition, please view it with Puzzle.tex.\n\n");

   // Open the output file for the TeX version.
   fp = fopen ("Letters.txt", "w");

   for (i = 1; i <= 5; i++) {
      for (j = 1; j <= 5; j++) {
         x = c[j];
         y = c[6-i];
         if (X[i][j] == '@') {
            X[i][j] = ' ';
            fprintf (fp, "\\putrectangle corners at %5.3f %5.3f and %5.3f %5.3f\n",
               x+0.1, y+0.1, x-0.1, y-0.1);
            fprintf (fp, "\\plot %5.3f %5.3f  %5.3f %5.3f  %5.3f %5.3f  %5.3f %5.3f  %5.3f %5.3f /\n",
               x-0.1+d, y-0.1+d, x+0.1-d, y-0.1+d, x+0.1-d, y+0.1-d, x-0.1+d, y+0.1-d, x-0.1+d, y-0.1+d);
         }
         fprintf (fp, "\\put {\\bf %1c} at %5.3f %5.3f\n", X[i][j], x, y);
      }
   }

   fclose (fp);

   return;

}

///////////////////////////////////////////////////////////////////////////////////
// Generate a letter randomly. Uses ASCII convention: '@'=64, 'A'=65, ..., 'Z'= 90
///////////////////////////////////////////////////////////////////////////////////
char RandomLetter () {

   return (char) RandomInteger (64, 90);

}





