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

// This file declares and defines various functions that are used
// in all Metropolis applications.  

// These are some standard C function libraries.
#include <stdlib.h>  // "standard library"
#include <math.h>    // various math functions, such as exp()
#include <stdio.h>   // various "input/output" functions
#include <time.h>    // functions for timing computations
#include <string.h>  // various string functions

// These functions are found below.
double MTUniform (void);
int    RandomInteger (int, int);
void   Pause (void);
void   Exit (void);
double Time (void);
int    GetInteger (const char *);
double GetDouble (const char *);

////////////////////////////////////////////////////////////////////////////////
// MERSENNE TWISTER
// By M. Matsumoto and T. Nishimura (1998).
// "Mersenne Twister: a 623-dimensionally equidistributed uniform pseudo-random
//   number generator".
// ACM Transactions of Modeling and Computer Simulation 8(1):3-30.
// Any coding errors introduced are my own (C.D. Howard).

// An unsigned integer is represented by 32 bits in base 2.  The largest
// unsigned integer is: 2^32 - 1 = 4,294,967,295 (base 10);
//          = ffffffff (base 16);
//          = 1111 1111 1111 1111 1111 1111 1111 1111 (base 2).
// The digits in hexadecimal (base 16) are 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, a, b,
//                                         c, d, e, f.
////////////////////////////////////////////////////////////////////////////////
double MTUniform () {

   static unsigned int X[1248], m[2], seeded = 0, k;
   unsigned int N, Y;

   // Seed the RNG when a new seed is passed or it has not yet been initialized.
   if (!seeded) {
      // If no seed is specified, default is 1.
      X[0] = GetInteger ("\nPlease seed the Mersenne Twister with a positive integer... ");
      // Now seed X[1],X[2],...,X[623] with your favorite LCG.
      for (k = 1; k < 624; k++) {
         X[k] = 22695477 * X[k-1] + 1;
      }
      m[0] = 0; m[1] = 0x9908b0df;
      // The counter "k" is now 624.
      seeded = 1;
   }

   // Concatenate the first bit of X[k-624] with the last 31 bits of X[k-623],
   //    "|" is "bit-wise or".
   Y = (X[k-624] & 0x80000000) | (X[k-623] & 0x7fffffff);

   // Now apply an invertible linear transformation Y and bit-wise add X[k-227].
   X[k] = ((Y >> 1) ^ m[Y & 1] ^ X[k-227]);

   // Re-load X[0],X[1],...,X[623] as you go.
   X[k-624] = X[k];

   // Apply the tempering function (also an invertible linear transformation).
   N = X[k];
   N ^= (N >> 11);
   N ^= (N << 7) & 0x9d2c5680;
   N ^= (N << 15) & 0xefc60000;
   N ^= (N >> 18);

   // Increment the counter; shift vectors when appropriate.
   k ++;
   if (k == 1248) k = 624;

   // Now 0 <= N <= 4,294,967,295; scale it to be on the interval (0,1).
   return ( (N + 0.5) / 4294967296.0 );

}

////////////////////////////////////////////////////////////////////////////////
// Generate an integer uniformly from among {a,...,b}. ///////////////////////
////////////////////////////////////////////////////////////////////////////////
int RandomInteger (int a, int b) {

   return ((int) ((b-a+1) * MTUniform () + a));

}


////////////////////////////////////////////////////////////////////////////////
// This function waits for a user-input "enter", then CONTINUES the program.
////////////////////////////////////////////////////////////////////////////////
void Pause () {

   char input[100];

   printf ("\n");
   printf ("Hit Enter to continue program... ");
   fgets (input, 9, stdin);

   return;

}

////////////////////////////////////////////////////////////////////////////////
// This function waits for a user-input "enter", then EXITS the program.
// It prevents the window from closing up before the output can be viewed.
////////////////////////////////////////////////////////////////////////////////
void Exit () {

   char input[100];

   printf ("\n");
   printf ("Hit Enter to exit program... ");
   fgets (input, 9, stdin);

   exit (0);

}

////////////////////////////////////////////////////////////////////////////////
// This function computes elapsed computation time in seconds.
////////////////////////////////////////////////////////////////////////////////
double Time () {

   static clock_t time;
   static int initialized = 0;

   // With the first call to this function, "time" is initialized.
   if (!initialized) {
      time = clock ();
      initialized = 1;
   }

   // With subsequent calls, elapsed time since the first call is returned.
   return ((double) (clock() - time)) / CLOCKS_PER_SEC;

}

////////////////////////////////////////////////////////////////////////////////
// This function gets an integer value typed in by the user.
////////////////////////////////////////////////////////////////////////////////
int GetInteger (const char *question) {

   char input[100];
   int n;

   // Print the question.
   printf (question);

   // Get the answer.
   fgets (input, 99, stdin);
   sscanf (input, "%d", &n);

   return (n);

}

////////////////////////////////////////////////////////////////////////////////
// This function gets a double precision value typed in by the user.
////////////////////////////////////////////////////////////////////////////////
double GetDouble (const char *question) {

   char input[100];
   double x;

   // Print the question.
   printf (question);

   // Get the answer.
   fgets (input, 99, stdin);
   sscanf (input, "%lf", &x);

   return (x);

}


