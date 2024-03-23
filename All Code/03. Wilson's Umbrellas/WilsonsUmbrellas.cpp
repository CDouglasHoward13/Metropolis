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
// This code simulates Wilson's umbrellas from Section 3, Example 1.
////////////////////////////////////////////////////////////////////////////////

#include "MetropolisFunctions.h"


int main() {

   int home, school, n, N, wet;
   double U, p;

   // Get the probability of foul weather.
   printf ("I am simulating Wilson's umbrellas from Section 3, Example 1.\n");

   // Seed the RNG.
   MTUniform();

   // Get the probability of rain.
   p = GetDouble ("\nWhat is the probability of rain?... ");

   // Initially all his umbrellas are at home and none at school.
   home = GetInteger("\nHow many umbrellas does Wilson have?... ");
   school = 0;

   // Count how many walks Wilson gets wet.
   wet = 0;

   // Simulate 100 million round trips.
   N = 100000000;

   // Initialize the timer.
   Time ();

   // Begin the simulations
   for (n = 1; n <= N; n++) {

      if (n == N/10) {
         printf ("\nShould be done in %.1f seconds.\n", 10.0*Time());
      }   

      // Home ---> school.
      U = MTUniform();
      // If it rains:
      if (U <= p) {
         if (home > 0) {
            home --;      // Take an umbrella from home to school
            school ++;    // and stay dry.
         } else {
            wet ++;       // Wilson gets wet.
         }
      }

      // School ---> home.
      U = MTUniform();
      // If it rains:
      if (U <= p) {
         if (school > 0) {
            school --;    // Take an umbrella from school to home
            home ++;      // and stay dry.
         } else {
            wet ++;       // Wilson gets wet.
         }
      }

   }

   // Report the results. In N round-trips Wilson did 2*N walks.
   printf ("\nFraction of times Wilson got wet is %.4f\n", wet/(2.0*N));

   // Pause so execution window doesn't close up.
   Pause();

}

