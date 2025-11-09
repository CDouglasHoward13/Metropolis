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
// This code reads in the magic square in MSout.txt and verifies
// that it is, indeed, a magic square.
////////////////////////////////////////////////////////////////////////////////

char input[1000], *p;
int row, col, m, M[100][100], Numbers[10001], sum, good, order;
#include "MetropolisFunctions.h"

int main() {

   FILE *fp;

   // Initialize the Numbers list.
   for (m = 1; m <= 10000; m++) Numbers[m] = 0;

   fp = fopen ("MSout.txt", "r");

   // Get the magic square's order.
   fgets (input, 399, fp);
   sscanf (input, "%d", &order);
   printf ("Verifying the order %d magic square in MSout.txt.\n\n", order);

   // Read in the magic square row-by-row.
   for (row = 1; row <= order; row++) {

      // Get the next row from MSout.txt.
      fgets (input, 399, fp);
      p = input;
      for (col = 1; col <= order; col++) {
         // Next number.
         sscanf (p, "%d", &m);
         M[row][col] = m;
         // Record that the number m was found in the magic square.
         Numbers[m] += 1;
         // Scoot the pointer over to the next number in the row.
         p += 6;
      }

   }

   fclose (fp);

   // Verify that all numbers 1,2,...,order^2 are present in the square.
   good = 1;
   for (m = 1; m <= order*order; m++) {
      if (Numbers[m] != 1) good = 0;
   }

   if (good) {
      printf ("All numbers 1,2,...,%d are found in the magic square.\n", order*order);
      Pause ();
   }
   else {
      printf ("Some numbers are missing in the square!\n");
      Exit ();
   }

   printf ("Sum of rows:\n");
   for (row = 1; row <= order; row++) {
      sum = 0;
      for (col = 1; col <= order; col++) {
         sum += M[row][col];
      }
      printf ("%2d  %6d\n", row, sum);
   }
   Pause ();

   printf ("Sum of columns:\n");
   for (col = 1; col <= order; col++) {
      sum = 0;
      for (row = 1; row <= order; row++) {
         sum += M[row][col];
      }
      printf ("%2d  %6d\n", col, sum);
   }
   Pause ();

   printf ("Sum of two diagonals:\n");
   sum = 0;
   for (row = 1; row <= order; row++) {
      col = row;
      sum += M[row][col];
   }
   printf ("d1  %6d\n", sum);

   sum = 0;
   for (row = 1; row <= order; row++) {
      col = order+1 - row;
      sum += M[row][col];
   }
   printf ("d2  %6d\n", sum);
   Pause ();

}












