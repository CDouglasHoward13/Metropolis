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
// This code uses the Metropolis algorithm to reconstruct an image
//   that has been randomly degraded as described in Section 16.

////////////////////////////////////////////////////////////////////////////////

// Image arrays are global variables so that the above functions have access.
int **x,    // Degraded image, then reconstructed image.
    **d,    // Degraded image
    **best; // Lowest energy (best-found) image.

// E_min and lambda are used in the functions Metropolis() and Energy ().
double E_min, lambda;

// These functions are found below.
void   AllocateImageMemory (void);
void   ReportImage ();
void   CopyImage (int **, int **);
double Energy (void);
double DeltaEnergy (int, int);
void   GetDegradedImage (void);
void   Metropolis (void);

// These functions are common to all applications.
#include "MetropolisFunctions.h"

/////////////////////////////////////////////////////////////////////////////////////////
// Main program. ////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////
int main() {

   // Get the degraded image; put it in the "d" array and also the "best" array.
   GetDegradedImage ();

   // Report the degraded image.
   ReportImage ();

   // Reconstruct the image via Metropolis. Report the reconstruction process.
   Metropolis ();

}

/////////////////////////////////////////////////////////////////////////////////////////
// This function reads in the degraded image.   /////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////
void GetDegradedImage () {

   int i, j;
   double t, p, q;
   FILE *fp;
   char input[100];


   // Allocate array space for image pixels.
   AllocateImageMemory ();

   // Open the data file.
   fp = fopen ("DegradedImagePixels.txt", "r");

   // First line contains the years of degradation.
   fgets (input, 99, fp);
   sscanf (input+2, "%lf", &t);
   printf ("I'm reconstructing an image degraded %.0f years.\n", t);

   // Seed the RNG.
   MTUniform ();

   // Read in the blackened pixel coordinates in the degraded image.
   while (1) {

      if (fgets (input, 99, fp) == NULL) break;
      sscanf (input, "%d %d", &i, &j);
      d[i][j] = 1;

   }

   // Compute the probability "p" that any particular pixel is flipped
   //    at "t" years and the related quantity "lambda".
   p = 1.0/40000.0;
   q = 1.0 - p;
   p = (1.0 - pow(q-p, 40.0*t)) / 2.0;

   // This parameter is used in the energy computation, which is why "t" is
   //    needed as an input to Metropolis.
   lambda = 1.0 / (1.0 + log((1.0-p)/p));

   // Start the Markov chain in the degraded state.
   CopyImage (d, x);

   // Initially the degraded image itself is the best reconstruction so far.
   CopyImage (d, best);

}


////////////////////////////////////////////////////////////////////////////////
// Reconstruct the image via Metropolis.   /////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
void Metropolis () {

   int i0, j0, NextReport, AcceptTransition, n;
   double U, p, T, E, DeltaE, t, t1;

   // Get the temperature.
   T = GetDouble ("\nWhat is the temperature (best about 0.1)?... ");

   // Initalize the energy of the image E and the lowest energy found
   //    so far E_min.
   // Since we are not interested in the value of the minimal energy (only
   //    in the minimizing configuration), E can be set to any value.
   E_min = E = Energy ();

   // Report next image at Markov chain period number "NextReport".
   // Reports are generated at 1000, 10000, 100000, 1000000, 10000000, and
   //    100000000 steps of the Markov chain.
   NextReport = 1000;
   n = 0;

   printf ("\nI'll be done in 60 seconds. ");
   t = t1 = Time ();

   // Run the Markov chain for 60 seconds.
   // The image "x" will always be the current value of the Markov chain.
   while (t < 60.0) {

      // Increment Markov chain step counter.
      n ++;

      // Every five seconds indicate that it's still thinking.
      t = Time ();
      if (t > t1 + 5.0) {
         printf (". ");
         t1 = t;
      }

      // Generate proposed transition -- (i0,j0) is the proposed pixel to flip.
      i0 = RandomInteger (1, 200);
      j0 = RandomInteger (1, 200);

      // Compute the change in energy associated with a flip of that pixel...
      DeltaE = DeltaEnergy (i0, j0);

      // Start with the zero-temperature dynamics.
      AcceptTransition = 0;
      if (DeltaE <= 0) {
         AcceptTransition = 1;
      }

      // If T > 0 -- accept an increase in energy with the appropriate
      //   probability, called "p" below.
      else if (T > 0) {
         p = exp (-DeltaE / T);
         U = MTUniform ();
         if (U <= p) {
            AcceptTransition = 1;
         }   
      }

      // Effect the transition if indicated. Compute the new energy and if it
      //    is the best-found-so-far, record relevant information.
      if (AcceptTransition) {

         // Flip site (i0,j0).
         x[i0][j0] = 1 - x[i0][j0];

         // Update the current energy of the image.
         E += DeltaE;

         // Update lowest energy found so far, as appropriate.
         if (E < E_min) {

            // Record the new lowest energy found so far.
            E_min = E;

            // Record the new lowest energy image.
            CopyImage (x, best);

         }
         
      }

      // Periodically report the reconstructed image to an output file.
      if (n == NextReport) {

         // Report the image as currently reconstructed.
         ReportImage ();

         // Update the next image to be reported.
         NextReport *= 10;

      }

   } // This ends the Markov chain for loop.

   // Report the lowest energy reconstruction.
   ReportImage ();

   // Finish up.
   printf ("\n\n");
   printf ("%.1f million Markov chain steps completed in 60 seconds.\n\n", n/1000000.0);
   printf ("View the reconstruction process with ShowImageReconstruction.tex using Plain TeX.\n");

   // Pause and exit program.
   Exit ();

}

////////////////////////////////////////////////////////////////////////////////
// Allocate array space for the images /////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
void AllocateImageMemory () {

   // Allocate space for a 200 x 200 pixel image plus boundary.
   // Make the image all white -- this happens by default because
   // calloc() initializes arrays to 0. Allocate for x, d, and best.

   int i;

   // "x" image.
   x = (int **) calloc (202, sizeof (int));
   for (i = 0; i <= 201; i++) {
      x[i] = (int *) calloc (202, sizeof (int));
   }

   // "d" image.
   d = (int **) calloc (202, sizeof (int));
   for (i = 0; i <= 201; i++) {
      d[i] = (int *) calloc (202, sizeof (int));
   }

   // "best" image.
   best = (int **) calloc (202, sizeof (int));
   for (i = 0; i <= 201; i++) {
      best[i] = (int *) calloc (202, sizeof (int));
   }

   return;

}



////////////////////////////////////////////////////////////////////////////////
// Report the specified image to the specified output file. The coordinates ////
// of all blackened pixels are reported on a [-1,1] x [-1,1] scale. ////////////
////////////////////////////////////////////////////////////////////////////////
void ReportImage () {

   static n = 0;
   int i, j, **image;
   FILE *fp;
   char filename[10][100] = {"DegradedImage.txt", "1000.txt", "10000.txt",
                             "100000.txt", "1000000.txt", "10000000.txt",
                             "100000000.txt", "BestReconstruction.txt"};

   // Usually report image "x", "best" is the last reported.
   image = (n < 7 ? x : best);

   // Open the appropriate output file.
   fp = fopen (filename[n], "w");

   // Now report the coordinates of the blackened pixels to a file.
   for (i = 1; i <= 200; i++) {
      for (j = 1; j <= 200; j++) {
         if (image[i][j]) {
            fprintf (fp, "%d %d\n", i, j);
         }
      }
   }
   fclose (fp);

   // Next image.
   n ++;

   return;

}

////////////////////////////////////////////////////////////////////////////////
// Copy image "from" into image "to" ///////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
void CopyImage (int **from, int **to) {

   int i, j;

   for (i = 1; i <= 200; i++) {
      for (j = 1; j <= 200; j++) {
         to[i][j] = from[i][j];
      }
   }

   return;

}


/////////////////////////////////////////////////////////////////////////////////
// Compute the energy of configuration "x". /////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
double Energy () {

   int i, j;
   double b, D = 0, B = 0;

   for (i = 1; i <= 200; i++) {
      for (j = 1; j <= 200; j++) {

         // First compute the number of boundary segments around (i,j).
         b  =    (x[i][j+1] != x[i][j])   // boundary with northern neighbor?
               + (x[i+1][j] != x[i][j])   // eastern neighbor?
               + (x[i][j-1] != x[i][j])   // southern neighbor?
               + (x[i-1][j] != x[i][j]);  // western neighbor?
         // Multiply by 1/2 to reflect double counting.
         B += 0.5 * b;

         // Now see if pixel (i,j) disagrees with the original image.
         D += (x[i][j] != d[i][j]);

      }
   }

   return lambda*B + (1.0-lambda)*D;

}

////////////////////////////////////////////////////////////////////////////////
// Compute the change in energy when pixel (i0,j0) in image x is flipped. //////
////////////////////////////////////////////////////////////////////////////////
double DeltaEnergy (int i0, int j0) {

   double b, DeltaD, DeltaB, DeltaE;

   // - First compute Delta D.
   DeltaD = (x[i0][j0] == d[i0][j0] ? 1 : -1);

   // - Now compute Delta B.
   b  =    (x[i0][j0+1] != x[i0][j0])   // boundary with northern neighbor?
         + (x[i0+1][j0] != x[i0][j0])   // eastern neighbor?
         + (x[i0][j0-1] != x[i0][j0])   // southern neighbor?
         + (x[i0-1][j0] != x[i0][j0]);  // western neighbor?
   DeltaB = 4 - 2*b;

   // - Now compute Delta E.
   DeltaE = lambda*DeltaB + (1.0-lambda)*DeltaD;

   return DeltaE;

}




