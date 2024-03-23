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
// This program creates an image and degrades it for a specified number of
// years as described in Section 16.
// Output files OriginalImagePixels.txt and DegradedImagePixels.txt are created.
// DegradedImagePixels.txt is the data input for ImageReconstruction.cpp.
////////////////////////////////////////////////////////////////////////////////

// These functions are found below.
void AllocateImageMemory ();
void MakeImage (int);
void ReportImage (int, int);

// These functions are common to all applications.
#include "MetropolisFunctions.h"

// The image is a global variable.
int **x;

////////////////////////////////////////////////////////////////////////////////
// Main program.
////////////////////////////////////////////////////////////////////////////////
int main() {

   int i, j, n, t, which;

   // Allocate array space for the images.
   AllocateImageMemory ();

   // Get information from the user.
   printf ("I will randomly degrade an image for you.\n\n");
   which = GetInteger ("Which image should I use: bull's eye (1) or smiley face (2)?... ");
   n = GetInteger ("\nHow many years of degradation should I do (100 <= years <= 2000)?... ");
   MTUniform(); // Seed the RNG.

   // Made the desired image.
   MakeImage (which);

   // Report the undegraded image.
   ReportImage (which, n);

   // Degrade the image over time, flipping 40 pixels per year.
   for (t = 1; t <= 40*n; t++) {

      i = RandomInteger (1, 200);
      j = RandomInteger (1, 200);

      // Flip the pixel at (i,j).
      x[i][j] = 1 - x[i][j];

   }

   // Report the degraded image.
   ReportImage (which, n);

   printf ("\nView the image degradation with ShowImageDegredation.tex using Plain TeX.\n");

   Pause ();

}

////////////////////////////////////////////////////////////////////////////////
// Report the current state of the image.
////////////////////////////////////////////////////////////////////////////////
void ReportImage (int which, int n) {

   int i, j;
   static int k = 0;
   FILE *fp;

   // Report the original or degraded image black pixels.
   if (k == 0) {
      fp = fopen ("OriginalImagePixels.txt", "w");
      if (which == 1) {
         fprintf (fp, "%% Bull's eye image undegraded.\n");
      } else {
         fprintf (fp, "%% Smiley face image undegraded.\n");
      }
   }

   else {
      fp = fopen ("DegradedImagePixels.txt", "w");
      fprintf (fp, "%% %d Years of degradation.\n", n);
   }


   for (i = 1; i <= 200; i++) {
      for (j = 1; j <= 200; j++) {
         if (x[i][j]) {
            fprintf (fp, "%d  %d\n", i, j);
         }
      }
   }
   fclose (fp);

   k = 1;

   return;

}

////////////////////////////////////////////////////////////////////////////////
// Generate the desired image and put it in the "x[i][j]" array. ///////////////
////////////////////////////////////////////////////////////////////////////////
void MakeImage (int which) {

   int i, j;
   double x0, y0, d;

   // Make the selected image; x[i][j] = 1 if black, 0 if white.
   for (i = 1; i <= 200; i++) {
      for (j = 1; j <= 200; j++) {

         // (x0,y0) are the coordinates of pixel (i,j).  They are scaled to
         //   fall between -1 and +1.
         x0 = (i-100) / 100.0;
         y0 = (j-100) / 100.0;

         // Bull's eye.
         if (which == 1) {
            d = sqrt(x0*x0 + y0*y0);
            if (d <= .199 || (.4 <= d && d <= .599) || (.8 <= d && d <= .999) ) {
               x[i][j] = 1;
            }
            else   {
               x[i][j] = 0;
            }
         }

         // Smiley face.
         else {
            d = sqrt(x0*x0 + y0*y0);
            // Circle around face:
            if (0.85 <= d && d <= 0.999) {
               x[i][j] = 1;
            }
            // Right eye:
            else if (sqrt(pow(x0-.5,2) +pow(y0-0.25,2)) < .129) {
               x[i][j] = 1;
            }
            // Left eye:
            else if (sqrt(pow(x0+.5,2) +pow(y0-0.25,2)) < .129) {
               x[i][j] = 1;
            }
            //Smile:
            else if ( (sqrt(pow(x0,2) +pow(y0-(-0.10),2)) < .5) &&
                      (sqrt(pow(x0,2) +pow(y0-0.15,2)) > .499) &&
                      (y0 < -0.10)
                    ) {
               x[i][j] = 1;
            }
            else   {
               x[i][j] = 0;
            }
         }

      } // "j" loop.
   } // "i" loop.

   return;

}

////////////////////////////////////////////////////////////////////////////////
// Allocate array space for the images /////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
void AllocateImageMemory () {

   // Allocate space for a 200 x 200 pixel image plus a white boundary (top,
   //    bottom, left, and right).
   // The white boundary is when i or j are 0 or 201. This is not part of the image.

   int i;

   x = (int **) calloc (202, sizeof (int));
   for (i = 0; i <= 201; i++) {
      x[i] = (int *) calloc (202, sizeof (int));
   }

   return;

}








