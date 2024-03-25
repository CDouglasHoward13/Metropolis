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
// This code does clustering analysis for Math and Reading/Writing SAT scores
// for 500 students as explained in Section 19.  It partitions the students
// into 20 clusters of students having similar SAT scores.
////////////////////////////////////////////////////////////////////////////////

// Global variables.
int N, K, i0, to, from, *cluster, *best, *count;
double  *x, *y, *Xsum, *Xbar, *X2sum, *Ysum, *Ybar, *Y2sum, *var, Emin;

// These functions are common to all applications.
#include "MetropolisFunctions.h"
                                                                                    
// These functions are all found below.
void   AllocateMemory ();
void   Report ();
void   ComputeSums ();
int    ClosestCenter (int);
double Distance2 (double, double, double, double);
double Energy ();
void   Proposal ();
void   Restore ();
void   UpdateData (int, int, int);
void   Metropolis ();
void   GetData ();
void   CopyConfiguration (int *, int *);
void   VoronoiTessellation ();
void   Prune (double *, double *, double *, double *, double, double, double, double);


////////////////////////////////////////////////////////////////////////////////
// Main program.
////////////////////////////////////////////////////////////////////////////////
int main() {

   // There 500 are students in the data set.
   N = 500;

   // Seeking to split them into 20 clusters with similar SAT scores.
   K = 20;

   // Allocate array space.
   AllocateMemory ();

   // Read in the data set.
   GetData ();

   // Apply Metropolis and report the results
   Metropolis ();

}

////////////////////////////////////////////////////////////////////////////////
// Apply the Metropolis algorithm.
////////////////////////////////////////////////////////////////////////////////
void Metropolis () {

   double deltaE, T, T0, E, t, t1, tmin;
   int AcceptTransition;

   // Initial temperature. This is much higher than for other applications,
   //    because the energy function is much greater (minimal is around 623,000).
   T0 = 300.0;

   // Record initial energy. At this stage the initial energy is the minimum
   //    ovserved so far. tmin records the time of the lowest energy-so-far.
   ComputeSums ();
   tmin = 0;
   Emin = E = Energy ();

   // Seed the RNG.
   printf ("I'm partitioning the SAT data into 20 clusters of like scores.\n");
   MTUniform ();

   printf ("\nI'll be done in 60 seconds. ");
   t = t1 = Time ();

   // Run the Markov chain for 60 seconds.
   while (t < 60.0) {

      // Every five seconds indicate that it's still thinking.
      t = Time ();
      if (t > t1 + 5.0) {
         printf (". ");
         t1 = t;
      }

      // Simulated annealing temperature schedule. T goes linearly
      //   from T0 to 0 as time elapses from 0 to 60 seconds.
      T = T0 * (60.0-t)/60.0;

      // Propose a change to the configuration.
      Proposal ();

      // Compute the change in energy.
      deltaE = Energy () - E;

      // Determine if proposed transition is accepted.
      AcceptTransition = 0;

      if (deltaE <= 0) {
         AcceptTransition = 1;
      }

      // So far it's just the zero temperature dynamics.
      // If T > 0, accept a transition to higher energy with some probability < 1.
      else if (T > 0) {
         if (MTUniform () <= exp (-deltaE/T)) {
            AcceptTransition = 1;
         }
      }

      // Keep the changed configuration and update the energy if the
      //   proposed transition is accepted. If the new energy is the best-so-far,
      //   update Emin and record the configuration.
      if (AcceptTransition) {
         E += deltaE;
         if (E < Emin) {
            Emin = E;
            CopyConfiguration (cluster, best);
            tmin = t;
         }
      }

      // Otherwise, restore the configuration.
      else {
         Restore ();
      }

   }

   // Time is up.  Copy the best configuration back into cluster[].
   // Re-compute and report the lowest-found energy.
   CopyConfiguration (best, cluster);
   ComputeSums ();
   Emin = Energy ();
   printf ("\n\n");
   printf ("At time %.1f, best-found energy is %.0f\n\n", tmin, Emin);
   printf ("View the clustering results in Clusters.txt.\n\n");
   printf ("View a scatter plot using plain TeX with ScatterPlot.tex.\n");

   // Report the Metropolis clusters.
   Report ();

}

////////////////////////////////////////////////////////////////////////////////
// Read in the data set.
////////////////////////////////////////////////////////////////////////////////
void GetData () {

   int i, k;
   char input[100];
   FILE *fp;

   // Read in the data set:
   fp = fopen ("SATs.txt", "r");
   for (i = 1; i <= N; i++) {

      // Student number i's scores:
      fgets (input, 99, fp);
      sscanf (input, "%lf %lf", x+i, y+i);

      // Assign a cluster to student number i:
      k = cluster[i] = (int) ((x[i]+y[i]-400)/60) + 1;

      // Increment that cluster's student count.
      count[k] ++;

   }

   fclose (fp);


}
////////////////////////////////////////////////////////////////////////////////
// Compute cluster counts and sums.
////////////////////////////////////////////////////////////////////////////////
void ComputeSums () {

   int i, k;

   // Initialize to zero.
   for (k = 1; k <= K; k++) {
      count[k] = Xsum[k] = X2sum[k] = Ysum[k] = Y2sum[k] = 0;
   }

   // Calculate population and coordinate sums for each cluster.
   for (i = 1; i <= N; i++) {

      k = cluster[i];  // This is the cluster that student number i is in.

      count[k] += 1;

      Xsum[k]  += x[i];

      X2sum[k] += x[i]*x[i];

      Ysum[k]  += y[i];

      Y2sum[k] += y[i]*y[i];

   }


}


////////////////////////////////////////////////////////////////////////////////
// Propose a change to the current configuration.
////////////////////////////////////////////////////////////////////////////////
void Proposal () {


   // Pick a student to move to a different cluster.

   i0 = RandomInteger (1, N);


   // Which cluster is it moving from?

   from = cluster[i0];


   // Pick a different cluster to move i0 to.

   to = from;

   while (to == from) {

      to = RandomInteger (1, K);

   }


   // Update data to reflect moving i0 from from to to. This function is

   //    20 times faster than updating the sums with ComputeSums ().

   UpdateData (i0, from, to);


}


////////////////////////////////////////////////////////////////////////////////

// Restore the current configuration. Update relevant data.
////////////////////////////////////////////////////////////////////////////////
void Restore () {

   UpdateData (i0, to, from);

}


////////////////////////////////////////////////////////////////////////////////

// Move student i from cluster f to cluster t. Update related info.

////////////////////////////////////////////////////////////////////////////////

void UpdateData (int i, int f, int t) {


   // t is i's new cluster.

   cluster[i] = t;


   // Update cluster counts.

   count[t] ++;

   count[f] --;


   // Decrement cluster f's data counts.

   Xsum[f]  -= x[i];
   X2sum[f] -= x[i]*x[i];
   Ysum[f]  -= y[i];
   Y2sum[f] -= y[i]*y[i];

   // Increment cluster t's data counts.
   Xsum[t]  += x[i];
   X2sum[t] += x[i]*x[i];
   Ysum[t]  += y[i];
   Y2sum[t] += y[i]*y[i];

}

////////////////////////////////////////////////////////////////////////////////

// Report the clusters to Clusters.txt.
////////////////////////////////////////////////////////////////////////////////
void Report () {

   int i, k;
   FILE *fp;

   fp = fopen ("Clusters.txt", "w");
   for (k = 1; k <= K; k++) {
      // Cluster k's center is at...
      Xbar[k] = Xsum[k]/count[k];
      Ybar[k] = Ysum[k]/count[k];
      // Now report the cluster.
      fprintf (fp, "\n");
      fprintf (fp, "Cluster number %d, with center at %.1f %.1f:\n", k, Xbar[k], Ybar[k]);
      for (i = 1; i <= N; i++) if (cluster[i] == k) {
         fprintf (fp, "%.0f %.0f\n", x[i], y[i]);
      }
   }
   fclose (fp);

   // Cluster centers.
   fp = fopen ("Centers.txt", "w");
   for (k = 1; k <= K; k++) {
      fprintf (fp, "%8.2f %8.2f\n", Xbar[k], Ybar[k]);
   }
   fclose (fp);

   // Report any abnormalities to the screen.
   for (i = 1; i <= N; i++) {
      if (ClosestCenter(i) != cluster[i]) {
         printf ("Student %d is not in the best cluster.\n", i);
      }
   }

   // Generate the Voronoi tessellation showing the clusters.
   VoronoiTessellation ();

   Pause ();

}


////////////////////////////////////////////////////////////////////////////////

// Find the center closest to each student.

////////////////////////////////////////////////////////////////////////////////

int ClosestCenter (int i) {


   int kmin, k;

   double d2min, d2;


   d2min = 1e+20;

   for (k = 1; k <= K; k++) {

      d2 = Distance2 (x[i],y[i], Xbar[k],Ybar[k]); 

      if (d2 < d2min) {

         d2min = d2;

         kmin = k;

      }

   }


   return kmin;


}


////////////////////////////////////////////////////////////////////////////////

// Compute the squared-distance from (x1,y1) to (x2,y2).
////////////////////////////////////////////////////////////////////////////////
double Distance2 (double x1, double y1, double x2, double y2) {

   return pow(x1-x2,2) + pow(y1-y2,2);

}


////////////////////////////////////////////////////////////////////////////////

// Calculate the energy function.
////////////////////////////////////////////////////////////////////////////////
double Energy () {

   int k;
   double E = 0.0;

   for (k = 1; k <= K; k++) {

      // If cluster k is non-empty, compute its variance and add it to E.
      if (count[k] != 0) {
         E += X2sum[k] - Xsum[k]*Xsum[k]/count[k] + Y2sum[k] - Ysum[k]*Ysum[k]/count[k];
      }

      // Otherwise, augment the energy by a huge number.
      else {
         E += 1e+10;
      }

   }

   return E;

}

////////////////////////////////////////////////////////////////////////////////
// Copy list a into list b.
////////////////////////////////////////////////////////////////////////////////
void CopyConfiguration (int *a, int *b) {

   int i;

   for (i = 1; i <= N; i++) {
      b[i] = a[i];
   }

}

////////////////////////////////////////////////////////////////////////////////
// Allocate array space.
////////////////////////////////////////////////////////////////////////////////
void AllocateMemory () {

   // All lists below are initialized to 0s.
   Xsum    = (double *) calloc (K+1, sizeof (double));
   Xbar    = (double *) calloc (K+1, sizeof (double));
   X2sum   = (double *) calloc (K+1, sizeof (double));
   Ysum    = (double *) calloc (K+1, sizeof (double));
   Ybar    = (double *) calloc (K+1, sizeof (double));
   Y2sum   = (double *) calloc (K+1, sizeof (double));
   var     = (double *) calloc (K+1, sizeof (double));
   x       = (double *) calloc (N+1, sizeof (double));
   y       = (double *) calloc (N+1, sizeof (double));
   count   = (int *)    calloc (K+1, sizeof (int));
   cluster = (int *)    calloc (N+1, sizeof (int));
   best    = (int *)    calloc (N+1, sizeof (int));

}



// THE TWO FUNCTIONS BELOW HERE CREATE THE VORONOI TESSELLATION GRAPH GENERATED BY
//   THE CENTERS {(Xbar[i],Ybar[i]): 1 <= i <= K}.


////////////////////////////////////////////////////////////////////////////////
// Generate the line segments forming the Voronoi tessellation.
////////////////////////////////////////////////////////////////////////////////
void VoronoiTessellation () {

   int i, j, k, point;
   double a, b, m, n1, n2, vx, vy, wx, wy;
   FILE *fp;

   // The algorithm below assumes the data lies in the square [n1,n2] x [n1,n2].
   // For this application, the data lies in the square [200,800] x [200,800], so...
   n1 = 200;
   n2 = 800;

   fp = fopen ("VoronoiGraph.txt", "w");

   // Loop through all pairs of points.
   // Find all Voronoi segments surrounding point i --- (Xbar[i],Ybar[i]).
   for (i = 1; i <= K; i++) {

      // Determine what portion of the line bisecting the line segment from point i to point j,
      //  if any, is such a segment.
      for (j = 1; j < i; j++) {

         // First compute the endpoints of the perpendicular bisector of the segment
         //   joining (Xbar[i],Ybar[i]) and (Xbar[j],Ybar[j]), restricted to the
         //   square [n1,n2] x [n1,n2].

         // (a,b) is the midpoint of that segment.
         a = (Xbar[i] + Xbar[j]) / 2;
         b = (Ybar[i] + Ybar[j]) / 2;

         // m is the slope of the perpendicular bisector.
         m = - (Xbar[i] - Xbar[j]) / (Ybar[i] - Ybar[j]);

         // Make the bisector vw lie between x=n1 and x=n2. Its equation is y = m*(x-a) + b.
         vx = n1;
         vy = m*(vx-a) + b;
         wx = n2;
         wy = m*(wx-a) + b;

         // Now make the bisector lie below y=n2. Keep only the portion of it that is
         //   closer to (0,0) than to (0,2*n2).
         Prune (&vx,&vy, &wx,&wy,  0,0, 0,2*n2);

         // Now make the bisector lie above y=n1. Keep only the portion of it that is
         //   closer to (0,n1+(n2-n1)) than to (0,n1-(n2-n1)).
         Prune (&vx,&vy, &wx,&wy,  0,n2, 0,2*n1-n2);

         // Loop through all points k different from i and j. Keep the portion of vw that
         //   is closer to point (x[i],y[i]) than to point (x[k],y[k]).
         point = 0;
         for (k = 1; k <= K && !point; k++) if (k != i && k != j) {
            Prune (&vx,&vy, &wx,&wy,  Xbar[i],Ybar[i], Xbar[k],Ybar[k]);
            if (Distance2 (wx,wy, vx,vy) < 0.00001) {
               point = 1;
            }
         }

         // If bisector segment vw has not been reduced to a point, report it to the
         //  output file.
         if (!point) {
            fprintf (fp, "\\plot %8.4f %8.4f  %8.4f %8.4f /\n", vx, vy, wx, wy);
         }

      }

   }
   fclose (fp);

}

////////////////////////////////////////////////////////////////////////////////
// This function "prunes" the line segment vw = ((vx,vy),(wx,wy)) updating the
//  endpoints to keep only the portion of vw that is closer to (x1,y1) than
//  to (x2,y2). Working with pointers (*vx, e.g.) also changes their values in
//  the upstairs function VoronoiTessellation().
////////////////////////////////////////////////////////////////////////////////
void Prune (double *vx, double *vy,
            double *wx, double *wy,
            double x1,  double y1,
            double x2,  double y2) {

   double m, mprime, a, b, x, y;

   // Find where the perpendicular bisector of the segment (x1,y1) to (x2,y2) intersects
   //    the line through v=(vx,vy) and w=(wx,wy). Call this point (x,y).

   // The equation of the perpendicular bisector is y = m*x + (b-m*a), where...
   a = (x1 + x2) / 2;
   b = (y1 + y2) / 2;
   m = - (x1 - x2) / (y1 - y2);

   // The equation of the line going through v and w is y = mprime*x + (vy - mprime*vx),
   //  where...
   mprime = (*wy - *vy) / (*wx - *vx);

   // Setting m*x + (b-m*a) = mprime*x + (vy-mprime*vx) gives...
   x = ((*vy - mprime*(*vx)) - (b-m*a)) / (m - mprime); 

   // And then...
   y = mprime*x + (*vy - mprime*(*vx));

   // Keep only the part of the segment joining v and w that is closer to (x1,y1) than (x2,y2).
   if (Distance2(*wx,*wy, x2,y2) < Distance2(*wx,*wy, x1,y1)) {
      *wx = x;
      *wy = y;
   }
   if (Distance2(*vx,*vy, x2,y2) < Distance2(*vx,*vy, x1,y1)) {
      *vx = x;
      *vy = y;
   }
   // If both v and w are closer to (x1,y1), there is no updating.
   // If both v and w are closer to (x2,y2), the updated segment becomes a point -- (x,y).

   return;

}



