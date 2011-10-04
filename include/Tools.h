#ifndef TOOLS_H
#define TOOLS_H

#include <iostream>
#include <cstdlib>

/**
 * @author Brecht Verstichel
 * @date 26-09-2011\n\n
 * This is a class that contains some static objects (lists, matrices) needed by all the classes,
 * and it would be stupid to store it in every class.
 */

class Tools{

   public:

      static void init(int,int);

      static void clear();

      static double g3j(int,int,int,int,int,int);

      static double g6j(int,int,int,int,int,int);

      static double gC(int,int,int,int);

      static double gD(int,int,int,int);

      static dDPM &gunit();

   private:

      //!nr of sites
      static int M;

      //!nr of particles
      static int N;

      //!the projected unit matrix
      static dDPM *unit;

      //!3j symbols
      static double ******_3j;

      //!6j symbols
      static double ******_6j;

      //!symbol needed for I2 and Q2 condition, see notes
      static double ****C;

      //!symbol needed for Q2 condition, see notes
      static double ****D;

};

#endif
