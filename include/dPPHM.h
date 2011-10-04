#ifndef dPPHM_H
#define dPPHM_H

#include <iostream>
#include <fstream>

using std::ostream;

#include "xTPM.h"

class SUP;

/**
 * @author Brecht Verstichel
 * @date 04-08-2011\n\n
 * This class dPPHM is a class written for the I2 and Q1 conditions, which is a PPHM matrix with the spatial third index diagonal.
 */
class dPPHM {

   /**
    * Output stream operator overloaded
    * @param output The stream to which you are writing (e.g. cout)
    * @param ddpm_p the dPPHM you want to print
    */
   friend ostream &operator<<(ostream &output,const dPPHM &ddpm_p);

   public:

      //constructor
      dPPHM();

      //copy constructor
      dPPHM(const dPPHM &);

      //destructor
      virtual ~dPPHM();

      xTPM &operator[](int);

      const xTPM &operator[](int) const;

      double operator()(int l,int S,int S_ab,int a,int b,int S_cd,int c,int d) const;

      static int get_inco(int S,int S_ab,int a,int b);

      int gN() const;

      int gM() const;

      double trace() const;

      void dscal(double);

      dPPHM &operator=(const dPPHM &);

      dPPHM &operator+=(const dPPHM &);

      dPPHM &operator-=(const dPPHM &);

      dPPHM &operator=(double);

      dPPHM &daxpy(double,const dPPHM &);

      double ddot(const dPPHM &) const;

      void L_map(const dPPHM &,const dPPHM &);

      void sqrt(int);

      void invert();

      dPPHM &mprod(const dPPHM &,const dPPHM &);

      void symmetrize();

      void fill_Random();

      void I(const dDPM &);

      void Q(const dDPM &);

      static void init(int,int);

   private:

      //!nr of particles
      static int N;

      //!dimension of sp hilbert space
      static int M;

      //!double array of xTPM objects
      xTPM **dpphm;

      //!static array holding the 6j symbols needed
      static double **_6j;

};

#endif
