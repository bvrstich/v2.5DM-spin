#ifndef dDPM_H
#define dDPM_H

#include <iostream>
#include <fstream>

using std::ostream;

#include "rTPM.h"

/**
 * @author Brecht Verstichel
 * @date 04-08-2011\n\n
 * This class dDPM is a class written for the I1 and Q2 conditions, which is a DPM matrix with the spacial third index diagonal.
 * The v2.5DM program will use this as the central variable. 
 */
class dDPM {

   /**
    * Output stream operator overloaded
    * @param output The stream to which you are writing (e.g. cout)
    * @param ddpm_p the dDPM you want to print
    */
   friend ostream &operator<<(ostream &output,const dDPM &ddpm_p);

   public:

      //constructor
      dDPM();

      //copy constructor
      dDPM(const dDPM &);

      //destructor
      virtual ~dDPM();

      rTPM &operator[](int);

      const rTPM &operator[](int) const;

      double operator()(int l,int S,int S_ab,int a,int b,int S_cd,int c,int d) const;

      static double get_inco(int &l,int S,int &S_ab,int &a,int &b);

      int gN() const;

      int gM() const;

      double trace() const;

      void dscal(double);

      dDPM &operator=(const dDPM &);

      dDPM &operator+=(const dDPM &);

      dDPM &operator-=(const dDPM &);

      dDPM &operator=(double);

      dDPM &daxpy(double,const dDPM &);

      double ddot(const dDPM &) const;

      void L_map(const dDPM &,const dDPM &);

      void sqrt(int);

      void invert();

      dDPM &mprod(const dDPM &,const dDPM &);

      void symmetrize();

      void fill_Random();

      void proj_W();

      void test_proj() const;

      static void init(int,int);

      static void clear();

   private:

      //!nr of particles
      static int N;

      //!dimension of sp hilbert space
      static int M;

      //!double array of rTPM objects
      rTPM **ddpm;

      //!static array holding the 6j symbols needed
      static double **_6j;

};

#endif
