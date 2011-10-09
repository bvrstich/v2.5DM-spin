#ifndef dPHHM_H
#define dPHHM_H

#include <iostream>
#include <fstream>

using std::ostream;

#include "rxPHM.h"

class SUP;

/**
 * @author Brecht Verstichel
 * @date 04-08-2011\n\n
 * This class dPHHM is a class written for the I2 and Q1 conditions, which is a PPHM matrix with the spatial third index diagonal.
 */
class dPHHM {

   /**
    * Output stream operator overloaded
    * @param output The stream to which you are writing (e.g. cout)
    * @param ddpm_p the dPHHM you want to print
    */
   friend ostream &operator<<(ostream &output,const dPHHM &ddpm_p);

   public:

      //constructor
      dPHHM();

      //copy constructor
      dPHHM(const dPHHM &);

      //destructor
      virtual ~dPHHM();

      rxPHM &operator[](int);

      const rxPHM &operator[](int) const;

      double operator()(int,int,int,int,int,int,int,int) const;

      static int get_inco(int,int,int,int);

      int gN() const;

      int gM() const;

      double trace() const;

      void dscal(double);

      dPHHM &operator=(const dPHHM &);

      dPHHM &operator+=(const dPHHM &);

      dPHHM &operator-=(const dPHHM &);

      dPHHM &operator=(double);

      dPHHM &daxpy(double,const dPHHM &);

      double ddot(const dPHHM &) const;

      void L_map(const dPHHM &,const dPHHM &);

      void sqrt(int);

      void invert();

      dPHHM &mprod(const dPHHM &,const dPHHM &);

      void symmetrize();

      void fill_Random();

      void G1(const dDPM &);

      void G2(const dDPM &);

      void test_sym() const;

      static void init(int,int);

   private:

      //!nr of particles
      static int N;

      //!dimension of sp hilbert space
      static int M;

      //!double array of rxPHM objects
      rxPHM **dphhm;

};

#endif
