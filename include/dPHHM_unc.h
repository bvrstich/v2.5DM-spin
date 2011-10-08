#ifndef dPHHM_unc_H
#define dPHHM_unc_H 
#include <iostream>
#include <fstream>

using std::ostream;

#include "rxPHM_unc.h"

class dDPM_unc;

/**
 * @author Brecht Verstichel
 * @date 13-05-2011\n\n
 * This class dPHHM_unc is a class written for the I1 and Q2 conditions, which is a DPM matrix with the third indices diagonal.
 * The v2.5DM program will use this as the central variable. 
 */
class dPHHM_unc {

   /**
    * Output stream operator overloaded
    * @param output The stream to which you are writing (e.g. cout)
    * @param dphhm_p the dPHHM_unc you want to print
    */
   friend ostream &operator<<(ostream &output,const dPHHM_unc &dphhm_p);

   public:

      //constructor
      dPHHM_unc();

      //copy constructor
      dPHHM_unc(const dPHHM_unc &);

      //destructor
      virtual ~dPHHM_unc();

      rxPHM_unc &operator[](int);

      const rxPHM_unc &operator[](int) const;

      int gN() const;

      int gM() const;

      double trace() const;

      void dscal(double);

      dPHHM_unc &operator=(const dPHHM_unc &);

      dPHHM_unc &operator+=(const dPHHM_unc &);

      dPHHM_unc &operator-=(const dPHHM_unc &);

      dPHHM_unc &operator=(double);

      dPHHM_unc &daxpy(double,const dPHHM_unc &);

      double ddot(const dPHHM_unc &) const;

      void L_map(const dPHHM_unc &,const dPHHM_unc &);

      void sqrt(int);

      void invert();

      dPHHM_unc &mprod(const dPHHM_unc &,const dPHHM_unc &);

      void symmetrize();

      void fill_Random();

      void print_eig();

      void uncouple(const dPHHM &);

      void out_sp(const char *) const;

      static void init(int,int);

   private:

      //!nr of particles
      static int N;

      //!dimension of sp hilbert space
      static int M;

      //!double array of rxPHM_unc objects
      rxPHM_unc **dphhm;

};

#endif
