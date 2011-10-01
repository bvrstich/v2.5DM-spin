#ifndef dPPHM_unc_H
#define dPPHM_unc_H 
#include <iostream>
#include <fstream>

using std::ostream;

#include "xTPM_unc.h"

class dDPM_unc;

/**
 * @author Brecht Verstichel
 * @date 13-05-2011\n\n
 * This class dPPHM_unc is a class written for the I1 and Q2 conditions, which is a DPM matrix with the third indices diagonal.
 * The v2.5DM program will use this as the central variable. 
 */
class dPPHM_unc {

   /**
    * Output stream operator overloaded
    * @param output The stream to which you are writing (e.g. cout)
    * @param dpphm_p the dPPHM_unc you want to print
    */
   friend ostream &operator<<(ostream &output,const dPPHM_unc &dpphm_p);

   public:

      //constructor
      dPPHM_unc();

      //copy constructor
      dPPHM_unc(const dPPHM_unc &);

      //destructor
      virtual ~dPPHM_unc();

      xTPM_unc &operator[](int);

      const xTPM_unc &operator[](int) const;

      int gN() const;

      int gM() const;

      double trace() const;

      void dscal(double);

      dPPHM_unc &operator=(const dPPHM_unc &);

      dPPHM_unc &operator+=(const dPPHM_unc &);

      dPPHM_unc &operator-=(const dPPHM_unc &);

      dPPHM_unc &operator=(double);

      dPPHM_unc &daxpy(double,const dPPHM_unc &);

      double ddot(const dPPHM_unc &) const;

      void L_map(const dPPHM_unc &,const dPPHM_unc &);

      void sqrt(int);

      void invert();

      dPPHM_unc &mprod(const dPPHM_unc &,const dPPHM_unc &);

      void symmetrize();

      void fill_Random();

      void print_eig();

      void uncouple(const dPPHM &);

      void out_sp(const char *) const;

      void I(const dDPM_unc &);

      static void init(int,int);

   private:

      //!nr of particles
      static int N;

      //!dimension of sp hilbert space
      static int M;

      //!double array of xTPM_unc objects
      xTPM_unc **dpphm;

};

#endif
