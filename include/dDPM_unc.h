#ifndef dDPM_unc_H
#define dDPM_unc_H 
#include <iostream>
#include <fstream>

using std::ostream;

#include "rxTPM_unc.h"

/**
 * @author Brecht Verstichel
 * @date 13-05-2011\n\n
 * This class dDPM_unc is a class written for the I1 and Q2 conditions, which is a DPM matrix with the third indices diagonal.
 * The v2.5DM program will use this as the central variable. 
 */
class dDPM_unc {

   /**
    * Output stream operator overloaded
    * @param output The stream to which you are writing (e.g. cout)
    * @param ddpm_p the dDPM_unc you want to print
    */
   friend ostream &operator<<(ostream &output,const dDPM_unc &ddpm_p);

   public:

      //constructor
      dDPM_unc();

      //copy constructor
      dDPM_unc(const dDPM_unc &);

      //destructor
      virtual ~dDPM_unc();

      rxTPM_unc &operator[](int);

      const rxTPM_unc &operator[](int) const;

      int gN() const;

      int gM() const;

      double trace() const;

      void dscal(double);

      dDPM_unc &operator=(const dDPM_unc &);

      dDPM_unc &operator+=(const dDPM_unc &);

      dDPM_unc &operator-=(const dDPM_unc &);

      dDPM_unc &operator=(double);

      dDPM_unc &daxpy(double,const dDPM_unc &);

      double ddot(const dDPM_unc &) const;

      void L_map(const dDPM_unc &,const dDPM_unc &);

      void sqrt(int);

      void invert();

      dDPM_unc &mprod(const dDPM_unc &,const dDPM_unc &);

      void symmetrize();

      void fill_Random();

      void print_eig();

      void uncouple(const dDPM &);

      static void init(int,int);

   private:

      //!nr of particles
      static int N;

      //!dimension of sp hilbert space
      static int M;

      //!double array of rxTPM_unc objects
      rxTPM_unc **ddpm;

};

#endif
