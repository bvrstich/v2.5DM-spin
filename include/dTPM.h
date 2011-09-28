#ifndef dTPM_H
#define dTPM_H

#include <iostream>
#include <fstream>

using std::ostream;

#include "SPM.h"

/**
 * @author Brecht Verstichel
 * @date 23-06-2011\n\n
 * This class dTPM is an array of SPM matrices, it will contain the bar of the different blocks of the dDPM and consorts.
 */
class dTPM {

   /**
    * Output stream operator overloaded
    * @param output The stream to which you are writing (e.g. cout)
    * @param dtpm_p the dTPM you want to print
    */
   friend ostream &operator<<(ostream &output,const dTPM &dtpm_p);

   public:

      //constructor
      dTPM();

      //copy constructor
      dTPM(const dTPM &);

      //destructor
      virtual ~dTPM();

      SPM &operator[](int);

      const SPM &operator[](int) const;

      int gN() const;

      int gM() const;

      void freebar(double,const dDPM &);

      double trace() const;

      void dscal(double);

      dTPM &operator=(const dTPM &);

      dTPM &operator+=(const dTPM &);

      dTPM &operator-=(const dTPM &);

      dTPM &operator=(double);

      dTPM &daxpy(double,const dTPM &);

      double ddot(const dTPM &) const;

      void L_map(const dTPM &,const dTPM &);

      void sqrt(int);

      void invert();

      dTPM &mprod(const dTPM &,const dTPM &);

      void symmetrize();

      void fill_Random();

      static void init(int,int);

   private:

      //!nr of particles
      static int N;

      //!dimension of sp hilbert space
      static int M;

      //!double array of SPM objects
      SPM **dtpm;

};

#endif
