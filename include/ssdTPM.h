#ifndef ssdTPM_H
#define ssdTPM_H

#include <iostream>
#include <fstream>

using std::ostream;

#include "SPM.h"

/**
 * @author Brecht Verstichel
 * @date 23-06-2011\n\n
 * This class ssdTPM is an array of SPM matrices, it will contain a special "bar of the different blocks of the dDPM and consorts.
 * It is a called ssdTPM because it has no spin dependence anymore, unlike the true dTPM, so spinsummed (ss) dTPM.
 */
class ssdTPM {

   /**
    * Output stream operator overloaded
    * @param output The stream to which you are writing (e.g. cout)
    * @param dtpm_p the ssdTPM you want to print
    */
   friend ostream &operator<<(ostream &output,const ssdTPM &dtpm_p);

   public:

      //constructor
      ssdTPM();

      //copy constructor
      ssdTPM(const ssdTPM &);

      //destructor
      virtual ~ssdTPM();

      SPM &operator[](int);

      const SPM &operator[](int) const;

      int gN() const;

      int gM() const;

      void bar(double,const dDPM &);

      double trace() const;

      void dscal(double);

      ssdTPM &operator=(const ssdTPM &);

      ssdTPM &operator+=(const ssdTPM &);

      ssdTPM &operator-=(const ssdTPM &);

      ssdTPM &operator=(double);

      ssdTPM &daxpy(double,const ssdTPM &);

      double ddot(const ssdTPM &) const;

      void L_map(const ssdTPM &,const ssdTPM &);

      void sqrt(int);

      void invert();

      ssdTPM &mprod(const ssdTPM &,const ssdTPM &);

      void symmetrize();

      void fill_Random();

      static void init(int,int);

   private:

      //!nr of particles
      static int N;

      //!dimension of sp hilbert space
      static int M;

      //!double array of SPM objects
      SPM **ssdtpm;

};

#endif
