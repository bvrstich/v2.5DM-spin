#ifndef dDPV_H
#define dDPV_H

#include <iostream>
#include <fstream>

using std::ostream;

#include "Vector.h"

/**
 * @author Brecht Verstichel
 * @date 30-05-2011\n\n
 * This class dDPV is a class written the eigenvalues of the dDPM object.
 */
class dDPV {

   /**
    * Output stream operator overloaded
    * @param output The stream to which you are writing (e.g. cout)
    * @param ddpv_p the dDPV you want to print
    */
   friend ostream &operator<<(ostream &output,const dDPV &ddpv_p);

   public:

      //standard constructor
      dDPV(dDPM &);
      
      //copy constructor
      dDPV(const dDPV &);

      //destructor
      virtual ~dDPV();

      //geef N terug
      int gN() const;

      //geef M terug
      int gM() const;

      BlockVector<rTPM> &operator[](int);

      const BlockVector<rTPM> &operator[](int) const;

      dDPV &operator=(const dDPV &);

      dDPV &operator=(double );

      dDPV &operator+=(const dDPV &);

      dDPV &operator-=(const dDPV &);

      dDPV &daxpy(double alpha,const dDPV &);

      double sum() const;

      double log_product() const;

      double ddot(const dDPV &) const;

      void dscal(double alpha);

      double min() const;
      
      double max() const;

      double lsfunc(double) const;

      static void init(int,int);

   private:
      
      //!double pointer to vector<rTPM> that contain all the eigenvalues of a dDPM object
      BlockVector<rTPM> **ddpv;

      //!nr of particles
      static int N;

      //!dimension of sp hilbert space
      static int M;

};

#endif
