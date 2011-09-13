#ifndef SUP_H
#define SUP_H

#include <iostream> 
#include <fstream> 

using std::ostream;

/**
 * @author Brecht Verstichel
 * @date 09-06-2011\n\n
 * This class, SUP is a blockmatrix over the carrierspace's of active N-representability conditions. 
 * This class contains two dDPM objects containing the I1 and Q2 conditions, two dPPHM objects containing I2 and Q1
 * and two dPHHM objects containing G1 and G2.
 */
class SUP{
  
   /**
    * Output stream operator overloaded
    * @param output The stream to which you are writing (e.g. cout)
    * @param X_p the SUP you want to print
    */
   friend ostream &operator<<(ostream &output,const SUP &X_p);

   public:

      //constructor
      SUP();

      //copy constructor
      SUP(const SUP &);

      //destructor
      ~SUP();

      //overload += operator
      SUP &operator+=(const SUP &);

      //overload -= operator
      SUP &operator-=(const SUP &);

      //overload equality operator
      SUP &operator=(const SUP &);

      //overload equality operator
      SUP &operator=(double &);

      dDPM &gI1();

      const dDPM &gI1() const;

      int gN() const;

      int gM() const;

      int gdim() const;

      double ddot(const SUP &) const;

      void invert();

      void dscal(double alpha);

      //positieve of negatieve vierkantswortel uit een supermatrix
      void sqrt(int option);

      void L_map(const SUP &,const SUP &);

      void daxpy(double alpha,const SUP &);

      double trace() const;

      SUP &mprod(const SUP &,const SUP &);

      void fill(const dDPM &);

      void fill_Random();

      static void init(int,int);

   private:

      //!pointer to a dDPM object, containts the W matrix that has to be positive semidefinite: I1 condition
      dDPM *I1;

      //!number of spatial orbs
      static int M;

      //!nr of particles
      static int N;

      //!total dimension of the SUP matrix
      static int dim;

};

#endif
