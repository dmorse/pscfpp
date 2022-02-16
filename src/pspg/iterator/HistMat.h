#ifndef PSPG_HIST_MAT_H
#define PSPG_HIST_MAT_H

#include <util/containers/DMatrix.h>

namespace Pscf {
namespace Pspg
{
   using namespace Util;

   //upper triangular
   class HistMat : public DMatrix<double> 
   {
   public:
      HistMat();

      HistMat(int n);

      ~HistMat();

      //allocate an internal matrix of size maxHist + 1 to keep d(k)d'(k) values
      //requires maxHist to be constant within simulation.
      void allocate(int capacity);

      //reset the state of the HistMat
      //the matrix remains allocated with all values set to zero
      //boolean flags all reset to initial state
      void reset();

      void clearColumn(int nHist);

      //insert new d(k)d'(k) value
      void evaluate(float elm, int nHist, int columnId);

      double makeUmn(int i,int j, int nHist);

      double makeVm(int i, int nHist);

   private:
      int maxHist_; //keeps a copy of max history
      bool filled_ = false;
   };

   inline HistMat::HistMat() 
      : maxHist_(0){}

   inline HistMat::HistMat(int n)
    : maxHist_(n){}

   inline HistMat::~HistMat() {}

   inline void HistMat::allocate(int capacity) {
      maxHist_ = capacity - 1;
      DMatrix<double>::allocate(capacity, capacity);
   }

   inline void HistMat::reset() {
      filled_ = false;
      for (int i = 0; i < maxHist_ + 1; ++i) {
         for (int j = 0; j < maxHist_ + 1; ++j) {
            (*this)(i, j) = 0;
         }
      }
   }

   inline void HistMat::clearColumn(int nHist) {
      //if the matrix is not entirely filled
      if (nHist <= maxHist_ && !filled_) {
         //set to zero
         for (int i = maxHist_; i > maxHist_ - nHist - 1; --i) {
            (*this)(maxHist_ - nHist, i) = 0;
         }
         if (nHist == maxHist_) {
            filled_ = true;
         }
      }
      else { //out of space
         for (int i = maxHist_; i > 0; i--) {
            for (int j = maxHist_; j > i - 1; j--) {
               (*this)(i, j) = (*this)(i - 1, j - 1);
            }
         }
         for (int i = 0; i < maxHist_ + 1; i++) {
            (*this)(0, i) = 0;
         }
      }
   }

   inline void HistMat::evaluate(float elm, int nHist, int columnId) {
      if (nHist < maxHist_) { //if matrix is not filled, something tricky
         (*this)(maxHist_ - nHist, columnId + maxHist_ - nHist) = (double)elm;
      }
      else { //just fill the first row
         (*this)(0, columnId) = (double)elm;
      }
   }

   inline double HistMat::makeUmn(int i, int j, int nHist) {
      int offset;
      if (nHist < maxHist_) { //matrix is not filled, soemthing tricky
         offset = maxHist_ - nHist;
      }
      else { //matrix is filled
         offset = 0;
      }
      return (*this)(offset, offset) + (*this)(offset + 1 + i , offset + 1 + j) -
         (*this)(offset, offset + 1 + i) - (*this)(offset, offset + 1 + j);
   }

   inline double HistMat::makeVm(int i, int nHist) {
      int offset;
      if (nHist < maxHist_) {
         offset = maxHist_ - nHist;
      }
      else {
         offset = 0;
      }
      return (*this)(offset, offset) - (*this)(offset, offset + i + 1);
   }

}
}
#endif
