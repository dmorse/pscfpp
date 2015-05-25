/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "SymmTensorAverage.h"   // class header
#include <util/space/Tensor.h>

#include <math.h>

namespace Util
{

   /*
   * Default constructor.
   */
   SymmTensorAverage::SymmTensorAverage(int blockFactor)
    : ParamComposite(),
      nSamplePerBlock_(0),
      iBlock_(0)
   {
      setClassName("SymmTensorAverage");
      int i, j, k;
      k = 0;
      for (i = 0; i < Dimension; ++i) {
         for (j = 0; j <= i; ++j) {
            accumulators_[k].setBlockFactor(blockFactor);
            ++k;
         }
      }
   }

   /*
   * Destructor.
   */
   SymmTensorAverage::~SymmTensorAverage()
   {}

   /*
   * Set nSamplePerBlock parameter.
   */
   void SymmTensorAverage::setNSamplePerBlock(int nSamplePerBlock)
   {
      if (nSamplePerBlock < 0) {
         UTIL_THROW("Attempt to set nSamplePerBlock < 0");
      }
      nSamplePerBlock_ = nSamplePerBlock;
      int i, j, k;
      k = 0;
      for (i = 0; i < Dimension; ++i) {
         for (j = 0; j <= i; ++j) {
            accumulators_[k].setNSamplePerBlock(nSamplePerBlock);
            ++k;
         }
      }
   }

   /*
   * Read nSamplePerBlock from file.
   */
   void SymmTensorAverage::readParameters(std::istream& in)
   {
      read<int>(in, "nSamplePerBlock", nSamplePerBlock_);
      if (nSamplePerBlock_ < 0) {
         UTIL_THROW("Invalid input: nSamplePerBlock < 0");
      }
      int i, j, k;
      k = 0;
      for (i = 0; i < Dimension; ++i) {
         for (j = 0; j <= i; ++j) {
            accumulators_[k].setNSamplePerBlock(nSamplePerBlock_);
            ++k;
         }
      }
   }

   /*
   * Load internal state from archive.
   */
   void SymmTensorAverage::loadParameters(Serializable::IArchive &ar)
   {
      loadParameter<int>(ar, "nSamplePerBlock", nSamplePerBlock_);
      if (nSamplePerBlock_ < 0) {
         UTIL_THROW("Loading value nSamplePerBlock < 0");
      }
      ar & iBlock_;
      int i, j, k;
      k = 0;
      for (i = 0; i < Dimension; ++i) {
         for (j = 0; j <= i; ++j) {
            ar & accumulators_[k];
            ++k;
         }
      }
   }

   /*
   * Save internal state to archive.
   */
   void SymmTensorAverage::save(Serializable::OArchive &ar)
   {  ar & *this; }

   /*
   * Reset all accumulators and counters to zero.
   */
   void SymmTensorAverage::clear()
   {
      iBlock_ = 0;
      int i, j, k;
      k = 0;
      for (i = 0; i < Dimension; ++i) {
         for (j = 0; j <= i; ++j) {
            accumulators_[k].clear();
            ++k;
         }
      }
   }

   /*
   * Add a sampled value to the ensemble.
   */
   void SymmTensorAverage::sample(const Tensor& value)
   {
      int i, j, k;
      k = 0;
      for (i = 0; i < Dimension; ++i) {
         for (j = 0; j <= i; ++j) {
            accumulators_[k].sample(value(i,j));
            ++k;
         }
      }
      if (nSamplePerBlock_) {
         if (iBlock_ == nSamplePerBlock_) {
            iBlock_ = 0;
         }
         ++iBlock_;
      }
   }

   /*
   * Access accumulator associated with one component.
   */
   const Average& SymmTensorAverage::operator () (int i, int j)
   {
      int k;
      if (j > i) {
        k = i;
        i = j;
        j = k;
      }
      k = i*(i+1)/2 + j;
      return accumulators_[k];
   }

}
