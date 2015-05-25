#ifndef UTIL_AUTOCORR_STAGE_TPP
#define UTIL_AUTOCORR_STAGE_TPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "AutoCorrStage.h"

#include <util/accumulators/setToZero.h>
#include <util/accumulators/product.h>
#include <util/format/Int.h>
#include <util/format/write.h>

#include <complex>
#include <math.h>

using std::complex;

namespace Util
{

   /*
   * Constructor for rootPtr AutoCorrStage, with stageId = 0.
   */
   template <typename Data, typename Product>
   AutoCorrStage<Data, Product>::AutoCorrStage()
    : bufferCapacity_(64),
      maxStageId_(10),
      blockFactor_(2),
      buffer_(),
      corr_(),
      nCorr_(),
      nSample_(0),
      blockSum_(),
      nBlockSample_(0),
      stageInterval_(1),
      childPtr_(0),
      rootPtr_(0),
      stageId_(0)
   {
      rootPtr_ = this;
      setToZero(blockSum_);
   }

   /*
   * Private constructor for stages with stageId > 0.
   */
   template <typename Data, typename Product>
   AutoCorrStage<Data, Product>::AutoCorrStage(long stageInterval,
                                               int stageId, int maxStageId,
                                               AutoCorrStage<Data, Product>* rootPtr,
                                               int blockFactor)
    : bufferCapacity_(rootPtr->bufferCapacity()),
      maxStageId_(maxStageId),
      blockFactor_(blockFactor),
      buffer_(),
      corr_(),
      nCorr_(),
      nSample_(0),
      blockSum_(),
      nBlockSample_(0),
      stageInterval_(stageInterval),
      childPtr_(0),
      rootPtr_(rootPtr),
      stageId_(stageId)
   {
      allocate();
      setToZero(blockSum_);
   }

   /*
   * Destructor.
   */
   template <typename Data, typename Product>
   AutoCorrStage<Data, Product>::~AutoCorrStage()
   {
      if (childPtr_) {
         delete childPtr_;
      }
   }

   /*
   * Set buffer capacity and initialize.
   */
   template <typename Data, typename Product>
   void AutoCorrStage<Data, Product>::setParam(int bufferCapacity,
                                               int maxStageId, int blockFactor)
   {
      bufferCapacity_ = bufferCapacity;
      maxStageId_ = maxStageId;
      blockFactor_ = blockFactor;
      allocate();
   }

   /*
   * Set previously allocated to initial empty state.
   */
   template <typename Data, typename Product>
   void AutoCorrStage<Data, Product>::clear()
   {
      setToZero(blockSum_);
      nSample_ = 0;
      nBlockSample_ = 0;
      if (bufferCapacity_ > 0) {
         for (int i=0; i < bufferCapacity_; ++i) {
            setToZero(corr_[i]);
            nCorr_[i] = 0;
         }
         buffer_.clear();
      }
   }

   /*
   * Add a sampled value to the ensemble.
   */
   template <typename Data, typename Product>
   void AutoCorrStage<Data, Product>::sample(Data value)
   {
      if (bufferCapacity_ <= 0) {
         bufferCapacity_ = 64;
         blockFactor_ = 2;
         maxStageId_ = 10;
         allocate();
      }

      // Increment global accumulators
      buffer_.append(value);
      for (int i=0; i < buffer_.size(); ++i) {
         corr_[i] += product(buffer_[i], value);
         ++nCorr_[i];
      };
      ++nSample_;

      // Increment block accumulators
      blockSum_ += value;
      ++nBlockSample_;

      if (nBlockSample_ == blockFactor_) {
         if (stageId_ < maxStageId_) {
            if (!childPtr_) {
               long nextStageInterval = stageInterval_*blockFactor_;
               int  nextStageId = stageId_ + 1;
               childPtr_ = new AutoCorrStage(nextStageInterval, nextStageId,
                                             maxStageId_, rootPtr_, blockFactor_);
               rootPtr_->registerDescendant(childPtr_);
            }
            // Add block average as first value in child sequence
            blockSum_ /= double(blockFactor_);
            childPtr_->sample(blockSum_);
         }
         // Reset block accumulators
         setToZero(blockSum_);
         nBlockSample_ = 0;
      }

   }

   /*
   * Serialize this AutoCorr.
   */
   template <typename Data, typename Product>
   template <class Archive>
   void
   AutoCorrStage<Data, Product>::serialize(Archive& ar,
                                           const unsigned int version)
   {
      // Parameters (required for initialization)
      ar & bufferCapacity_;
      ar & maxStageId_;
      ar & blockFactor_;

      // Serialize private state info.
      serializePrivate(ar, version);
   }

   /*
   * Serialize private info for this AutoCorrStage.
   */
   template <typename Data, typename Product>
   template <class Archive>
   void
   AutoCorrStage<Data, Product>::serializePrivate(Archive& ar,
                                          const unsigned int version)
   {
      // State of simple nonblock correlator
      ar & buffer_;
      ar & corr_;
      ar & nCorr_;
      ar & nSample_;
      isValid();

      ar & blockSum_;
      ar & nBlockSample_;

      // Constructor sets rootPtr_, stageInterval_ and stageId_

      // Does this stage have a child?
      int hasChild;
      if (Archive::is_saving()) {
         hasChild = (childPtr_ == 0) ? 0 : 1;
      }
      ar & hasChild;

      // Serialize child (if any)
      if (hasChild) {
         if (Archive::is_loading()) {
            long nextStageInterval = stageInterval_*blockFactor_;
            int  nextStageId       = stageId_ + 1;
            childPtr_ = new AutoCorrStage(nextStageInterval,
                                          nextStageId, maxStageId_,
                                          rootPtr_, blockFactor_);
            rootPtr_->registerDescendant(childPtr_);
         }
         ar & (*childPtr_);
      } else {
         if (Archive::is_loading()) {
            childPtr_ = 0;
         }
      }

   }

   /*
   * Return capacity of history buffer.
   */
   template <typename Data, typename Product>
   int AutoCorrStage<Data, Product>::bufferCapacity() const
   {  return bufferCapacity_;  }

   /*
   * Return current size of the history buffer.
   */
   template <typename Data, typename Product>
   int AutoCorrStage<Data, Product>::bufferSize() const
   {  return buffer_.size(); }

   /*
   * Return the number of sampled values.
   */
   template <typename Data, typename Product>
   long AutoCorrStage<Data, Product>::nSample() const
   {  return nSample_; }

   /*
   * Return the number of measured values per sample at this stage.
   */
   template <typename Data, typename Product>
   long AutoCorrStage<Data, Product>::stageInterval() const
   {  return stageInterval_; }

   /*
   * Calculate and output autocorrelation function, assuming zero average.
   */
   template <typename Data, typename Product>
   void AutoCorrStage<Data, Product>::output(std::ostream& outFile)
   {
      Product aveSq;
      setToZero(aveSq);
      output(outFile, aveSq);
   }

   /*
   * Calculate and output autocorrelation function.
   */
   template <typename Data, typename Product>
   void AutoCorrStage<Data, Product>::output(std::ostream& outFile, Product aveSq)
   {
      int min;
      if (stageId_ == 0) {
         min = 0;
      } else {
         min = bufferCapacity_ / blockFactor_;
      }

      Product autocorr;
      for (int i = min; i < buffer_.size(); ++i) {
         autocorr = corr_[i]/double(nCorr_[i]);
         autocorr -= aveSq;
         outFile << Int(i*stageInterval_) << " ";
         write<Product>(outFile, autocorr);
         outFile << std::endl;
      }
      if (childPtr_) {
         childPtr_->output(outFile, aveSq);
      }
   }

   /*
   * Return autocorrelation at a given lag time, assuming zero average.
   */
   template <typename Data, typename Product>
   Product AutoCorrStage<Data, Product>::autoCorrelation(int t) const
   {
      Product aveSq;
      setToZero(aveSq);
      return autoCorrelation(t, aveSq);
   }

   /*
   * Return autocorrelation at a given lag time.
   */
   template <typename Data, typename Product>
   Product AutoCorrStage<Data, Product>::autoCorrelation(int t, Product aveSq) const
   {
      assert(t < buffer_.size());
      Product autocorr = corr_[t]/double(nCorr_[t]);
      autocorr -= aveSq;
      return autocorr;
   }

   /*
   *  Return correlation time, in Data samples, assuming zero average.
   */
   template <typename Data, typename Product>
   double AutoCorrStage<Data, Product>::corrTime() const
   {
      Product aveSq;
      setToZero(aveSq);
      return corrTime(aveSq);
   }

   /*
   *  Return correlation time in unit of sampling interval.
   */
   template <typename Data, typename Product>
   double AutoCorrStage<Data, Product>::corrTime(Product aveSq) const
   {
      // Compute variance from value C(t=0)
      Product variance = corr_[0];
      variance /= double(nCorr_[0]);
      variance -= aveSq;

      // Compute integral of C(t)/C(0)
      Product autocorr, sum;
      setToZero(sum);
      int  size = buffer_.size();
      for (int i = 1; i < size/2; ++i) {
         autocorr = corr_[i];
         autocorr /= double(nCorr_[i]);
         autocorr -= aveSq;
         sum += autocorr;
      }
      sum /= variance;
      return sum;
   }

   // Private member functions

   /*
   * Allocate arrays and CyclicBuffer, and initialize.
   */
   template <typename Data, typename Product>
   void AutoCorrStage<Data, Product>::allocate()
   {
      if (bufferCapacity_ > 0) {
         corr_.allocate(bufferCapacity_);
         nCorr_.allocate(bufferCapacity_);
         buffer_.allocate(bufferCapacity_);
      }
      AutoCorrStage<Data, Product>::clear();
   }

   /*
   * Are capacities consistent?
   */
   template <typename Data, typename Product>
   bool AutoCorrStage<Data, Product>::isValid()
   {
      bool valid = true;
      if (bufferCapacity_ != corr_.capacity()) valid = false;
      if (bufferCapacity_ != nCorr_.capacity()) valid = false;
      if (bufferCapacity_ != buffer_.capacity()) valid = false;
      if (!valid) {
         UTIL_THROW("Invalid AutoCorrStage");
      }
      return valid;
   }

   template <typename Data, typename Product>
   void
   AutoCorrStage<Data, Product>::registerDescendant(AutoCorrStage<Data, Product>* ptr)
   {}

}
#endif
