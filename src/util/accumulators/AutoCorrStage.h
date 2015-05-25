#ifndef UTIL_AUTO_CORR_STAGE_H
#define UTIL_AUTO_CORR_STAGE_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/containers/RingBuffer.h>   // member
#include <util/containers/DArray.h>       // member
#include <util/global.h>

namespace Util
{

   /**
   * Hierarchical auto-correlation function algorithm.
   *
   * This class calculates an autocorrelation function for a sequence
   * x(i) of values of a variable or object of type Data. The resulting
   * autocorrelation function is and array of values of type Product,
   * where C(j) = <x(i-j), x(i)>. Here <A,B> denotes an inner product 
   * of type Product for objects A and B of type Data.
   *
   * The meaning of the inner product is defined for various data 
   * types b the overloaded function Product product(Data, Data) that 
   * is defined for double, complex and Vector data in the product.h 
   * file.
   *
   * The zero value for variables of type Data is returned by the 
   * overloaded function void setToZero(Data) method defined in the 
   * setToData.h file.
   *
   * This class implements a hierarchical algorithm to calculate C(j).
   * The algorithm is implemented by a linked list of AutoCorrStage 
   * objects. Each object in this list is assigned an integer chainId.  
   * The "primary" AutoCorrStage object in this list, with chainId=0, 
   * calculates the autocorrelation for a primary sequence of primary
   * Data values that are passed to the sample method of this object. 
   * For each n > 0, the object with chainId = n calculates the 
   * autocorrelation function for a sequence of values in which each
   * value is an average of a block of blockFactor**n consecutive 
   * values of the primary sequence or, equivalently, an average of 
   * blockFactor consecutive values of the sequence maintained by the 
   * parent object with chainId = n-1. Additional stages are added to 
   * this list dynamically as needed.
   *
   * \ingroup Accumulators_Module
   */
   template <typename Data, typename Product>
   class AutoCorrStage 
   {

   public:

      /**
      * Constructor
      *
      * This constructor creates a primary AutoCorrStage object with
      * stageId = 0 and stageInterval = 1. A private constructor is
      * used to recursively create descendant stages as needed.
      */
      AutoCorrStage();

      /**
      * Destructor.
      *
      * Recursively destroy all descendant stages.
      */
      virtual ~AutoCorrStage();

      /**
      * Set all parameters and allocate to initialize state.
      *
      * \param bufferCapacity max. number of values stored in buffer
      * \param maxStageId maximum stage index (0=primary)
      * \param blockFactor ratio of block sizes of subsequent stages
      */
      void setParam(int bufferCapacity=64, int maxStageId=0, 
                    int blockFactor=2);

      /**
      * Sample a value.
      *
      * \param value current Data value
      */
      virtual void sample(Data value);

      /**
      * Clear accumulators and destroy descendants.
      */
      void clear();

      /**
      * Serialize to/from an archive.
      *
      * \param ar      archive
      * \param version archive version id
      */
      template <class Archive>
      void serialize(Archive& ar, const unsigned int version);

      //@}
      ///\name Accessors
      //@{

      /**
      * Output the autocorrelation function, assuming zero mean.
      *
      * This calls output(std::ostream out, Product aveSq) with 
      * a zero value for aveSq.
      *
      * \param out output stream.
      */
      void output(std::ostream& out);

      /**
      * Output the autocorrelation function.
      *
      * \param out output stream.
      * \param aveSq square of average <x>^2 subtracted from <x(t)x(0)>
      */
      void output(std::ostream& out, Product aveSq);

      /**
      * Return capacity of history buffer.
      */
      int bufferCapacity() const;

      /**
      * Return current size of history buffer.
      */
      int bufferSize() const;

      /**
      * Return the number of sampled values.
      */
      long nSample() const;

      /**
      * Return the number of primary values per block at this stage.
      */
      long stageInterval() const;

      /**
      * Return autocorrelation at a given time, assuming zero average.
      *
      * This calls autoCorrelations(t, aveSq) with a zero value for
      * for aveSq.
      *
      * \param t the lag time, in Data samples
      */
      Product autoCorrelation(int t) const;

      /**
      * Return autocorrelation at a given lag time.
      *
      * \param t the lag time, in Data samples
      * \param aveSq square average <x>^2 subtracted from <x(t)x(0)>
      */
      Product autoCorrelation(int t, Product aveSq) const;

      /**
      * Estimate of autocorrelation time, in samples.
      *
      * This variant assumes a zero average. 
      */
      double corrTime() const;

      /**
      * Numerical integration of autocorrelation function.
      *
      * \param aveSq square average <x>^2 subtracted from <x(t)x(0)>
      */
      double corrTime(Product aveSq) const;

      //@}

   protected:

      // Physical capacity (# of elements) of buffer, corr, and nCorr.
      int bufferCapacity_;

      /// Maximum allowed stage index (controls maximum degree of blocking).
      int maxStageId_;

      /// Number of values per block (ratio of intervals for successive stages).
      int blockFactor_;

      /**
      * Allocate memory and initialize to empty state.
      */
      void allocate();

      /**
      * Does this have a child AutoCorrStage?
      */
      bool hasChild() const;

      /**
      * Return the child AutoCorrStage by reference.
      */
      AutoCorrStage& child();

      /**
      * Register the creation of a descendant stage.
      *
      * This should be called only by a root stage.
      *
      * \param ptr pointer to a descendant AutoCorrStage.
      */
      virtual void registerDescendant(AutoCorrStage<Data, Product>* ptr);

      /**
      * Serialize private data members, and descendants.
      *
      * \param ar      archive
      * \param version archive version id
      */
      template <class Archive>
      void serializePrivate(Archive& ar, const unsigned int version);

   private:

      // Ring buffer containing sequence of Data values
      RingBuffer<Data> buffer_;

      // Array in which corr_[j] = sum of values of <x(i-j), x(i)>
      DArray<Product> corr_;

      // Array in which nCorr_[i] = number of products added to corr_[i]
      DArray<int> nCorr_;

      // Sum of all previous values of x(t)
      Data sum_;

      /// Number of sampled values.
      long nSample_;

      /// Sum of sampled values in the current block.
      Data blockSum_;

      /// Number of values in the current block.
      long nBlockSample_;

      /// Number of measured values per sampled value at this stage.
      long stageInterval_;

      /// Pointer to child stage, if any.
      AutoCorrStage<Data, Product>* childPtr_;

      /// Pointer to root stage.
      AutoCorrStage<Data, Product>* rootPtr_;

      /// Stage index
      int stageId_;

      /**
      * Constructor for child objects (private).
      *
      * \param stageInterval  number of primary values per sample
      * \param stageId integer id for this stage
      * \param maxStageId integer id for this stage
      * \param rootPtr pointer to root AutoCorrStage
      * \param maxStageId index of the highest stage (primary=0)
      * \param blockFactor ratio of block sizes of subsequent stages
      */
      AutoCorrStage(long stageInterval, int stageId, int maxStageId,
                    AutoCorrStage<Data, Product>* rootPtr, int blockFactor);

      /**
      * Copy constructor - private and not implemented to prohibit copying.
      */
      AutoCorrStage(const AutoCorrStage& other);

      /**
      * Assignment - private and not implemented to prohibit assignment.
      */
      AutoCorrStage& operator = (const AutoCorrStage& other);

      /**
      * Is the internal state valid?
      */
      bool isValid();

   };

   // Inline methods

   /*
   * Does this object have a child?  (protected)
   */
   template <typename Data, typename Product>
   inline bool AutoCorrStage<Data, Product>::hasChild() const
   {  return bool(childPtr_); }

   /*
   * Return child object by reference. (protected)
   */
   template <typename Data, typename Product>
   inline AutoCorrStage<Data, Product>& AutoCorrStage<Data, Product>::child()
   { return *childPtr_; }

}
#endif
