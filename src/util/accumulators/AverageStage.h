#ifndef UTIL_AVERAGE_STAGE_H
#define UTIL_AVERAGE_STAGE_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/global.h>

namespace Util
{

   /**
   * Evaluate average with hierarchical blocking error analysis.
   *
   * This class implements an algorithm to evaluate the average of a 
   * sequence, using a hierarchical blocking algorithm to estimate the
   * error on the average. The algorithm is based on the calculation of
   * variances for sequences of block averages for multiple levels of
   * block sizes, as described in the following reference:
   *
   * ``Error estimates on averages of correlated data", H. Flyvbjerg 
   *  and H.G. Petersen, J. Chem. Phys. 91, pgs. 461-466 (1989).
   *
   * The blocking algorithm is implemented here by a creating a linked 
   * list of AverageStage objects, each of which is responsible for
   * computing the variance on block averages using a different level
   * of blocking. Each object in this list is assigned an integer chainId.
   * The first AverageStage object in the list, with chainId=0, calculates 
   * the average and variance for a "primary" sequence of measured values 
   * that are passed as parameters to its sample method. This first object
   * is normally an instance of the Average class, which is a subclass of
   * AverageStage that implements features that are only required by the 
   * primary stage. This object has a pointer to a child AverageStage with 
   * chainId=1 that calculates the variance of a secondary sequence in which 
   * each value is the average of blockFactor consecutive values in the 
   * primary sequence. The object with chainId=1 in turn has has a pointer 
   * to a child object with chainId=2 that calculates the variance of a 
   * sequence in which each value is the average of a block of blockFactor**2 
   * consecutive values of the primary sequence, and so on. In general, the 
   * object with chainId=n, calculates the variance of a sequence in which 
   * each value is an average of blockFactor**n values of the primary 
   * sequence. Each value in the sequence analyzed by the object with 
   * chainId=n+1 is calculated by the parent object with chainId=n, by 
   * calculating an average of a block of blockFactor consecutive values
   * of its own sequence and passing this block average as a parameter
   * the sample() function of the object with chainId=n+1. New stages in
   * this linked list are instantiated and to the list as needed as the 
   * length of the primary sequence grows: When an object with chainId=n 
   * has been passed a sequence of exactly blockFactor values, it creates 
   * a child AverageStage object with chainId=n+1 and passes the average 
   * of these first blockFactor values to the sample function of the 
   * child object as the first value in its sequence.
   *
   * A value of the integer parameter blockFactor is passed to the 
   * constructor of the primary AverageStage object. This parameter is 
   * set to blockFactor=2 by default. Its value may be reset using the 
   * setBlockFactor() function before any data is sampled, but may not
   * be changed thereafter.
   *
   * \ingroup Accumulators_Module
   */
   class AverageStage
   {

   public:

      /**
      * Constructor
      *
      * This constructor creates a primary AverageStage object with
      * stageId = 0 and stageInterval = 1. A private constructor is
      * used to recursively create children of this object.
      *
      * \param blockFactor ratio of block sizes of subsequent stages
      */
      AverageStage(int blockFactor = 2);

      /**
      * Destructor.
      *
      * Recursively destroy all children.
      */
      virtual ~AverageStage();

      /**
      * Reset the value of blockFactor.
      *
      * \throw Exception if called when nSample > 0.
      */
      void setBlockFactor(int blockFactor);

      /**
      * Initialize all accumulators and recursively destroy all children.
      */
      virtual void clear();

      /**
      * Add a sampled value to the ensemble.
      *
      * \param value sampled value
      */
      virtual void sample(double value);

      /**
      * Add a sampled value to the ensemble.
      *
      * \param ar       input or output archive
      * \param version  file version id
      */
      template <class Archive>
      void serialize(Archive& ar, const unsigned int version);

      ///\name Accessors
      //@{

      /**
      * Return the average of all sampled values.
      */
      double average() const;

      /**
      * Return the variance of all sampled values.
      */
      double variance() const;

      /**
      * Return the standard deviation of all sampled values.
      *
      * \return sqrt(variance())
      */
      double stdDeviation() const;

      /**
      * Return a naive estimate for the std deviation of the average.
      *
      * \return sqrt(variance()/nSample())
      */
      double error() const;

      /**
      * Return the number of sampled values in this sequence.
      */
      long nSample() const;

      /**
      * Return the number of sampled values per block at this stage.
      */
      long stageInterval() const;

      //@}

   protected:

      /**
      * Does this object have a child AverageStage for block averages?
      */
      bool hasChild() const;

      /**
      * Return the child AverageStage by reference.
      */
      AverageStage& child();

   private:

      /// Sum of all sampled values.
      double sum_;

      /// Sum of squares of all sampled values.
      double sumSq_;

      /// Sum of sampled values in the current block.
      double blockSum_;

      /// Number of sampled values.
      long nSample_;

      /// Number of values in the current block.
      long nBlockSample_;

      /// Number of measured values per sampled value at this stage.
      long stageInterval_;

      /// Pointer to child stage, if any.
      AverageStage* childPtr_;

      /// Pointer to root stage. Null if this is the root stage.
      AverageStage* rootPtr_;

      /// Stage index
      int stageId_;

      /// Number of samples per block.
      int blockFactor_;

      /**
      * Constructor for child objects (private).
      *
      * \param stageInterval  number of values per sample at this stage
      * \param stageId  integer id for this stage
      * \param rootPtr  pointer to root AverageStage
      * \param blockFactor  ratio of block sizes of subsequent stages
      */
      AverageStage(long stageInterval, int stageId, 
                   AverageStage* rootPtr, int blockFactor);

      /**
      * Copy constructor (private and not implemented).
      */
      AverageStage(const AverageStage& other);

      /**
      * Assignment (private and not implemented).
      */
      AverageStage& operator = (const AverageStage& other);

      /**
      * Register the creation of a descendant stage.
      *
      * This should be called only by a root stage.
      *
      * \param descendantPtr pointer to a descendant AverageStage.
      */
      virtual void registerDescendant(AverageStage* descendantPtr);

   };

   // Inline methods

   /*
   * Does this object have a child?  (protected)
   */
   inline bool AverageStage::hasChild() const
   { return bool(childPtr_); }

   /*
   * Return child object by reference. (protected)
   */
   inline AverageStage& AverageStage::child()
   { return *childPtr_; }

   // Method template

   /*
   * Serialize this stage.
   */
   template <class Archive>
   void AverageStage::serialize(Archive& ar, const unsigned int version)
   {
      ar & sum_;
      ar & sumSq_;
      ar & blockSum_;
      ar & nSample_;
      ar & nBlockSample_;
      ar & stageInterval_;
      ar & blockFactor_;

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
            int  nextStageId = stageId_ + 1;
            childPtr_ = new AverageStage(nextStageInterval, nextStageId,
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

}
#endif
