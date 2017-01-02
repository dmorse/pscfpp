#ifndef UTIL_AVERAGE_H
#define UTIL_AVERAGE_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/global.h>
#include <util/accumulators/AverageStage.h>  // base class
#include <util/param/ParamComposite.h>       // base class

#include <vector>

namespace Util
{

   /**
   * Calculates the average and variance of a sampled property.
   *
   * Average calculates block and global averages of a sampled value 
   * and its square, from which it obtains a global average and 
   * variance for a sequence. A hierarchical blocking algorithm is 
   * used to estimate the error on the average. No error estimate 
   * is provided for the variance.
   *
   * The sample function of also optionally calculates block averages,
   * which can be useful for reducing how frequently values are logged
   * to a file. The parameter nSamplePerBlock is the number of samples 
   * per block average. This is initialized to zero. A zero value 
   * disables calculation of block averages. An overloaded method of 
   * the sample function that takes an std::ostream file as an argument 
   * outputs block averages to file as blocks are completed.
   *
   * The hierarchical blocking algorithm is implemented using a linked
   * list of Util::AverageStage objects. See documentation of that class
   * for further details, and a literature reference.
   *
   * \ingroup Accumulators_Module
   */
   class Average : public AverageStage, public ParamComposite
   {

   public:

      /**
      * Constructor
      *
      * \param blockFactor ratio of block sizes for subsequent stages.
      */
      Average(int blockFactor = 2);

      /**
      * Destructor
      */
      virtual ~Average();

      /**
      * Read parameter nSamplePerBlock from file and initialize.
      *
      * See setNSamplePerBlock() for discussion of value.
      *
      * \param in input stream
      */
      void readParameters(std::istream& in);

      /**
      * Set nSamplePerBlock.
      *
      * If nSamplePerBlock > 0, the sample function will increment block
      * averages, and reset the average every nSamplePerBlock samples.
      *
      * If nSamplePerBlock == 0, block averaging is disabled. This is the
      * default (i.e., the initial value set in the constructor).
      *
      * \param nSamplePerBlock number of samples per block average output
      */
      void setNSamplePerBlock(int nSamplePerBlock);

      /**
      * Load internal state from an archive.
      *
      * \param ar input/loading archive
      */
      virtual void loadParameters(Serializable::IArchive &ar);

      /**
      * Save internal state to an archive.
      *
      * \param ar output/saving archive
      */
      virtual void save(Serializable::OArchive &ar);

      /**
      * Serialize this Average to or from an archive.
      *
      * \param ar       input or output archive
      * \param version  file version id
      */
      template <class Archive>
      void serialize(Archive& ar, const unsigned int version);

      /**
      * Clear all accumulators, set to empty initial state.
      */
      void clear();

      /**
      * Add a sampled value to the ensemble.
      *
      * \param value sampled value
      */
      void sample(double value);

      /**
      * Add a sampled value to ensemble, and output block averages.
      *
      * \param value  sampled value
      * \param out  output stream to which to write block averages
      */
      void sample(double value, std::ostream& out);

      /**
      * Output final statistical properties to file.
      *
      * This function outputs the average value, an estimate of the
      * error on the average, the variance. It also outputs a sequence
      * of naive values for the error on the average obtained from 
      * sequences of block averages, with different levels of blocking.
      * The naive estimate obtained from each stage is calculated as if
      * subsequent values were uncorrelated. This gives sqrt(variance/nSample),
      * where variance is the variance of the sequence of block averages
      * processed by that stage, and nSample is the number of such block
      * averages thus far. The final estimate of the error on the average
      * is obtained by trying to identify several stages of block
      * averaging that yield statistically indistinguishable naive 
      * estimates.
      *
      * \param out output stream
      */
      void output(std::ostream& out) const;

      /**
      * Return estimated error on average from blocking analysis.
      */
      double blockingError() const;

      /**
      * Get number of samples per block average.
      *
      * A zero value indicates that block averaging is disabled.
      */
      int nSamplePerBlock() const;

      /**
      * Get number of samples in current block average.
      *
      * Return 0 if block averaging disabled, if !nSamplePerBlock.
      */
      int iBlock() const;

      /**
      * Is the current block average complete?
      * 
      * \return (iBlock > 0) && (iBlock == nSamplePerBlock)
      */
      bool isBlockComplete() const;
   
      /*
      * Get current block average value.
      *
      * \throw Exception if block is empty, or blocking disabled.
      */
      double blockAverage() const;

   private:

      /**
      * Add pointer to a descendant to an array.
      *
      *\param descendantPtr Pointer to descendant AverageStage.
      */
      virtual void registerDescendant(AverageStage* descendantPtr);

      /// Array of pointers to descendants.
      std::vector<AverageStage*> descendants_;

      /// Sum of values in current output block
      double blockSum_;

      /// Number of samples in current output block.
      int iBlock_;

      /// Number of sampled values per output block.
      int nSamplePerBlock_;

      /// Private and not implemented to prohibit copying.
      Average(const Average& other);

      /// Private and not implemented to prohibit assignment.
      Average& operator = (const Average& other);

   };

   // Inline method definitions

   #if 0
   /*
   * Add a sampled value to the ensemble, and output block averages.
   */
   inline void Average::sample(double value, std::ostream& out)
   {  sample(value, &out); }
   #endif

   /*
   * Get nSamplePerBlock, number of samples per block average.
   */
   inline int Average::nSamplePerBlock() const
   {  return nSamplePerBlock_; }

   /*
   * Get iBlock, number of samples in current block average.
   */
   inline int Average::iBlock() const
   {  return iBlock_; }

   /*
   * Is the current block average complete?
   */
   inline bool Average::isBlockComplete() const
   { 
      return (iBlock_ && (iBlock_ == nSamplePerBlock_));
   }

   /*
   * Get current block average.
   */
   inline double Average::blockAverage() const
   {
      if (iBlock_ == 0) {
         UTIL_THROW("Attempt to get block average with no data");
      }
      return blockSum_/double(iBlock_);
   }

   /*
   * Serialize this Average.
   */
   template <class Archive>
   void Average::serialize(Archive& ar, const unsigned int version)
   {
      AverageStage::serialize(ar, version);
      ar & blockSum_;
      ar & iBlock_;
      ar & nSamplePerBlock_;
   }

}
#endif
