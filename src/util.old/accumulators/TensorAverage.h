#ifndef UTIL_TENSOR_AVERAGE_H
#define UTIL_TENSOR_AVERAGE_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/global.h>
#include <util/param/ParamComposite.h>  // base class
#include <util/accumulators/Average.h>  // member template argument
#include <util/containers/FArray.h>     // member template

#include <vector>

namespace Util
{

   class Tensor;

   /**
   * Calculates averages of all components of a Tensor-valued variable.
   *
   * TensorAverage is a simple container for an array of Average objects,
   * each of which calculates averages and error estimates for one
   * component of a Tensor.
   *
   * \ingroup Accumulators_Module
   */
   class TensorAverage : public ParamComposite
   {

   public:

      /**
      * Constructor
      *
      * \param blockFactor ratio of block sizes for subsequent stages.
      */
      TensorAverage(int blockFactor = 2);

      /**
      * Destructor
      */
      virtual ~TensorAverage();

      /**
      * Set nSamplePerBlock.
      *
      * If nSamplePerBlock > 0, the sample function will increment block
      * averages, and reset the average every nSamplePerBlock samples.
      *
      * If nSamplePerBlock == 0, block averaging is disabled. This is the
      * default (i.e., the initial value set in the constructor).
      *
      * \param nSamplePerBlock  number of samples per block average output
      */
      void setNSamplePerBlock(int nSamplePerBlock);

      /**
      * Read parameter nSamplePerBlock from file and initialize.
      *
      * See setNSamplePerBlock() for discussion of value.
      *
      * \param in  input stream
      */
      void readParameters(std::istream& in);

      /**
      * Load internal state from an archive.
      *
      * \param ar  input/loading archive
      */
      virtual void loadParameters(Serializable::IArchive &ar);

      /**
      * Save internal state to an archive.
      *
      * \param ar  output/saving archive
      */
      virtual void save(Serializable::OArchive &ar);

      /**
      * Serialize this to or from an archive.
      *
      * \param ar  input or output archive
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
      * \param value  sampled value
      */
      void sample(const Tensor& value);

      /**
      * Access the Average object for one tensor component.
      *
      * \param i  first index of associated tensor component
      * \param j  second index of associated tensor component
      * \return  Average object associated with element (i, j)
      */
      const Average& operator () (int i, int j);

      /**
      * Get number of samples per block average.
      *
      * Returns zero if block averaging is disabled.
      *
      * \return  number of samples per block (or 0 if disabled).
      */
      int nSamplePerBlock() const;

      /**
      * Get number of samples in current block average.
      *
      * Returns 0 if block averaging is disabled (i.e., nSamplePerBlock == 0).
      *
      * \return  number of samples in current block (or 0 if disabled)
      */
      int iBlock() const;

      /**
      * Is the current block average complete?
      *
      * Returns true iff blocking is enabled and iBlock == nSamplePerBlock
      *
      * \return (iBlock > 0) && (iBlock == nSamplePerBlock)
      */
      bool isBlockComplete() const;

   private:

      /// Array of average accumulators, one per tensor component.
      FArray<Average, Dimension*Dimension> accumulators_;

      /// Number of sampled values per output block.
      int nSamplePerBlock_;

      /// Number of samples in current output block.
      int iBlock_;

      /// Private and not implemented to prohibit copying.
      TensorAverage(const TensorAverage& other);

      /// Private and not implemented to prohibit assignment.
      TensorAverage& operator = (const TensorAverage& other);

   };

   // Inline method definitions

   /*
   * Get nSamplePerBlock, number of samples per block average.
   */
   inline int TensorAverage::nSamplePerBlock() const
   {  return nSamplePerBlock_; }

   /*
   * Get iBlock, number of samples in current block average.
   */
   inline int TensorAverage::iBlock() const
   {  return iBlock_; }

   /*
   * Is the current block average complete?
   */
   inline bool TensorAverage::isBlockComplete() const
   { return (iBlock_ && (iBlock_ == nSamplePerBlock_)); }

   /*
   * Serialize this Average.
   */
   template <class Archive>
   void TensorAverage::serialize(Archive& ar, const unsigned int version)
   {
      ar & nSamplePerBlock_;
      ar & iBlock_;
      int i, j, k;
      k = 0;
      for (i = 0; i < Dimension; ++i) {
         for (j = 0; j < Dimension; ++j) {
            ar & accumulators_[k];
            ++k;
         }
      }
   }

}
#endif
