#ifndef UTIL_MEAN_SQ_DISP_ARRAY_H
#define UTIL_MEAN_SQ_DISP_ARRAY_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

// Needed in header
#include <util/param/ParamComposite.h>
#include <util/containers/RingBuffer.h>
#include <util/containers/Array.h>

// Needed in implementation
#include <util/accumulators/MeanSqDispArray.h>
#include <util/accumulators/setToZero.h>
#include <util/space/Vector.h>
#include <util/format/Int.h>
#include <util/format/Dbl.h>

#include <complex>
using std::complex;

namespace Util
{

   /**
   * Mean-squared displacement (MSD) vs. time for an ensembles of sequences.
   *   
   * This class calculates the mean-squared difference <|x(i) - x(i-j)|^2> for
   * an ensemble of statistically equivalent sequences x(i) of values of a
   * variable of type Data. The meaning of |a - b|^2 is defined for int, 
   * double, and Vector data by explicit specializations of the private 
   * method double sqDiff(Data&, Data).
   * 
   * \ingroup Accumulators_Module
   */
   template <typename Data>
   class MeanSqDispArray : public ParamComposite 
   {
   
   public:
      
      /* 
      * Constructor.
      */
      MeanSqDispArray();

      /* 
      * Destructor.
      */
      ~MeanSqDispArray();

      /**
      * Read parameters, allocate memory and clear history.
      *
      * Reads parameters nEnsemble and capacity, allocates memory,
      * and then calls clear().
      *
      * \param in input parameter stream
      */
      void readParameters(std::istream& in);
   
      /**
      * Set parameters, allocate memory, and clear history.
      *
      * Sets parameters nEnsemble and capacity, allocates
      * memory and then calls clear().
      *
      * \param ensembleCapacity number of sequence in ensemble
      * \param bufferCapacity   number of variable values per sequence
      */
      void setParam(int ensembleCapacity, int bufferCapacity);
   
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
      * Serialize this MeanSqDispArray to/from an archive.
      *
      * \param ar       input or output archive
      * \param version  file version id
      */
      template <class Archive>
      void serialize(Archive& ar, const unsigned int version);

      /**
      * Set actual number of sequences in ensemble.
      *
      * \pre readParameters() or setParam() must have been called previously
      * \pre nEnsemble <= ensembleCapacity
      *
      * \param nEnsemble actual number of sequences in ensemble
      */
      void setNEnsemble(int nEnsemble);
   
      /**
      * Reset to empty state.
      */
      void clear();
   
      /**
      * Sample an array of current values.
      *
      * \param values Array of current values
      */
      void sample(const Array<Data>& values);
   
      /**
      * Output the autocorrelation function
      */
      void output(std::ostream& out);
  
      /** 
      * Return capacity of the history buffer for each sequence.
      */
      int bufferCapacity() 
      {  return bufferCapacity_; }
  
      /** 
      * Return number of sequences in the ensemble.
      */
      int nEnsemble() 
      {  return nEnsemble_; }
  
      /** 
      * Return number of values sampled from each sequence thus far.
      */
      int nSample() 
      {  return nSample_; }
   
   private:

      /// Array of Ring buffers containing many sequences of stored values.
      DArray< RingBuffer<Data> >  buffers_; 
   
      /// Array in which sqDiffSums[j] = sum of values <|x(i) - x(i-j)|^2> .
      DArray<double>  sqDiffSums_;
   
      /// Array in which nValues_[i] = number of values added to sqDiffSums_[i].
      DArray<int>  nValues_;
   
      /// Maximum number of sequences in the ensemble.
      int  ensembleCapacity_;
      
      /// Physical capacity (# of elements) of each buffer, corr, and nValues.
      int  bufferCapacity_;
      
      /// Total number of sequences in the ensemble.
      int  nEnsemble_;
   
      /// Total number of previous values of x(t) per sequence.
      int  nSample_;

      /**
      * Allocate memory and call clear.
      *
      * Called within readParameters and setParam, after setting parameters.
      *
      * Precondition: nEnsemble_ and bufferCapacity_ must have been set.
      */
      void allocate();
   
      /**
      * Calculate square of the difference.
      *
      * Explicit specializations of this function are provided to define
      * the norm of the difference for double and Vector Data.
      *
      * \param  data1 current  data value.
      * \param  data2 previous data value.
      * \return Square difference |data1 - data2|^2
      */
      double sqDiff(const Data& data1, const Data& data2);
   
   };


   /*
   * Default constructor.
   */
   template <typename Data>
   MeanSqDispArray<Data>::MeanSqDispArray() 
    : buffers_(),
      sqDiffSums_(),
      nValues_(),
      ensembleCapacity_(0),
      bufferCapacity_(0),
      nEnsemble_(0),
      nSample_(0)
   {  setClassName("MeanSqDispArray"); }

   /*
   * Default destructor.
   */
   template <typename Data>
   MeanSqDispArray<Data>::~MeanSqDispArray() 
   {}
 
   /*
   * Read parameters from file.
   */
   template <typename Data>
   void MeanSqDispArray<Data>::readParameters(std::istream& in)
   {
      read<int>(in, "ensembleCapacity", ensembleCapacity_);
      read<int>(in, "bufferCapacity", bufferCapacity_);
      allocate();
      nEnsemble_ = ensembleCapacity_;
   }
   
   /*
   * Set parameters and initialize.
   */
   template <typename Data> 
   void 
   MeanSqDispArray<Data>::setParam(int ensembleCapacity, int bufferCapacity)
   {
      ensembleCapacity_ = ensembleCapacity;
      bufferCapacity_ = bufferCapacity;
      allocate();
      // Set number of sequences to maximum capacity as a default.
      nEnsemble_  = ensembleCapacity;
   }

   /*
   * Set or reset nEnsemble.
   */
   template <typename Data>
   void MeanSqDispArray<Data>::setNEnsemble(int nEnsemble)
   {
      if (ensembleCapacity_ == 0) 
         UTIL_THROW("No memory has been allocated: ensembleCapacity_ == 0"); 
      if (nEnsemble > ensembleCapacity_) 
         UTIL_THROW("nEnsemble > ensembleCapacity_");
      nEnsemble_ = nEnsemble;
   }

   /*
   * Load internal state from archive.
   */
   template <typename Data>
   void MeanSqDispArray<Data>::loadParameters(Serializable::IArchive &ar)
   {
      loadParameter<int>(ar, "ensembleCapacity", ensembleCapacity_);
      loadParameter<int>(ar, "bufferCapacity", bufferCapacity_);
      ar & nEnsemble_;
      ar & nValues_;
      ar & nSample_;
      ar & buffers_; 
      ar & sqDiffSums_;
   }

   /*
   * Serialize this MeanSqDispArray.
   */
   template <typename Data>
   template <class Archive>
   void MeanSqDispArray<Data>::serialize(Archive& ar, const unsigned int version)
   {
      ar & ensembleCapacity_;
      ar & bufferCapacity_;
      ar & nEnsemble_;
      ar & nValues_;
      ar & nSample_;
      ar & buffers_; 
      ar & sqDiffSums_;
   }

   /*
   * Save internal state to archive.
   */
   template <typename Data>
   void MeanSqDispArray<Data>::save(Serializable::OArchive &ar)
   { ar & *this; }
   
   /* 
   * Set previously allocated to initial empty state.
   */
   template <typename Data>
   void MeanSqDispArray<Data>::clear()
   {   
      nSample_ = 0;

      if (bufferCapacity_ > 0) {

         int i;
         for (i = 0; i < bufferCapacity_; ++i) {
            setToZero(sqDiffSums_[i]);
            nValues_[i] = 0;
         }
   
         for (i = 0; i < ensembleCapacity_; ++i) {
            buffers_[i].clear();
         }

      }

   }
   
   /*
   * Allocate arrays and an array of CyclicBuffer objects (private method).
   */
   template <typename Data>
   void MeanSqDispArray<Data>::allocate()
   { 
      if (bufferCapacity_ > 0) { 
 
         // Allocate accumulator arrays
         sqDiffSums_.allocate(bufferCapacity_);
         nValues_.allocate(bufferCapacity_);
   
         // Allocate buffers
         buffers_.allocate(ensembleCapacity_);
         for (int i=0; i < ensembleCapacity_; ++i) {
            buffers_[i].allocate(bufferCapacity_);
         }

      }
      clear();
   }
   
   /*
   * Sample a single value from a time sequence.
   */
   template <typename Data>
   inline void MeanSqDispArray<Data>::sample(const Array<Data>& values)
   {
      int i, j;

      ++nSample_;

      for (i = 0; i < nEnsemble_; ++i) {
         buffers_[i].append(values[i]);
      }
      
      int bufferSize = buffers_[0].size();
      for (j = 0; j < bufferSize; ++j) {
         for (i = 0; i < nEnsemble_; ++i) {
            sqDiffSums_[j] += sqDiff(buffers_[i][j], values[i]);
         }
         ++nValues_[j];
      };
   }
  
   /**
   * Square difference for integer data = double(|data1 - data2|^2) .
   *
   * \param data1 first integer
   * \param data2 second integer
   */ 
   template <>
   inline double 
   MeanSqDispArray<int>::sqDiff(const int& data1, const int& data2) 
   {
      int diff; 
      diff = data1 - data2;
      return double(diff*diff);
   }

   /**
   * Square difference for double data = |data1 - data2|^2.
   *
   * \param data1 first value
   * \param data2 second value
   */ 
   template <>
   inline double 
   MeanSqDispArray<double>::sqDiff(const double& data1, const double& data2) 
   {
      double diff; 
      diff = data1 - data2;
      return diff*diff;
   }

   /**
   * Square difference for Vector data = |data1 - data2|^2 .
   *
   * \param data1 first vector
   * \param data2 second vector
   */
   template <>
   inline double 
   MeanSqDispArray<Vector>::sqDiff(const Vector& data1, const Vector& data2) 
   {
      Vector diff;
      diff.subtract(data1, data2);
      return diff.square();
   }

   /*
   * Final output
   */
   template <typename Data>
   void MeanSqDispArray<Data>::output(std::ostream& out) 
   {
      double msd;
   
      // Calculate and output mean-squared difference
      int bufferSize = buffers_[0].size();
      for (int i = 0; i < bufferSize; ++i) {
         msd = sqDiffSums_[i]/double(nValues_[i]*nEnsemble_);
         out << Int(i) << Dbl(msd) << std::endl;
      }
      
   }
   
}
#endif
