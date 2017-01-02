#ifndef UTIL_AUTO_CORR_ARRAY_H
#define UTIL_AUTO_CORR_ARRAY_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/ParamComposite.h>       // base class
#include <util/containers/DArray.h>          // member template
#include <util/containers/RingBuffer.h>      // member template parameter

#include <util/accumulators/setToZero.h>
#include <util/accumulators/product.h>
#include <util/containers/Array.h>      
#include <util/space/Vector.h>
#include <util/format/Int.h>
#include <util/format/Dbl.h>

#include <complex>
using std::complex;

namespace Util
{

   template <typename> class Array;
   
   /**
   * Auto-correlation function for an ensemble of sequences.
   *   
   * This class calculates an autocorrelation function for a ensemble of
   * statistically equivalent sequences x(i) of values of a variable of 
   * type Data. The resulting autocorrelation function is an array of 
   * values of type Product, where C(j) = <x(i-j), x(i)>. Here <A,B> 
   * denotes an inner product of type Product for objects A and B of 
   * type Data.
   *
   * The meaning of <A,B>  for two Data values is defined for various
   * data types by the overloaded functions product(Data, Data) defined
   * in file "product.h" . These functions define a product as an
   * arithmetic product for floating point numbers, and use the 
   * following definitions for complex numbers and Vector objects:
   *
   *     double  product(double, double)   = A*B
   *     complex product(complex, complex) = conjug(A)*B
   *     double  product(Vector,  Vector)  = A.dot(B)
   *
   * The meaning of setting a variable to zero is defined for various
   * types of data by the overloaded functions setToZero(Data&) that 
   * are defined in file setToZero.h.
   *
   * \ingroup Accumulators_Module
   */
   template <typename Data, typename Product>
   class AutoCorrArray : public ParamComposite
   {
   
   public:
   
      /// Default constructor.
      AutoCorrArray();

       /// Default destructor.
      ~AutoCorrArray();
  
      /**
      * Read parameters, allocate memory and clear history.
      *
      * Reads parameters ensembleCapacity and bufferCapacity, allocates
      * memory, sets nEnsemble=ensembleCapacity, and calls clear().
      *
      * \param in input parameter stream
      */
      virtual void readParameters(std::istream& in);
   
      /**
      * Allocate memory, and clear history.
      *
      * Sets parameters ensembleCapacity and bufferCapacity, allocates
      * memory, sets nEnsemble=ensembleCapacity, and calls clear().
      *
      * \param ensembleCapacity maximum number of sequences in ensemble
      * \param bufferCapacity   maximum number of values in each history
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
      * Set actual number of sequences in ensemble.
      *
      * \pre readParam() or setParam() must have been called previously
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
      * Serialize this AutoCorrArray to/from an archive.
      *
      * \param ar       input or output archive
      * \param version  file version id
      */
      template <class Archive>
      void serialize(Archive& ar, const unsigned int version);

      /** 
      * Return maximum number of samples in history for each sequence.
      */
      int bufferCapacity() const;
  
      /** 
      * Return nEnsemble.
      */
      int nEnsemble() const;
  
      /** 
      * Return the total number of samples per sequence thus far.
      */
      int nSample() const;
  
      /** 
      * Return average of sampled values.
      */
      Data average() const;
  
      /**
      * Numerical integral of autocorrelation function 
      */
      double corrTime() const; 
      
   private:
   
      /// Array of Ring buffers, each containing a sequence of values.
      DArray< RingBuffer<Data> > buffers_; 
   
      /// Array in which corr[j] = sum of values of <x(i-j), x(i)>
      DArray<Product> corr_;
   
      /// Array in which nCorr[i] = number of values added to corr[i]
      DArray<int> nCorr_;
   
      /// Sum of all previous values of x(t).
      Data sum_;
   
      /// Maximum number of sequences for which memory is allocated
      int  ensembleCapacity_;
 
      /// Capacity (# of elements) of corr, nCorr and each buffer.
      int  bufferCapacity_;
      
      /// Total number of sequences
      int  nEnsemble_;
   
      /// Total number of previous values of x(t) per sequence
      int  nSample_;

      /**
      * Allocate memory and call clear.
      *
      * Called by readParam and setParam. 
      *
      * Precondition: nEnsemble_ and bufferCapacity_ must be set
      */
      void allocate();
   
   };

   /*
   * Default constructor.
   */
   template <typename Data, typename Product>
   AutoCorrArray<Data, Product>::AutoCorrArray() 
    : buffers_(),
      corr_(),
      nCorr_(),
      ensembleCapacity_(0),
      bufferCapacity_(0),
      nEnsemble_(0),
      nSample_(0)
   {
      setClassName("AutoCorrArray"); 
      setToZero(sum_); 
   }

   /*
   * Destructor.
   */
   template <typename Data, typename Product>
   AutoCorrArray<Data, Product>::~AutoCorrArray() 
   {}
 
   /*
   * Read parameters from file.
   */
   template <typename Data, typename Product>
   void AutoCorrArray<Data, Product>::readParameters(std::istream& in)
   {
      read<int>(in, "ensembleCapacity", ensembleCapacity_);
      read<int>(in, "bufferCapacity",   bufferCapacity_);
      nEnsemble_ = ensembleCapacity_;
      allocate();
   }
   
   /*
   * Set parameters and initialize.
   */
   template <typename Data, typename Product>
   void AutoCorrArray<Data, Product>::setParam(int ensembleCapacity, int bufferCapacity)
   {
      ensembleCapacity_ = ensembleCapacity;
      bufferCapacity_ = bufferCapacity;
      allocate();
      nEnsemble_ = ensembleCapacity;
   }

   /*
   * Load internal state from archive.
   */
   template <typename Data, typename Product>
   void AutoCorrArray<Data, Product>::loadParameters(Serializable::IArchive &ar)
   {
      loadParameter<int>(ar, "ensembleCapacity", ensembleCapacity_);
      loadParameter<int>(ar, "bufferCapacity",   bufferCapacity_);
      ar & nEnsemble_;
      ar & buffers_; 
      ar & corr_;
      ar & nCorr_;
      ar & sum_;
      ar & nSample_;
   }

   /*
   * Save internal state to archive.
   */
   template <typename Data, typename Product>
   void AutoCorrArray<Data, Product>::save(Serializable::OArchive &ar)
   { ar & *this; }

   /*
   * Set or reset nEnsemble.
   */
   template <typename Data, typename Product>
   void AutoCorrArray<Data, Product>::setNEnsemble(int nEnsemble)
   {
      if (ensembleCapacity_ == 0) {
         UTIL_THROW("No memory has been allocated: ensembleCapacity_ == 0"); 
      }
      if (nEnsemble > ensembleCapacity_) {
         UTIL_THROW("nEnsemble > ensembleCapacity_");      
      }
      nEnsemble_ = nEnsemble;
   }

   /* 
   * Set accumulator to initial empty state.
   */
   template <typename Data, typename Product>
   void AutoCorrArray<Data, Product>::clear()
   {   
      setToZero(sum_);
      nSample_ = 0;
      if (bufferCapacity_ > 0) {
         int i;
         for (i = 0; i < bufferCapacity_; ++i) {
            setToZero(corr_[i]);
            nCorr_[i] = 0;
         }
         for (i = 0; i < ensembleCapacity_; ++i) {
            buffers_[i].clear();
         }
      }
   }
   
   /*
   * Allocate arrays and CyclicBuffer (private method).
   */
   template <typename Data, typename Product>
   void AutoCorrArray<Data, Product>::allocate()
   { 
      if (bufferCapacity_ > 0) { 
         // Allocate autocorrelation accumulators
         corr_.allocate(bufferCapacity_);
         nCorr_.allocate(bufferCapacity_);
   
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
   template <typename Data, typename Product>
   void AutoCorrArray<Data, Product>::sample(const Array<Data>& values)
   {
      int i, j;
      ++nSample_;

      for (i = 0; i < nEnsemble_; ++i) {
         sum_ += values[i];
         buffers_[i].append(values[i]);
      }
      
      int bufferSize = buffers_[0].size();
      for (j=0; j < bufferSize; ++j) {
         for (i = 0; i < nEnsemble_; ++i) {
            corr_[j] += product(buffers_[i][j], values[i]);
         }
         ++nCorr_[j];
      };
   }
   
   /*
   * Return capacity of history buffer for each sequence.
   */
   template <typename Data, typename Product>
   int AutoCorrArray<Data, Product>::bufferCapacity() const
   {  return bufferCapacity_; }

   /*
   * Return number of independent sequences.
   */
   template <typename Data, typename Product>
   int AutoCorrArray<Data, Product>::nEnsemble() const
   {  return nEnsemble_; }

   /*
   * Return number of sampled values.
   */
   template <typename Data, typename Product>
   int AutoCorrArray<Data, Product>::nSample() const
   {  return nSample_; }

   /*
   * Return average of sampled values.
   */
   template <typename Data, typename Product>
   Data AutoCorrArray<Data, Product>::average() const
   {
      Data ave = sum_;
      ave /= double(nSample_*nEnsemble_);
      return ave;
   }

   /*
   * Final output.
   */
   template <typename Data, typename Product>
   void AutoCorrArray<Data, Product>::output(std::ostream& outFile) 
   {
      Product autocorr;
      // Product aveSq;
      // Data    ave = average();
   
      // Calculate and output autocorrelation
      // aveSq = product(ave, ave);
      int bufferSize = buffers_[0].size();
      for (int i = 0; i < bufferSize; ++i) {
         autocorr = corr_[i]/double(nCorr_[i]*nEnsemble_);
         //autocorr = autocorr - aveSq;
         //outFile << Int(i) << Dbl(autocorr) << Int(nCorr_[i]) << std::endl;
         outFile << Int(i) << Dbl(autocorr) << std::endl;
      }
   }
   
   /*
   *  Return correlation time in unit of sampling interval  
   */
   template <typename Data, typename Product>
   double AutoCorrArray<Data, Product>::corrTime() const
   {
      Data    ave;
      Product aveSq, variance, autocorr, sum;
      int     bufferSize = buffers_[0].size();
   
      // Calculate average of sampled values
      ave  = sum_;
      ave /= double(nSample_*nEnsemble_);
      aveSq = product(ave, ave);
      variance = corr_[0]/double(nCorr_[0]*nEnsemble_);
      variance = variance - aveSq;
   
      // Sum over autocorrelation function
      setToZero(sum);
      for (int i = 1; i < bufferSize/2; ++i) {
         autocorr = corr_[i]/double(nCorr_[i]*nEnsemble_);
         autocorr = autocorr - aveSq;
         sum += autocorr;
      }
      sum /= variance;
      return sum; 
   }
  
   /*
   * Serialize this AutoCorrArray.
   */
   template <typename Data, typename Product>
   template <class Archive>
   void AutoCorrArray<Data, Product>::serialize(Archive& ar, 
                                                const unsigned int version)
   {
      ar & ensembleCapacity_;
      ar & bufferCapacity_;
      ar & nEnsemble_;
      ar & buffers_; 
      ar & corr_;
      ar & nCorr_;
      ar & sum_;
      ar & nSample_;
   }

}
#endif
