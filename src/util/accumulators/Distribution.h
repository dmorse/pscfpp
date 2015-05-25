#ifndef UTIL_DISTRIBUTION_H
#define UTIL_DISTRIBUTION_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/ParamComposite.h>  // base class
#include <util/containers/DArray.h>     // member template
#include <util/math/feq.h>              // Used in serialize template

namespace Util
{

   /**
   * A distribution (or histogram) of values for a real variable.
   *
   * \ingroup Accumulators_Module
   */
   class Distribution : public ParamComposite 
   {
   
   public:
  
      /** 
      * Default constructor.
      */
      Distribution();
   
      /** 
      * Copy constructor.
      *
      * \param other Distribution to be copied.
      */
      Distribution(const Distribution& other);
   
      /** 
      * Assignment operator.
      *
      * \param other Distribution to be assigned.
      */
      Distribution& operator = (const Distribution& other);
   
      /** 
      * Destructor.
      */
      virtual ~Distribution();
   
      /**
      * Read parameters from file and initialize.
      *
      * Read values of min, max, and nBin from file.
      * Allocate histogram array and clear all accumulators.
      *
      * \param in input parameter file stream
      */
      virtual void readParameters(std::istream& in);
  
      /** 
      * Set parameters and initialize.
      *
      * \param min  lower bound of range
      * \param max  upper bound of range
      * \param nBin number of bins in range [min, max]
      */
      void setParam(double min, double max, int nBin);
   
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
      * Serialize this Distribution to/from an archive.
      *
      * \param ar       input or output archive
      * \param version  file version id
      */
      template <class Archive>
      void serialize(Archive& ar, const unsigned int version);

      /**
      * Sample a value.
      *
      * \param value current value
      */
      void sample(double value);
   
      /**
      * Clear (i.e., zero) previously allocated histogram.
      */
      virtual void clear();
   
      /**
      * Output the distribution to file. 
      *
      * \param out output stream
      */
      void output(std::ostream& out);
  
      /**
      * Return the index of the bin for a value.
      *
      * \param value sampled value
      */
      int binIndex(double value) const;
  
      /** 
      * Get minimum value in range of histogram.
      */
      double min() const;
  
      /** 
      * Get maximum value in range of histogram.
      */
      double max() const;
    
      /** 
      * Get binWidth, the width of each bin.
      */
      double binWidth() const;
  
      /** 
      * Get the number of bins
      */
      int nBin() const;

      #ifdef UTIL_MPI
      /**
      * Reduce (add) distributions from multiple MPI processors.
      *
      * \param communicator MPI communicator
      * \param root rank of MPI root processor for reduction
      */
      void reduce(MPI::Intracomm& communicator, int root);
      #endif
        
   protected:
   
      DArray<long> histogram_;  ///< Histogram of occurences, one element per bin.
      double  min_;             ///< minimum value.
      double  max_;             ///< maximum value.
      double  binWidth_;        ///< width of bin = (max_-min_)/nBin_ .
      int     nBin_;            ///< number of bins.
      int     nSample_;         ///< Number of sampled values in Histogram.
      int     nReject_;         ///< Number of sampled values that were out of range.
   
   };

   // inline method definitions 
  
   /*
   * Return the index of the bin for a value.
   */
   inline int Distribution::binIndex(double value) const
   { return int( (value - min_)/binWidth_ ); }
  
   /*
   * Get minimum value in range of histogram.
   */
   inline double Distribution::min() const
   {  return min_; }
  
   /*
   * Get maximum value in range of histogram.
   */
   inline double Distribution::max() const
   {  return max_; }
 
   /*
   * Get binWidth, the width of each bin.
   */
   inline double Distribution::binWidth() const
   {  return binWidth_; }
  
   /*
   * Get the number of bins
   */
   inline int Distribution::nBin() const
   {  return nBin_; }

   /*
   * Serialize this Distribution.
   */
   template <class Archive>
   void Distribution::serialize(Archive& ar, const unsigned int version)
   {
      ar & min_;        
      ar & max_;    
      ar & nBin_;     
      ar & nSample_;   
      ar & nReject_;    
      ar & binWidth_;  
      ar & histogram_; 

      // Validate
      if (histogram_.capacity() != nBin_) {
         UTIL_THROW("Inconsistent histogram capacity");
      }
      if (!feq(binWidth_, (max_ - min_)/double(nBin_))) {
         UTIL_THROW("Inconsistent binWidth");
      }
   }

}
#endif
