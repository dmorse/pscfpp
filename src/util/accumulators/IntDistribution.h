#ifndef UTIL_INT_DISTRIBUTION_H
#define UTIL_INT_DISTRIBUTION_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/ParamComposite.h>
#include <util/containers/DArray.h>

namespace Util
{

   /**
   * A distribution (or histogram) of values for an int variable.
   * 
   * \ingroup Accumulators_Module
   */
   class IntDistribution : public ParamComposite 
   {
   
   public:
  
      /** 
      * Default constructor
      */
      IntDistribution();
   
      /** 
      * Copy constructor.
      *
      * \param other object to be copied
      */
      IntDistribution(const IntDistribution& other);
   
      /** 
      * Assignment operator.
      *
      * \param other object to be assigned
      */
      IntDistribution& operator = (const IntDistribution& other);
   
      /** 
      * Destructor
      */
      virtual ~IntDistribution();
   
      /**
      * Read parameters from file and initialize.
      *
      * Read values of min, max, and nBin from file.
      * Allocate histogram array and clear all accumulators.
      *
      * \param in input parameter file stream
      */
      void readParameters(std::istream& in);
  
      /** 
      * Set parameters and initialize.
      *
      * \param min  lower bound of range
      * \param max  upper bound of range
      */
      void setParam(int min, int max);
   
      /**
      * Load state from an archive.
      *
      * \param ar binary loading (input) archive.
      */
      virtual void loadParameters(Serializable::IArchive& ar);

      /**
      * Save state to an archive.
      *
      * \param ar binary saving (output) archive.
      */
      virtual void save(Serializable::OArchive& ar);
  
      /**
      * Serialize to/from an archive. 
      * 
      * \param ar      archive
      * \param version archive version id
      */
      template <class Archive>
      void serialize(Archive& ar, const unsigned int version);

      /**
      * Clear (i.e., zero) previously allocated histogram.
      */
      void clear();
   
      /**
      * Sample a value.
      *
      * \param value current value
      */
      void sample(int value);
   
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
      int binIndex(int value);
  
      /** 
      * Get minimum value in range of histogram.
      */
      int min() const;
  
      /** 
      * Get maximum value in range of histogram.
      */
      int max() const;
    
      /** 
      * Get the number of bins
      */
      int nBin() const;
        
   protected:
   
      DArray<long> histogram_;  ///< Histogram array.
      int   min_;               ///< minimum value.
      int   max_;               ///< maximum value.
      int   nBin_;              ///< number of bins.
      int   nSample_;           ///< Number of sampled values in Histogram.
      int   nReject_;           ///< Number of sampled values that were out of range.
   
   };

   // inline method definitions 
  
   /*
   * Return the index of the bin for a value.
   */
   inline int IntDistribution::binIndex(int value) 
   { return (value - min_); }
  
   /*
   * Get minimum value in range of histogram.
   */
   inline int IntDistribution::min() const
   {  return min_; }
  
   /*
   * Get maximum value in range of histogram.
   */
   inline int IntDistribution::max() const
   {  return max_; }
 
   /*
   * Get the number of bins
   */
   inline int IntDistribution::nBin() const
   {  return nBin_; }
        
   /*
   * Serialize this Distribution.
   */
   template <class Archive>
   void IntDistribution::serialize(Archive& ar, const unsigned int version)
   {
      ar & min_;        
      ar & max_;    
      ar & nBin_;     
      ar & nSample_;   
      ar & nReject_;    
      ar & histogram_; 
   }

}
#endif
