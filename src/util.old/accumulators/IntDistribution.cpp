/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "IntDistribution.h"
#include <util/format/Int.h>
#include <util/global.h>

namespace Util
{

   /* 
   * Constructor.
   */
   IntDistribution::IntDistribution() 
    : histogram_(),
      min_(0),
      max_(0),
      nBin_(0),
      nSample_(0),
      nReject_(0)
   {  setClassName("IntDistribution"); }
   
   /* 
   * Copy constructor.
   */
   IntDistribution::IntDistribution(const IntDistribution& other) 
    : histogram_(),
      min_(other.min_),
      max_(other.max_),
      nBin_(other.nBin_),
      nSample_(other.nSample_),
      nReject_(other.nReject_)
   {
      if (nBin_ > 0) {
         assert(nBin_ == max_ - min_ + 1);
         assert(other.histogram_.capacity() != 0);
         histogram_.allocate(nBin_);
         for (int i=0; i < nBin_; ++i) {
            histogram_[i] = other.histogram_[i];
         }
      } else {
         assert(other.histogram_.capacity() == 0);
         assert(min_ == 0);
         assert(max_ == 0);
         assert(nSample_ == 0);
         assert(nReject_ == 0);
      }
   }
   
   /* 
   * Assignment operator.
   */
   IntDistribution& IntDistribution::operator = (const IntDistribution& other) 
   {
      // Check for self assignment
      if (this == &other) return *this;

      // Check validity of other object
      if (other.nBin_ > 0) {
         assert(other.nBin_ == other.max_ - other.min_ + 1);
         assert(other.histogram_.capacity() != 0);
      } else {
         assert(other.nBin_ == 0);
         assert(other.histogram_.capacity() == 0);
         assert(other.min_ == 0);
         assert(other.max_ == 0);
         assert(other.nSample_ == 0);
         assert(other.nReject_ == 0);
      }

      // Assign primitive values
      min_      = other.min_;
      max_      = other.max_;
      nBin_     = other.nBin_;
      nSample_  = other.nSample_;
      nReject_  = other.nReject_;

      // Allocate histogram if necessary
      if (nBin_ > 0) {
         histogram_.allocate(nBin_);
         for (int i=0; i < nBin_; ++i) {
            histogram_[i] = other.histogram_[i];
         }
      }

      return *this;
   }
   
   /* 
   * Destructor.
   */
   IntDistribution::~IntDistribution() 
   {}

   /* 
   * Read min and max from file.
   */
   void IntDistribution::readParameters(std::istream& in)
   {
      read<int>(in, "min", min_);
      read<int>(in, "max", max_);
      nBin_ = max_ - min_ + 1;
      histogram_.allocate(nBin_);
      clear();
   }
   
   /*
   * Set parameters and initialize.
   */
   void IntDistribution::setParam(int min, int max)
   {
      min_      = min;
      max_      = max;
      nBin_     = max_ - min_ + 1;
      histogram_.allocate(nBin_);
      clear();
   }  
   
   /*
   * Load state from an archive.
   */
   void IntDistribution::loadParameters(Serializable::IArchive& ar)
   {  
      loadParameter<int>(ar, "min", min_);
      loadParameter<int>(ar, "max", max_);
      ar & nBin_;
      ar & nSample_;
      ar & nReject_;
      ar & histogram_;

      // Validate
      if (nBin_ != max_ - min_ + 1) {
         UTIL_THROW("Inconsistent values of nBin");
      }
      if (nBin_ != histogram_.capacity()) {
         UTIL_THROW("Inconsistent histogram capacity");
      }
   }

   /*
   * Save state to an archive.
   */
   void IntDistribution::save(Serializable::OArchive& ar)
   {  ar & *this; }

   /* 
   * Zero all accumulators.
   */
   void IntDistribution::clear()
   {  
      nSample_ = 0; 
      nReject_ = 0; 
      for (int i=0; i < nBin_; ++i) {
         histogram_[i] = 0;
      }
   }
   
   /* 
   * Add a value to the histogram
   */
   void IntDistribution::sample(int value)
   {
      int i;
      if (value >= min_ && value <= max_) {
         i = binIndex(value);
         histogram_[i] += 1;
         nSample_ += 1;
      } else {
         nReject_ += 1;
      }
   }
   
   /* 
   * Output histogram
   */
   void IntDistribution::output(std::ostream& out) 
   {
      for (int i=0; i < nBin_; ++i) {
         out << Int(i + min_) << Int(histogram_[i]) << std::endl;
      }
   }
  
   #if 0 
   /* 
   *
   */
   void IntDistribution::backup(FILE *file) 
   {
      fprintf(file, "nSample       %i \n", nSample_);
      fprintf(file, "nReject       %i \n", nReject_);
      for (int i=0; i < nBin_; ++i) {
         fprintf(file, "%li ", histogram_[i]);
      }
      fprintf(file, "\n");
   }
   
   /* 
   *
   */
   void IntDistribution::restore(FILE *file) {
      fscanf(file, "nSample       %i \n", &nSample_);
      fscanf(file, "nReject       %i \n", &nReject_);
      for (int i=0; i < nBin_; ++i) {
         fscanf(file, "%li ", &histogram_[i]);
      }
      fscanf(file, "\n");
   }
   #endif

}
