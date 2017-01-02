/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Distribution.h"
#include <util/format/Dbl.h>
#include <util/global.h>

namespace Util
{

   /* 
   * Default constructor.
   */
   Distribution::Distribution() 
    : histogram_(),
      min_(0.0),
      max_(0.0),
      binWidth_(0.0),
      nBin_(0),
      nSample_(0),
      nReject_(0)
   {  setClassName("Distribution"); }
   
   /* 
   * Copy constructor.
   */
   Distribution::Distribution(const Distribution& other) 
    : histogram_(),
      min_(other.min_),
      max_(other.max_),
      binWidth_(other.binWidth_),
      nBin_(other.nBin_),
      nSample_(other.nSample_),
      nReject_(other.nReject_)
   {
      if (nBin_ > 0) {
         assert(other.histogram_.capacity() != 0);
         histogram_.allocate(nBin_);
         for (int i=0; i < nBin_; ++i) {
            histogram_[i] = other.histogram_[i];
         }
      } else {
         assert(nBin_ == 0);
         assert(histogram_.capacity() == 0);
         assert(nSample_ == 0);
         assert(nReject_ == 0);
      }
   }
   
   /* 
   * Assignment operator.
   */
   Distribution& Distribution::operator = (const Distribution& other) 
   {
      // Check for self assignment
      if (this == &other) return *this;

      // Check validity of other object
      if (other.nBin_ > 0) {
         assert(other.histogram_.capacity() != 0);
      } else {
         assert(other.nBin_ == 0);
         assert(other.histogram_.capacity() == 0);
         assert(other.nSample_ == 0);
         assert(other.nReject_ == 0);
      }

      // Assign primitive values
      min_      = other.min_;
      max_      = other.max_;
      binWidth_ = other.binWidth_;
      nBin_     = other.nBin_;
      nSample_  = other.nSample_;
      nReject_  = other.nReject_;

      // Allocate and copy histogram, if necessary
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
   Distribution::~Distribution() 
   {}

   /* 
   * Read parameters and initialize.
   */
   void Distribution::readParameters(std::istream& in)
   {
      read<double>(in, "min", min_);
      read<double>(in, "max", max_);
      read<int>(in,   "nBin", nBin_);
      binWidth_  = (max_ - min_)/double(nBin_);
      histogram_.allocate(nBin_);
      clear();
   }
  
   /*
   * Set parameters and initialize.
   *
   * \param min  lower bound of range
   * \param max  upper bound of range
   * \param nBin number of bins in range [min, max]
   */
   void Distribution::setParam(double min, double max, int nBin)
   {
      min_   = min;
      max_   = max;
      nBin_  = nBin;
      binWidth_  = (max_ - min_)/double(nBin_);
      histogram_.allocate(nBin_);
      clear();
   }  
   
   /*
   * Load internal state from archive.
   */
   void Distribution::loadParameters(Serializable::IArchive &ar)
   {
      loadParameter<double>(ar, "min", min_);
      loadParameter<double>(ar, "max", max_);
      loadParameter<int>(ar, "nBin", nBin_);
      ar & nSample_;
      ar & nReject_;
      ar & binWidth_;
      ar & histogram_;

      // Validate
      if (histogram_.capacity() != nBin_) {
         UTIL_THROW("Inconsistent histogram capacity");
      }
      if (!feq(binWidth_, (max_ - min_)/double(nBin_))) {
         UTIL_THROW("Inconsistent binWidth_");
      }
   }

   /*
   * Save internal state to archive.
   */
   void Distribution::save(Serializable::OArchive &ar)
   { ar & *this; }
   
   /* 
   * Zero all accumulators.
   */
   void Distribution::clear()
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
   void Distribution::sample(double value)
   {
      int i;
      if (value > min_ && value < max_) {
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
   void Distribution::output(std::ostream& out) 
   {
      double x, rho;
      for (int i=0; i < nBin_; ++i) {
         x   = min_ + binWidth_*(double(i) + 0.5);
         rho = double(histogram_[i])/double(nSample_);
         rho = rho/binWidth_;
         out << Dbl(x, 18, 8) << Dbl(rho, 18, 8) << std::endl;
      }
   }

   #ifdef UTIL_MPI
   /*
   * Reduce (add) distributions from multiple MPI processors.
   */
   void Distribution::reduce(MPI::Intracomm& communicator, int root)
   {
  
      long* totHistogram = new long[nBin_]; 
      communicator.Reduce(histogram_.cArray(), totHistogram, nBin_, MPI::LONG, MPI::SUM, root);
      if (communicator.Get_rank() == root) {
         for (int i=0; i < nBin_; ++i) {
            histogram_[i] = totHistogram[i];
         }
      } else { 
         for (int i=0; i < nBin_; ++i) {
            histogram_[i] = 0.0;
         }
      }
      delete totHistogram;

      long totSample; 
      communicator.Reduce(&nSample_, &totSample, 1, MPI::LONG, MPI::SUM, root);
      if (communicator.Get_rank() == root) {
         nSample_ = totSample;
      } else {
         nSample_ = 0;
      }

      long totReject; 
      communicator.Reduce(&nReject_, &totReject, 1, MPI::LONG, MPI::SUM, root);
      if (communicator.Get_rank() == root) {
         nReject_ = totReject;
      } else {
         nReject_ = 0;
      }

   }
   #endif
        
}
