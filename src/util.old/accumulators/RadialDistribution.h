#ifndef UTIL_RADIAL_DISTRIBUTION_H
#define UTIL_RADIAL_DISTRIBUTION_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/accumulators/Distribution.h>

namespace Util
{

   /**
   * Distribution (or histogram) of values for particle separations.
   * 
   * \ingroup Accumulators_Module
   */
   class RadialDistribution : public Distribution 
   {
   
   public:
   
      /**
      * Default constructor.
      */
      RadialDistribution();
   
      /**
      * Copy constructor.
      * 
      * \param other object to be copied.
      */
      RadialDistribution(const RadialDistribution& other);
   
      /**
      * Assignment.
      * 
      * \param other object to be assigned.
      */
      RadialDistribution& operator = (const RadialDistribution& other);

      // base class destructor is fine
   
      /**
      * Read values of min, max, and nBin from file.
      *
      * \param in input parameter file stream.
      */
      virtual void readParameters(std::istream& in);
   
      /** 
      * Set parameters and initialize.
      *
      * \param max  upper bound of range
      * \param nBin number of bins in range [min, max]
      */
      void setParam(double max, int nBin);
   
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
      * Serialize this RadialDistribution to/from an archive.
      *
      * \param ar       input or output archive
      * \param version  file version id
      */
      template <class Archive>
      void serialize(Archive& ar, const unsigned int version);

      /**
      * Clear all accumulators.
      */
      virtual void clear();

      /**
      * Mark the beginning of a "snapshot" (i.e., a sampled time step). 
      */
      void beginSnapshot();
   
      /**
      * Set the factor used to normalize the RDF before output.
      *
      * \param norm normalizing factor
      */
      void setNorm(double norm);
   
      /**
      * Set true to enable output of spatial integral of g(r).
      *
      * \param outputIntegral true to enable output of integral.
      */
      void setOutputIntegral(bool outputIntegral);

      /**
      * Output the distribution to file. 
      *
      * \param out pointer to output file
      */
      void output(std::ostream& out);
   
      /**
      * Get number of snapshots.
      */
      long nSnapshot();
 
   private:
   
      /**
      * Normalizing factor: concentration^2 * system volume.
      *
      * If all pairs are sampled twice, the norm should be the number of
      * primary particles in the system (i.e., the number of particles
      * that can appear as the first particle in the pair) times the 
      * concentration (# per volume) of secondary particles. If each 
      * pair is sampled only once, the norm should be half that.
      *
      */
      double norm_;
   
      /**
      * Number of "snapshots" (i.e., simulation steps) sampled.
      */
      long   nSnapshot_;
   
      /**
      * If set true, output volume integral of normalized RDF.
      */
      bool   outputIntegral_;
   
   };

   // Inline methods

   /*
   * Get number of snapshots.
   */
   inline long RadialDistribution::nSnapshot()
   { return nSnapshot_; }

   /*
   * Serialize this Distribution.
   */
   template <class Archive>
   void RadialDistribution::serialize(Archive& ar, const unsigned int version)
   {
      Distribution::serialize(ar, version);
      ar & norm_;
      ar & nSnapshot_;
      ar & outputIntegral_;
   }

}
#endif
