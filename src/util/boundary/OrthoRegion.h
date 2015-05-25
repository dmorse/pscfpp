#ifndef UTIL_ORTHO_REGION_H
#define UTIL_ORTHO_REGION_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/space/Vector.h>

namespace Util
{

   /**
   * A region with orthogonal edges parallel to the x, y, and z axes.
   *
   * This struct is intended to be used as a private or protected 
   * component, i.e., a private or protected member or base class of 
   * another class.
   *
   * \ingroup Boundary_Module
   */
   struct OrthoRegion
   {

   public:

      /// Minimum coordinates: Require r[i] >= minima_[i].
      Vector minima_;

      /// Maximum coordinates: Require r[i] <  maxima_[i].
      Vector maxima_;

      /// OrthoRegion lengths:  lengths_[i] = maxima_[i] - minima_[i].
      Vector lengths_;

      /// Half region lengths: halfLengths_[i] = 0.5*lengths_[i].
      Vector halfLengths_;

      /// Volume: V = lengths_[0]*lengths_[1]*lengths_[2].
      double volume_;

      /**
      * Constructor.
      *
      * Initialize to minima_ = {0.0, 0.0, 0.0} and maxima_ = {1.0, 1.0, 1.0}
      */
      OrthoRegion();

      /**
      * Set lengths and volume to values consistent with minima and maxima.
      */
      void resetRegion();

      /**
      * Return true if valid, or throw Exception.
      */
      bool isValid();

      /**
      * Serialize to/from an archive.
      * 
      * \param ar       saving or loading archive
      * \param version  archive version id
      */
      template <class Archive>
      void serialize(Archive& ar, const unsigned int version);

   };

   /*
   * Serialize an OrthoRegion to/from an archive.
   */
   template <class Archive>
   void OrthoRegion::serialize(Archive& ar, const unsigned int version)
   {
      ar & minima_;
      ar & maxima_;
      ar & lengths_;
      ar & halfLengths_;
      ar & volume_;
   }


}
#endif
