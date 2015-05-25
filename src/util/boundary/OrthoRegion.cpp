/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "OrthoRegion.h"
#include <util/space/Dimension.h>
#include <util/math/feq.h>
#include <util/global.h>

#include <sstream>

namespace Util
{

   /* 
   * Constructor. 
   */
   OrthoRegion::OrthoRegion() 
   {
      for (int i = 0; i < Dimension; ++i) {
         minima_[i]      = 0.0;
         maxima_[i]      = 1.0;
         lengths_[i]     = 1.0;
         halfLengths_[i] = 0.5;
      };
      volume_ = 1.0;
   }

   /* 
   * Set lengths_, halfLengths_, and volume_ to consistent values.
   */
   void OrthoRegion::resetRegion() 
   {
      for (int i = 0; i < Dimension; ++i) {
         assert(maxima_[i] > minima_[i]);
         lengths_[i] = maxima_[i] - minima_[i];
         halfLengths_[i] = 0.5*lengths_[i];
      };
      volume_ = lengths_[0] * lengths_[1] * lengths_[2];
   }

   /* 
   * Check consistency of data.
   */
   bool OrthoRegion::isValid() 
   { 
      for (int i = 0; i < Dimension; ++i) {
         if (maxima_[i] <= minima_[i]) {
            UTIL_THROW("maxima_[i] <= minima_[i]");
         }
         if (!feq(lengths_[i], maxima_[i] - minima_[i])) {
            UTIL_THROW("lengths_[i] != maxima_[i] - minima_[i]");
         }
         if (!feq(halfLengths_[i], 0.5*lengths_[i])) {
            UTIL_THROW("halfLengths_[i] != 0.5*lengths_[i]");
         }
      }
      double product = lengths_[0]*lengths_[1]*lengths_[2];
      double diff = volume_ - product;
      if (!feq(volume_, product)) {
         std::stringstream buf;
         buf.precision(12);
         buf << "volume_ != product of lengths_" << std::endl;
         buf << "volume_     = " << volume_ << std::endl;
         buf << "lengths_[0] = " << lengths_[0] << std::endl;
         buf << "lengths_[1] = " << lengths_[1] << std::endl;
         buf << "lengths_[2] = " << lengths_[2] << std::endl;
         buf << "product     = " << product;
         buf << "diff        = " << diff;
         UTIL_THROW(buf.str().c_str());
      }
      return true;
   }

} 
