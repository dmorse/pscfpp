#ifndef UTIL_SERIALIZE_H
#define UTIL_SERIALIZE_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/global.h>

namespace Util
{

   /**
   * Serialize one object of type T.
   *
   * Default implementation calls serialize method of data object.
   * Can be overridden by any explicit specialization.
   *
   * \ingroup Serialize_Module
   *
   * \param ar      archive object
   * \param data    object to be serialized
   * \param version archive version id
   */
   template <class Archive, typename T>
   inline void serialize(Archive& ar, T& data, const unsigned int version)
   {  data.serialize(ar, version); }

   /**
   * Serialize an enumeration value.
   *
   * \ingroup Serialize_Module
   *
   * \param ar      archive object
   * \param data    object to be serialized
   * \param version archive version id
   */
   template <class Archive, typename T>
   inline void serializeEnum(Archive& ar, T& data, const unsigned int version = 0)
   {  
      unsigned int i;
      if (Archive::is_saving()) {
         i = (unsigned int)data;
      }
      ar & i;
      if (Archive::is_loading()) {
         data = (T)i;
      }
   }

   /**
   * Save a value, or save and check correctness on loading.
   *
   * \ingroup Serialize_Module
   *
   * \param ar      archive object
   * \param data    object to be serialized
   * \param label   label C-string for object.
   */
   template <class Archive, typename T>
   inline void serializeCheck(Archive& ar, T& data, const char* label = "")
   {
      T temp;
      if (Archive::is_saving()) {
         temp = data;
      }
      ar & temp;
      if (Archive::is_loading()) {
         if (temp != data) {
            UTIL_THROW("Inconsistent data");
         }
      }
   }

}
#endif
