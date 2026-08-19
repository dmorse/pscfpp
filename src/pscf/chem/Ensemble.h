#ifndef PSCF_ENSEMBLE_H
#define PSCF_ENSEMBLE_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/ParamComposite.h>   // base class
#include <util/global.h>

namespace Pscf {

   using namespace Util;

   /**
   * Statistical ensemble type for the number of molecules of one species.
   *
   * \ingroup Pscf_Chem_Module
   */
   enum class Ensemble {Unknown, Closed, Open};

   /**
   * istream extractor for an Ensemble enumeration.
   *
   * \param  in       input stream
   * \param  policy   Ensemble to be read
   * \return modified input stream
   */
   std::istream& operator >> (std::istream& in, Ensemble& policy);

   /**
   * ostream inserter for an Ensemble enumeration.
   *
   * Text representations of allowed values are "Open" and "Closed".
   *
   * \param  out      output stream
   * \param  policy   Ensemble to be written
   * \return modified output stream
   */
   std::ostream& operator << (std::ostream& out, Ensemble policy);

   /**
   * Serialize an Ensemble
   *
   * \param ar      archive object
   * \param policy  object to be serialized
   * \param version archive version id
   */
   template <class Archive>
   void serialize(Archive& ar, Ensemble& policy,
                  const unsigned int version)
   {  serializeEnum(ar, policy, version); }

}
#endif
