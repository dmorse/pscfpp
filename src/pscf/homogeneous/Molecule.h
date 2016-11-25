#ifndef PSCF_HOMOGENEOUS_MOLECULE_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/chem/Species.h>           // base class
#include <util/param/ParamComposite.h>   // base class

#include <pscf/homogeneous/Group.h>      // member template argument
#include <util/containers/Pair.h>        // member template
#include <util/containers/DArray.h>      // member template

#include <cmath>

namespace Pscf { 
namespace Homogeneous { 

   using namespace Util;

   /**
   * Descriptor of a molecular species in a homogeneous mixture.
   *
   * \ingroup Pscf_Homogeneous_Module
   */
   class Molecule : public Species, public ParamComposite
   {

   public:

      /**
      * Constructor.
      */
      Molecule();
 
      /**
      * Destructor.
      */
      ~Molecule();

      /**
      * Read and initialize.
      *
      * \param in input parameter stream
      */
      virtual void readParameters(std::istream& in);

      #if 0
      /**
      * Compute chemical potential or volume fraction.
      *
      */ 
      virtual void compute(const DArray<WField>& wFields);
      #endif
 
      /// \name Accessors (objects, by reference)
      //@{

      /**
      * Get a specified Group.
      *
      * \param id group index, 0 <= id < nGroup
      */
      Group& group(int id);

      //@}
      /// \name Accessors (by value)
      //@{

      /**
      * Number of monomer groups (monomer types).
      */
      int nGroup() const; 

      /**
      * Total molecule size  = volume / reference volume.
      */
      double size() const;

      //@}

   private:

      /// Array of Group objects in this polymer.
      DArray<Group> groups_;

      /// Number of groups in this polymer
      int nGroup_;

      /// Total size of all groups (in units of reference size).
      double size_;

   };

   /*
   * Number of groups.
   */
   inline int Molecule::nGroup() const
   {  return nGroup_; }

   /*
   * Total size of all groups = volume / reference volume
   */
   inline double Molecule::size() const
   {  return size_; }

   /*
   * Get a specified Group.
   */
   inline Group& Molecule::group(int id)
   {  return groups_[id]; }
 
}
}
#endif
