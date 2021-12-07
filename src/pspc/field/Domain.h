#ifndef PSPC_DOMAIN_H
#define PSPC_DOMAIN_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2019, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/ParamComposite.h>     // base class

#include <pspc/field/FieldIo.h>            // member
#include <pspc/field/FFT.h>                // member

#include <pscf/crystal/Basis.h>            // member
#include <pscf/mesh/Mesh.h>                // member
#include <pscf/crystal/UnitCell.h>         // member

#include <string>

namespace Pscf {
namespace Pspc
{

   using namespace Util;

   /**
   * Spatial domain and spatial discretization for a periodic structure.
   *
   * A Domain has (among other components):
   *
   *    - a UnitCell
   *    - a Mesh
   *    - a Basis
   *    - an Fft 
   *    - a FieldIo
   *    - a groupName string
   *
   * \ingroup Pspc_Field_Module
   */
   template <int D>
   class Domain : public ParamComposite
   {

   public:

      /// \name Construction, Initialization and Destruction
      //@{

      /**
      * Constructor.
      */
      Domain();

      /**
      * Destructor.
      */
      ~Domain();

      /**
      * Create association with a FileMaster, needed by FieldIo.
      *
      * \param fileMaster associated FileMaster object.
      */
      void setFileMaster(FileMaster& fileMaster);

      /**
      * Read body of parameter block (without opening and closing lines).
      *
      * \param in input parameter stream
      */
      virtual void readParameters(std::istream& in);

      /**
      * Read initialization data from header of an r-grid field file.
      *
      * \param in input parameter stream
      * \param nMonomer number of monomers in field file (output)
      */
      void readFieldHeader(std::istream& in, int& nMonomer);

      //@}
      /// \name Accessors 
      //@{

      /**
      * Get UnitCell (i.e., lattice type and parameters) by reference.
      */
      UnitCell<D>& unitCell();

      /**
      * Get spatial discretization mesh by reference.
      */
      Mesh<D>& mesh();

      /**
      * Get associated Basis object by reference.
      */
      Basis<D>& basis();

      /**
      * Get associated FFT object.
      */
      FFT<D>& fft();

      /**
      * Get associated FieldIo object.
      */
      FieldIo<D>& fieldIo();

      /** 
      * Get group name.
      */  
      std::string groupName() const;

      //@}

   private:

      // Private member variables

      /**
      * Crystallographic unit cell (crystal system and cell parameters).
      */
      UnitCell<D> unitCell_;

      /**
      * Spatial discretization mesh.
      */
      Mesh<D> mesh_;

      /**
      * Pointer to a Basis object
      */
      Basis<D> basis_;

      /**
      * FFT object to be used by iterator
      */
      FFT<D> fft_;

      /**
      * FieldIo object for field input/output operations
      */
      FieldIo<D> fieldIo_;

      /**
      * Group name.
      */
      std::string groupName_;

      /**
      * Has a FileMaster been set?
      */
      bool hasFileMaster_;

      /**
      * Has the domain been initialized?
      */
      bool isInitialized_;

   };

   // Inline member functions

   // Get the associated UnitCell<D> object.
   template <int D>
   inline UnitCell<D>& Domain<D>::unitCell()
   {  return unitCell_; }

   // Get the Mesh<D> object.
   template <int D>
   inline Mesh<D>& Domain<D>::mesh()
   {  return mesh_; }

   // Get the Basis<D> object.
   template <int D>
   inline Basis<D>& Domain<D>::basis()
   {  return basis_; }

   // Get the FFT<D> object.
   template <int D>
   inline FFT<D>& Domain<D>::fft()
   {  return fft_; }

   // Get the FieldIo<D> object.
   template <int D>
   inline FieldIo<D>& Domain<D>::fieldIo()
   {  return fieldIo_; }

   // Get the groupName string.
   template <int D>
   inline std::string Domain<D>::groupName() const
   {  return groupName_; }

   #ifndef PSPC_DOMAIN_TPP
   // Suppress implicit instantiation
   extern template class Domain<1>;
   extern template class Domain<2>;
   extern template class Domain<3>;
   #endif

} // namespace Pspc
} // namespace Pscf
#endif
