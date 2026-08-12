#ifndef PRDC_CL_C_FIELDS_H
#define PRDC_CL_C_FIELDS_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/containers/DArray.h>      // member template
#include <pscf/math/IntVec.h>            // member

// Forward declarations
namespace Pscf {
   namespace Prdc {
      template <int D> class UnitCell;
   }
}

namespace Pscf {
namespace Cp {

   using namespace Util;
   using namespace Prdc;

   /**
   * A container of complex-valued concentration fields (c fields).
   * 
   * <b> Overview </b>: A CFields container has an array of nMonomer chemical
   * fields (c fields) that are each associated with a monomer type. Fields 
   * can be accessed as const or non const references, or written to file.
   *
   * <b> Template parameters </b>: The template parameters represent:
   * 
   *     - D   : integer dimensionality of space, D=1,2, or 3
   *     - CFT : complex field type (e.g., CField<D>)
   *     - FIT : FieldIo type for field io operations (e.g., FieldIo<D>)
   * 
   * <b> Subclasses </b>: Partial specializations of CFields are used as
   * base classes for classes Cpc::CFields \<D \> and Rpg::CFields \<D\>:
   *
   *  - Subclass Cpc::CFields \<D\> is derived from a partial
   *    specialization of CFields with template parameters 
   *    CFT = Cpu::CFT \<D\> and FIT = Cpc::FIT \<D\> , and is used in
   *    the pscf_cpc CPU program.
   *
   *  - Subclass Rpg::CFields \<D\> is derived from a partial
   *    specialization of CFields with template parameters 
   *    CFT = Cuda::CFT \<D\> and FIT = Rpg::FIT \<D\> , and is used in
   *    the pscf_cpg GPU accelerated program.
   *
   *
   * \ingroup Cp_Field_Module
   */
   template <int D, class CFT, class FIT>
   class CFields 
   {

   public:

      /**
      * Constructor.
      */
      CFields();

      /**
      * Destructor.
      */
      ~CFields();

      /// \name Initialization and Memory Management
      ///@{

      /**
      * Create association with FIT (store pointer).
      *
      * \param fieldIo  associated FIT object
      */
      void setFieldIo(FIT const & fieldIo);

      /**
      * Set unit cell used when writing field files.
      *
      * This function creates a stored pointer to a UnitCell<D> that is
      * is used by the writeFields function, which write the unit cell
      * parameters from in this object to a field file header. This 
      * function may only be called once.
      *
      * \param cell  unit cell that is used by writeFields
      */
      void setWriteUnitCell(UnitCell<D> const & cell);

      /**
      * Allocate memory for fields.
      *
      * This function may only be called once.
      *
      * \param nMonomer  number of monomer types
      * \param dimensions  dimensions of spatial mesh
      */
      void allocate(int nMonomer, IntVec<D> const & dimensions);

      ///@}
      /// \name Field Output to File
      ///@{

      /**
      * Write fields to an output stream.
      *
      * \param out  output stream 
      */
      void writeFields(std::ostream& out) const;

      /**
      * Write fields to a named file.
      *  
      * \param filename  name of output file
      */
      void writeFields(std::string const & filename) const;

      ///@}
      /// \name Field Access (by reference)
      ///@{

      /**
      * Get the array of all fields (non-const reference).
      *
      * The array capacity is equal to the number of monomer types.
      */
      DArray< CFT > & fields();

      /**
      * Get the array of all fields (const reference).
      *
      * The array capacity is equal to the number of monomer types.
      */
      DArray< CFT > const & fields() const;

      /**
      * Get the field for one monomer type (non-const reference).
      *
      * \param monomerId integer monomer type index (0,..,nMonomer-1)
      */
      CFT & field(int monomerId);

      /**
      * Get the field for one monomer type (const reference).
      *
      * \param monomerId integer monomer type index (0,..,nMonomer-1)
      */
      CFT const & field(int monomerId) const;

      ///@}
      /// \name Boolean Flags
      ///@{

      /**
      * Has memory been allocated for fields ?
      */
      bool isAllocated() const;

      /**
      * Does this container have up-to-date field data ?
      */
      bool hasData() const;

      /**
      * Set the hasData boolean flag.
      *
      * This should be set true when fields are set to those computed
      * from the current w fields, and set false when any input to that
      * computation changes.
      */
      void setHasData(bool hasData);

      ///@}

   protected:

      /**
      * Get mesh dimensions in each direction, set on r-grid allocation.
      */
      IntVec<D> const & meshDimensions() const;

      /**
      * Get mesh size (number of grid points), set on r-grid allocation.
      */
      int meshSize() const;

      /**
      * Get number of monomer types.
      */
      int nMonomer() const;

      /**
      * Get associated FIT field IO object (const reference).
      */
      FIT const & fieldIo() const;

   private:

      /*
      * Array of complex valued fields.
      *
      * Element fields_[i] is a CFT (complex field type) object that 
      * contains values of the field associated with monomer i on the 
      * nodes of a regular mesh.
      */
      DArray< CFT > fields_;

      /*
      * Integer vector of grid dimensions.
      *
      * Element i is the number of grid points along direction i.
      */
      IntVec<D> meshDimensions_;

      /*
      * Total number grid points (product of all mesh dimensions).
      */
      int meshSize_;

      /*
      * Number of monomer types (number of fields).
      */
      int nMonomer_;

      /*
      * Pointer to unit cell written by the write functions.
      */
      UnitCell<D> const * writeUnitCellPtr_;

      /*
      * Pointer to an associated FIT (field IO) object.
      */
      FIT const * fieldIoPtr_;

      /*
      * Has memory been allocated for the fields? 
      */
      bool isAllocated_;

      /*
      * Does this container hold up-to-date field data? 
      */
      bool hasData_;

   };

   // Inline member functions

   // Get the array of all fields (const reference)
   template <int D, class CFT, class FIT> inline
   DArray< CFT > const & CFields<D,CFT,FIT>::fields() const
   {
      UTIL_ASSERT(isAllocated_);
      return fields_;
   }

   // Get the array of all fields (non-const reference)
   template <int D, class CFT, class FIT> inline
   DArray< CFT > & CFields<D,CFT,FIT>::fields() 
   {
      UTIL_ASSERT(isAllocated_);
      return fields_;
   }

   // Get a single field (const reference)
   template <int D, class CFT, class FIT> inline
   CFT const & CFields<D,CFT,FIT>::field(int id) const
   {
      UTIL_ASSERT(isAllocated_);
      return fields_[id];
   }

   // Get a single field (non-const reference)
   template <int D, class CFT, class FIT> inline
   CFT & CFields<D,CFT,FIT>::field(int id)
   {
      UTIL_ASSERT(isAllocated_);
      return fields_[id];
   }

   // Has memory been allocated for fields ?
   template <int D, class CFT, class FIT> inline 
   bool CFields<D,CFT,FIT>::isAllocated() const
   {  return isAllocated_; }

   // Is data up to date ?
   template <int D, class CFT, class FIT> inline 
   bool CFields<D,CFT,FIT>::hasData() const
   {  return hasData_; }

   // Protected inline member functions
   
   // Get mesh dimensions in each direction.
   template <int D, class CFT, class FIT> inline 
   IntVec<D> const & 
   CFields<D,CFT,FIT>::meshDimensions() const
   {  return meshDimensions_; }

   // Get mesh size (number of grid points).
   template <int D, class CFT, class FIT> inline 
   int CFields<D,CFT,FIT>::meshSize() const
   {  return meshSize_; }

   // Get number of monomer types.
   template <int D, class CFT, class FIT> inline 
   int CFields<D,CFT,FIT>::nMonomer() const
   {  return nMonomer_; }

   // Get the associated the field Io (FIT) object (const reference).
   template <int D, class CFT, class FIT> inline 
   FIT const & CFields<D,CFT,FIT>::fieldIo() const
   {
      UTIL_CHECK(fieldIoPtr_);
      return *fieldIoPtr_;
   }

} // namespace Cp
} // namespace Pscf
#endif
