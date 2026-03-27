#ifndef PRDC_CL_W_FIELDS_H
#define PRDC_CL_W_FIELDS_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/containers/DArray.h>      // member template
#include <pscf/math/IntVec.h>            // member

// Forward declarations
namespace Util {
   template <typename T> class Signal;
   template <> class Signal<void>;
}
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
   * A container of complex-valued w fields.
   * 
   * <b> Overview </b>: A WFields container has an array of chemical
   * fields (w fields) that are each associated with a monomer type. 
   * Fields can be modified by the setFields and readFields member 
   * functions, accessed as const references, or written to file.
   *
   * <b> Template parameters </b>: The template parameters represent:
   * 
   *     - D   : integer dimensionality of space, D=1,2, or 3
   *     - CFT : field type for complex data (e.g., CField<D>)
   *     - FIT : FieldIo type for field io operations (e.g., FieldIo<D>)
   * 
   * <b> Subclasses </b>: Partial specializations of WFields are used as
   * base classes for classes Cpc::WFields \<D \> and Rpg::WFields \<D\>:
   *
   *  - Subclass Cpc::WFields \<D\> is derived from a specialization of 
   *    this template with template arguments CFT = Cpu::CField \<D\> 
   *    and FIT = Cpc::FieldIo \<D\> , and is used in the pscf_cpc CPU 
   *    CL-FTS program.
   *
   *  - Subclass Rpg::WFields \<D\> is derived from an instantation of
   *    this template with template arguments CFT = Cuda::CField \<D\> 
   *    and FIT = Rpg::FieldIo \<D\> , and is used in the pscf_cpg GPU 
   *    accelerated CL-FTS program.
   *
   * <b> Signal </b>: A WFields owns an instance of class 
   * Util::Signal<void> that notifies all observers whenever the fields
   * owned by the WFields are modified by the set or readField function. 
   * This signal object may be accessed by reference using the signal() 
   * member function. The Util::Signal<void>::addObserver function may 
   * used to add an "observer" object and indicate a zero-parameter member 
   * function of each observer that will be called whenever the w fields 
   * are modified.
   *
   * \ingroup Cp_Field_Module
   */
   template <int D, class CFT, class FIT>
   class WFields 
   {

   public:

      /**
      * Constructor.
      */
      WFields();

      /**
      * Destructor.
      */
      ~WFields();

      /// \name Initialization and Memory Management
      ///@{

      /**
      * Create association with FIT (store pointer).
      *
      * \param fieldIo  associated FIT object
      */
      void setFieldIo(FIT const & fieldIo);

      /**
      * Set unit cell used when reading field files. 
      *
      * This function creates a stored pointer to a UnitCell<D> that is
      * is used by the readFields functions, which reset the unit cell
      * parameters in this object to those read from the field file
      * header. This function may only be called once.
      *
      * \param cell  unit cell that is modified by readFields 
      */
      void setReadUnitCell(UnitCell<D>& cell);

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
      /// \name Field and State Modifiers
      ///@{

      /**
      * Set values for all fields.
      *
      * On return, hasData is true.
      *
      * \param fields  array of new fields 
      */
      void setFields(DArray< CFT > const & fields);

      /**
      * Read all fields from an input file.
      * 
      * \param in  input stream from which to read fields
      */
      void readFields(std::istream& in);

      /**
      * Read all fields from a named file.
      *
      * \param filename  file from which to read fields
      */
      void readFields(std::string const & filename);

      /**
      * Clear data stored in this object without deallocating.
      */
      void clear();

      /**
      * Get the signal that notifies observers of w-field modification.
      */
      Signal<void>& signal();

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
      /// \name Field Access (by const reference)
      ///@{

      /**
      * Get the array of all fields.
      *
      * The array capacity is equal to the number of monomer types.
      */
      DArray< CFT > const & fields() const;

      /**
      * Get the field for one monomer type.
      *
      * \param monomerId integer monomer type index (0,..,nMonomer-1)
      */
      CFT const & field(int monomerId) const;

      ///@}
      /// \name Boolean Queries
      ///@{

      /**
      * Has memory been allocated for fields ?
      */
      bool isAllocated() const;

      /**
      * Has field data been set since it was last cleared?
      *
      * This flag is set true in the setFields and readFields functions,
      * and set false in the clear function.
      */
      bool hasData() const;

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
      * Pointer to unit cell modified by the read functions.
      */
      UnitCell<D> * readUnitCellPtr_;

      /*
      * Pointer to unit cell written by the write functions.
      */
      UnitCell<D> const * writeUnitCellPtr_;

      /*
      * Pointer to an associated FIT (i.e., FieldIo<D>) object.
      */
      FIT const * fieldIoPtr_;

      /*
      * Pointer to a Signal that is triggered by field modification.
      *
      * The Signal is constructed and owned by this container.
      */
      Signal<void>* signalPtr_;
 
      /*
      * Has memory been allocated for the fields? 
      */
      bool isAllocated_;

      /*
      * Has field data been set since it was last cleared?
      */
      bool hasData_;

      /*
      *  Assign one complex field (CFT) to another: lhs = rhs.
      *  
      *  \param lhs  left hand side of assignment
      *  \param rhs  right hand side of assignment
      */
      virtual void assignField(CFT& lhs, CFT const & rhs) const;

   };

   // Inline member functions

   // Mark field data stored in this object as invalid or outdated.
   template <int D, class CFT, class FIT> inline 
   void WFields<D,CFT,FIT>::clear()
   {  hasData_ = false; }

   // Get the array of all fields (const reference).
   template <int D, class CFT, class FIT> inline
   DArray< CFT > const & WFields<D,CFT,FIT>::fields() const
   {
      UTIL_ASSERT(isAllocated_);
      return fields_;
   }

   // Get a single field (const reference).
   template <int D, class CFT, class FIT> inline
   CFT const & WFields<D,CFT,FIT>::field(int id) const
   {
      UTIL_ASSERT(isAllocated_);
      return fields_[id];
   }

   // Has memory been allocated for fields ?
   template <int D, class CFT, class FIT> inline 
   bool WFields<D,CFT,FIT>::isAllocated() const
   {  return isAllocated_; }

   // Has field data been initialized ?
   template <int D, class CFT, class FIT> inline 
   bool WFields<D,CFT,FIT>::hasData() const
   {  return hasData_; }

   // Protected inline member functions
   
   // Get mesh dimensions in each direction.
   template <int D, class CFT, class FIT> inline 
   IntVec<D> const & 
   WFields<D,CFT,FIT>::meshDimensions() const
   {  return meshDimensions_; }

   // Get mesh size (number of grid points).
   template <int D, class CFT, class FIT> inline 
   int WFields<D,CFT,FIT>::meshSize() const
   {  return meshSize_; }

   // Get number of monomer types.
   template <int D, class CFT, class FIT> inline 
   int WFields<D,CFT,FIT>::nMonomer() const
   {  return nMonomer_; }

   // Get the associated FIT (FieldIo<D>) object (const reference).
   template <int D, class CFT, class FIT> inline 
   FIT const & WFields<D,CFT,FIT>::fieldIo() const
   {
      UTIL_CHECK(fieldIoPtr_);
      return *fieldIoPtr_;
   }

} // namespace Cp
} // namespace Pscf
#endif
