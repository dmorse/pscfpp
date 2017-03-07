#ifndef FD1D_FIELD_IO_H
#define FD1D_FIELD_IO_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <fd1d/SystemAccess.h>             // base class

namespace Pscf {
namespace Fd1d {

   /**
   * 
   *
   * \ingroup Pscf_Fd1d_Module
   */
   class FieldIo : public SystemAccess
   {

   public:

      typedef System::Field Field;

      /**
      * Default constructor.
      */
      FieldIo();

      /**
      * Constructor.
      */
      FieldIo(System& system);

      /**
      * Destructor.
      */
      ~FieldIo();


      /**
      * Read a set of fields, one per monomer type.
      *
      * \param fields  set of fields to read.
      * \param in  input stream 
      */
      void readFields(std::string const & filename, 
                      Array<Field> &  fields);

      /**
      * Read a set of fields, one per monomer type.
      *
      * \param fields  set of fields to read.
      * \param in  input stream 
      */
      void readFields(std::istream &in, Array<Field> &  fields);

      /**
      * Write a set of fields, one per monomer type.
      *
      * \param filename  output filename
      * \param fields  set of fields to read.
      */
      void writeFields(std::string const & filename, 
                       Array<Field> const &  fields);

      /**
      * Write a set of fields, one per monomer type.
      *
      * \param fields set of fields to written.
      * \param out  output stream 
      */
      void writeFields(std::ostream &out, Array<Field> const &  fields);

      /**
      * Interpolate a set of fields onto a new mesh.
      *
      * \param fields set of fields to be remeshed
      * \param nx number of grid points in new mesh
      * \param out output stream for remeshed field
      */
      void writeFields(DArray<System::WField>& fields, std::ostream& out);

      /**
      * Interpolate field onto a new mesh.
      *
      * \param field field to be remeshed
      * \param nx    number of grid points in new mesh
      * \param out   output stream for remeshed field
      */
      void remesh(DArray<System::WField>& fields, int nx, std::ostream& out);

      /**
      * Add points to the end of 
      *
      * \param field field to be remeshed
      * \param m     number of added grid points 
      * \param out   output stream for remeshed field
      */
      void extend(DArray<System::WField>& fields, int m, std::ostream& out);

   private:

      /**
      * Work array (size = # of monomer types).
      */
      DArray<double> w_;

   };

} // namespace Fd1d
} // namespace Pscf
#endif
