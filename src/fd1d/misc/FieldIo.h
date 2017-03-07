#ifndef FD1D_FIELD_IO_H
#define FD1D_FIELD_IO_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <fd1d/SystemAccess.h>          // base class
#include <util/containers/DArray.h>     // member
#include <util/containers/Array.h>      // function argument template

namespace Pscf {
namespace Fd1d {

   /**
   * Read and write fields to file.
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
      * This function uses the system FileMaster to opens and close a
      * input file named filename.
      *
      * \param filename name of input file
      * \param fields  array of fields to read, indexed by monomer id.
      */
      void readFields(std::string const & filename, 
                      Array<Field> &  fields);

      /**
      * Read a set of fields, one per monomer type.
      *
      * \pre File in must be open for reading.
      *
      * \param in  input stream, open for reading.
      * \param fields  array of fields to read, indexed by monomer id.
      */
      void readFields(std::istream &in, Array<Field> &  fields);

      /**
      * Write a set of fields, one per monomer type.
      *
      * This function uses the system FileMaster to opens and close a
      * output file named filename.
      *
      * \param filename  output filename
      * \param fields  set of fields to read.
      */
      void writeFields(std::string const & filename, 
                       Array<Field> const &  fields);

      /**
      * Write a set of fields, one per monomer type.
      *
      * \pre Stream out must be open for writing. 
      *
      * \param out  output stream 
      * \param fields  set of fields to written.
      */
      void writeFields(std::ostream &out, Array<Field> const &  fields);

      /**
      * Interpolate an array of fields onto a new mesh.
      *
      * \param filename name of output file for remeshed field
      * \param fields  field to be remeshed
      * \param nx  number of grid points in new mesh
      */
      void remesh(std::string const & filename, Array<Field> const & fields, int nx);

      /**
      * Interpolate an array of fields onto a new mesh.
      *
      * \param out  output stream for remeshed field
      * \param fields  field to be remeshed
      * \param nx  number of grid points in new mesh
      */
      void remesh(std::ostream& out, Array<Field> const & fields, int nx);

      /**
      * Add points to the end of mesh
      *
      * \param filename  name of output file for remeshed field
      * \param fields  field to be remeshed
      * \param m  number of added grid points
      */
      void extend(std::string const & filename, Array<Field> const & fields, int m);

      /**
      * Add points to the end of mesh
      *
      * \param out  output stream for extended field
      * \param fields  array of fields to be extended
      * \param m  number of added grid points
      */
      void extend(std::ostream& out, Array<Field> const & fields, int m);

   private:

      /**
      * Work array (size = # of monomer types).
      */
      DArray<double> w_;

   };

} // namespace Fd1d
} // namespace Pscf
#endif
