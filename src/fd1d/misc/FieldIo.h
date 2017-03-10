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
      * \param fields  array of fields to read, indexed by monomer id.
      * \param filename  name of input file
      */
      void readFields(Array<Field> &  fields, std::string const & filename);

      /**
      * Read a set of fields, one per monomer type.
      *
      * \pre File in must be open for reading.
      *
      * \param fields  array of fields to read, indexed by monomer id
      * \param in  input stream, open for reading.
      */
      void readFields(Array<Field>& fields, std::istream &in);

      /**
      * Write a set of fields, one per monomer type.
      *
      * This function uses the system FileMaster to opens and close a
      * output file named filename.
      *
      * \param filename  output filename
      * \param fields  array of fields to read, indexed by monomer id
      */
      void writeFields(Array<Field> const &  fields, 
                       std::string const & filename);

      /**
      * Write a set of fields, one per monomer type.
      *
      * \pre Stream out must be open for writing. 
      *
      * \param fields  set of fields to written.
      * \param out  output stream 
      */
      void writeFields(Array<Field> const & fields, std::ostream& out);

      /**
      * Write block concentration fields for all blocks.
      *
      * \param filename name of output file
      */
      void writeBlockCFields(std::string const & out);

      /**
      * Write block concentration fields for all blocks.
      *
      * \pre Stream out must be open for writing. 
      *
      * \param out  output stream 
      */
      void writeBlockCFields(std::ostream& out);

      /**
      * Write incoming q fields for a specified vertex.
      *
      * \param polymerId integer id of polymer species
      * \param vertexId integer id of vertex (end or junction)
      * \param filename name of output file
      */
      void writeVertexQ(int polymerId, int vertexId, std::string const & out);

      /**
      * Write incoming q fields for a specified vertex.
      *
      * \pre Stream out must be open for writing. 
      *
      * \param polymerId  integer id of polymer species
      * \param vertexId  integer id of vertex (end or junction)
      * \param out  output stream 
      */
      void writeVertexQ(int polymerId, int vertexId, std::ostream& out);

      /**
      * Interpolate an array of fields onto a new mesh and write to file.
      *
      * \param fields  field to be remeshed
      * \param nx  number of grid points in new mesh
      * \param filename name of output file for remeshed field
      */
      void remesh(Array<Field> const & fields, int nx,
                  std::string const & filename);

      /**
      * Interpolate an array of fields onto a new mesh.
      *
      * \param out  output stream for remeshed field
      * \param fields  field to be remeshed
      * \param nx  number of grid points in new mesh
      */
      void remesh(Array<Field> const & fields, int nx, std::ostream& out);

      /**
      * Add points to the end of mesh
      *
      * \param filename  name of output file for remeshed field
      * \param fields  field to be remeshed
      * \param m  number of added grid points
      */
      void extend(Array<Field> const & fields, int m, 
                  std::string const & filename);

      /**
      * Add points to the end of mesh
      *
      * \param out  output stream for extended field
      * \param fields  array of fields to be extended
      * \param m  number of added grid points
      */
      void extend(Array<Field> const & fields, int m, std::ostream& out);

   private:

      /**
      * Work array (size = # of monomer types).
      */
      DArray<double> w_;

   };

} // namespace Fd1d
} // namespace Pscf
#endif
