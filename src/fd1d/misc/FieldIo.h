#ifndef FD1D_FIELD_IO_H
#define FD1D_FIELD_IO_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2019, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/containers/DArray.h>     // member

namespace Util { class FileMaster; }

namespace Pscf {
namespace Fd1d {

   class Domain;
   class Mixture;
   using namespace Util;

   /**
   * Read and write fields to file.
   *
   * \ingroup Pscf_Fd1d_Module
   */
   class FieldIo 
   {

   public:

      typedef DArray<double> Field;

      /**
      * Constructor.
      */
      FieldIo();

      /**
      * Destructor.
      */
      ~FieldIo();

      /**
      * Get and store addresses of associated objects.
      *
      * \param domain  associated spatial domain
      * \param fileMaster  associated FileMaster (for file paths)
      */
      void associate(Domain const & domain,
                     FileMaster const & fileMaster);

      /**
      * Read a set of fields, one per monomer type.
      *
      * \pre File in must be open for reading.
      *
      * \param fields  array of fields to read, indexed by monomer id
      * \param in  input stream, open for reading.
      */
      void readFields(DArray<Field>& fields, std::istream &in);

      /**
      * Read a set of fields, one per monomer type.
      *
      * This function uses the system FileMaster to opens and close a
      * input file named filename.
      *
      * \param fields  array of fields to read, indexed by monomer id.
      * \param filename  name of input file
      */
      void readFields(DArray<Field> &  fields, std::string const & filename);

      /**
      * Write a set of fields, one per monomer type.
      *
      * \pre Stream out must be open for writing. 
      *
      * \param fields  set of fields to written.
      * \param out  output stream 
      */
      void writeFields(DArray<Field> const & fields, std::ostream& out);

      /**
      * Write a set of fields, one per monomer type.
      *
      * This function uses the system FileMaster to opens and close a
      * output file named filename.
      *
      * \param fields  array of fields to read, indexed by monomer id
      * \param filename  output filename
      */
      void writeFields(DArray<Field> const &  fields, 
                       std::string const & filename);

      /**
      * Write block concentration fields for all blocks.
      *
      * \pre Stream out must be open for writing. 
      *
      * \param mixture  associated Mixture MDE solver object
      * \param out  output stream 
      */
      void writeBlockCFields(Mixture const & mixture, std::ostream& out);

      /**
      * Write block concentration fields for all blocks.
      *
      * This function uses the system FileMaster to opens and close a
      * output file named filename.
      *
      * \param mixture  associated Mixture MDE solver object
      * \param filename name of output file
      */
      void writeBlockCFields(Mixture const & mixture,
                             std::string const & filename);

      /**
      * Write incoming q fields for a specified vertex.
      *
      * \pre Stream out must be open for writing. 
      *
      * \param mixture  associated Mixture MDE solver object
      * \param polymerId  integer id of polymer species
      * \param vertexId  integer id of vertex (end or junction)
      * \param out  output stream 
      */
      void writeVertexQ(Mixture const & mixture,
                        int polymerId, int vertexId, std::ostream& out);

      /**
      * Write incoming q fields for a specified vertex.
      *
      * \param mixture  associated Mixture MDE solver object
      * \param polymerId integer id of polymer species
      * \param vertexId integer id of vertex (end or junction)
      * \param filename name of output file
      */
      void writeVertexQ(Mixture const & mixture,
                        int polymerId, int vertexId, 
                        std::string const & filename);

      /**
      * Interpolate an array of fields onto a new mesh.
      *
      * \param fields  field to be remeshed
      * \param nx  number of grid points in new mesh
      * \param out  output stream for remeshed field
      */
      void remesh(DArray<Field> const & fields, int nx, std::ostream& out);

      /**
      * Interpolate an array of fields onto a new mesh and write to file.
      *
      * \param fields  field to be remeshed
      * \param nx  number of grid points in new mesh
      * \param filename name of output file for remeshed field
      */
      void remesh(DArray<Field> const & fields, int nx,
                  std::string const & filename);

      /**
      * Add points to the end of mesh
      *
      * \param fields  array of fields to be extended
      * \param m  number of added grid points
      * \param out  output stream for extended field
      */
      void extend(DArray<Field> const & fields, int m, std::ostream& out);

      /**
      * Add points to the end of mesh
      *
      * \param fields  field to be remeshed
      * \param m  number of added grid points
      * \param filename  name of output file for remeshed field
      */
      void extend(DArray<Field> const & fields, int m, 
                  std::string const & filename);

   private:

      /// Work array (capacity = # of monomer types).
      DArray<double> w_;

      // Pointers to associated objects.

      /// Pointer to spatial discretization domain.
      Domain const * domainPtr_;

      /// Pointer to Filemaster (holds paths to associated I/O files).
      FileMaster const * fileMasterPtr_;

      // Private accessor functions:

      /// Get spatial discretization domain by const reference.
      Domain const & domain() const
      {  
         UTIL_ASSERT(domainPtr_);  
         return *domainPtr_; 
      }

      /// Get FileMaster by reference.
      FileMaster const & fileMaster() const
      {  
         UTIL_ASSERT(fileMasterPtr_);  
         return *fileMasterPtr_; 
      }

   };


} // namespace Fd1d
} // namespace Pscf
#endif
