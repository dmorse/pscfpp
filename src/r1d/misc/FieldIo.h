#ifndef R1D_FIELD_IO_H
#define R1D_FIELD_IO_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/containers/DArray.h>     // member

namespace Util { class FileMaster; }

namespace Pscf {
namespace R1d {

   class Domain;
   class Mixture;
   using namespace Util;

   /**
   * Read and write fields to file.
   *
   * \ingroup Pscf_R1d_Module
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
      * This function uses the associated FileMaster to open an input
      * file named filename before reading, and closes the file after
      * reading.
      *
      * \param fields  array of fields to read, indexed by monomer id.
      * \param filename  name of input file
      */
      void readFields(DArray<Field>& fields, std::string const& filename);
      
      /**
      * Write a single field to an output stream.
      *
      * \pre Stream out must be open for writing. 
      *
      * \param field  field defined on r-space grid (input)
      * \param out  output stream
      * \param writeHeader  write file header iff this bool is true
      */
      void writeField(Field const& field, std::ostream& out, 
                      bool writeHeader = true) const;

      /**
      * Write a single field to a file.
      *
      * This function uses the associated FileMaster to open an output
      * file named filename before writing, and closes the file after
      * writing.
      *
      * \param field  field defined on r-space grid (input)
      * \param filename  output filename
      * \param writeHeader  write file header iff this bool is true
      */
      void writeField(Field const&  field, std::string const& filename, 
                      bool writeHeader= true) const;
      
      /**
      * Write a set of fields, one per monomer type, to an output stream.
      *
      * \pre Stream out must be open for writing. 
      *
      * \param fields  set of fields to written.
      * \param out  output stream 
      * \param writeHeader  write file header iff this bool is true
      */
      void writeFields(DArray<Field> const& fields, std::ostream& out, 
                       bool writeHeader= true);

      /**
      * Write a set of fields, one per monomer type, to a named file.
      *
      * This function uses the associated FileMaster to open an output
      * file named filename before writing, and closes the file after
      * writing.
      *
      * \param fields  array of fields to read, indexed by monomer id
      * \param filename  output filename
      * \param writeHeader  write header iff this bool is true
      */
      void writeFields(DArray<Field> const&  fields, 
                       std::string const& filename, bool writeHeader= true);

      /**
      * Write block concentration fields for all blocks to an output stream.
      *
      * \pre Stream out must be open for writing. 
      *
      * \param mixture  associated Mixture MDE solver object
      * \param out  output stream 
      */
      void writeBlockCFields(Mixture const& mixture, std::ostream& out);

      /**
      * Write block concentration fields for all blocks to a named file.
      *
      * This function uses the associated FileMaster to open an output
      * file named filename before writing, and closes the file after
      * writing.
      *
      *
      * \param mixture  associated Mixture MDE solver object
      * \param filename name of output file
      */
      void writeBlockCFields(Mixture const& mixture,
                             std::string const& filename);

      /**
      * Write product of incoming q fields for one vertex to stream.
      *
      * \pre Stream out must be open for writing. 
      *
      * \param mixture  associated Mixture MDE solver object
      * \param polymerId  integer id of polymer species
      * \param vertexId  integer id of vertex (end or junction)
      * \param out  output stream 
      */
      void writeVertexQ(Mixture const& mixture,
                        int polymerId, int vertexId, std::ostream& out);

      /**
      * Write incoming q fields for a specified vertex.
      *
      * \param mixture  associated Mixture MDE solver object
      * \param polymerId integer id of polymer species
      * \param vertexId integer id of vertex (end or junction)
      * \param filename name of output file
      */
      void writeVertexQ(Mixture const& mixture,
                        int polymerId, int vertexId, 
                        std::string const& filename);

      /**
      * Interpolate an array of fields onto a new mesh and write to stream.
      *
      * \param fields  field to be remeshed
      * \param nx  number of grid points in new mesh
      * \param out  output stream for remeshed field
      */
      void remesh(DArray<Field> const& fields, int nx, std::ostream& out);

      /**
      * Interpolate an array of fields onto a new mesh and write to file.
      *
      * \param fields  field to be remeshed
      * \param nx  number of grid points in new mesh
      * \param filename name of output file for remeshed field
      */
      void remesh(DArray<Field> const& fields, int nx,
                  std::string const& filename);

      /**
      * Add points to the end of a field mesh and write to stream.
      *
      * Values at the added mesh points are taken to be the same as
      * those at the last mesh point of the original mesh. 
      *
      * \param fields  array of fields to be extended
      * \param m  number of added grid points
      * \param out  output stream for extended field
      */
      void extend(DArray<Field> const& fields, int m, std::ostream& out);

      /**
      * Add points to the end of a field mesh and write to a file.
      *
      * Values at the added mesh points are taken to be the same as
      * those at the last mesh point of the original mesh. 
      *
      * \param fields  field to be remeshed
      * \param m  number of added grid points
      * \param filename  name of output file for remeshed field
      */
      void extend(DArray<Field> const& fields, int m, 
                  std::string const& filename);

   private:

      /// Work array (capacity = # of monomer types).
      mutable DArray<double> w_;

      // Pointers to associated objects.

      /// Pointer to spatial discretization domain.
      Domain const * domainPtr_;

      /// Pointer to Filemaster (holds paths to associated I/O files).
      FileMaster const * fileMasterPtr_;

      // Private accessor functions:

      /// Get spatial discretization domain by const reference.
      Domain const& domain() const
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


} // namespace R1d
} // namespace Pscf
#endif
