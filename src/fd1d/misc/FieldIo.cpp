/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2019, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "FieldIo.h"

#include <fd1d/domain/Domain.h>
#include <fd1d/solvers/Mixture.h>
#include <fd1d/solvers/Polymer.h>

#include <pscf/chem/Vertex.h>
#include <pscf/homogeneous/Clump.h>
//#include <pscf/inter/Interaction.h>

#include <util/misc/FileMaster.h>
#include <util/format/Str.h>
#include <util/format/Int.h>
#include <util/format/Dbl.h>

#include <string>

namespace Pscf {
namespace Fd1d
{

   using namespace Util;

   /*
   * Constructor.
   */
   FieldIo::FieldIo()
   {}

   /*
   * Destructor.
   */
   FieldIo::~FieldIo()
   {}

 
   void FieldIo::associate(Domain const &  domain, 
                           FileMaster const &  fileMaster)
   {
      domainPtr_ = &domain;
      fileMasterPtr_ = &fileMaster;
   }

   void FieldIo::readFields(DArray<Field> &  fields, 
                            std::string const & filename)
   {
      std::ifstream in;
      fileMaster().openInputFile(filename, in);
      readFields(fields, in);
      in.close();
   }

   void FieldIo::readFields(DArray<Field>& fields, std::istream& in)
   {
      // Read grid dimensions
      std::string label;
      in >> label;
      UTIL_CHECK(label == "nx");
      int nx;
      in >> nx;
      UTIL_CHECK(nx > 0);
      UTIL_CHECK(nx == domain().nx());
      in >> label;
      UTIL_CHECK (label == "nm");
      int nm;
      in >> nm;
      UTIL_CHECK(nm > 0);

      // Check dimensions of fields array
      UTIL_CHECK(nm == fields.capacity());
      for (int i = 0; i < nm; ++i) {
         UTIL_CHECK(nx == fields[i].capacity());
      }

      // Read fields
      int i, j, idum;
      for (i = 0; i < nx; ++i) {
         in >> idum;
         UTIL_CHECK(idum == i);
         for (j = 0; j < nm; ++j) {
            in >> fields[j][i];
         }
      }
   }

   void FieldIo::writeFields(DArray<Field> const &  fields, 
                             std::string const & filename)
   {
      std::ofstream out;
      fileMaster().openOutputFile(filename, out);
      writeFields(fields, out);
      out.close();
   }

   void FieldIo::writeFields(DArray<Field> const & fields, std::ostream& out)
   {
      int nm = fields.capacity();
      UTIL_CHECK(nm > 0);
      int nx = fields[0].capacity();
      UTIL_CHECK(nx == domain().nx());
      if (nm > 1) {
         for (int i = 0; i < nm; ++i) {
            UTIL_CHECK(nx == fields[i].capacity());
         }
      }
      out << "nx     "  <<  nx              << std::endl;
      out << "nm     "  <<  nm              << std::endl;

      // Write fields
      int i, j;
      for (i = 0; i < nx; ++i) {
         out << Int(i, 5);
         for (j = 0; j < nm; ++j) {
            out << "  " << Dbl(fields[j][i], 18, 11);
         }
         out << std::endl;
      }
   }

   void FieldIo::writeBlockCFields(Mixture const & mixture, 
                                   std::string const & filename)
   {
      std::ofstream out;
      fileMaster().openOutputFile(filename, out);
      writeBlockCFields(mixture, out);
      out.close();
   }

   /*
   * Write the concentrations associated with all blocks.
   */
   void FieldIo::writeBlockCFields(Mixture const & mixture, std::ostream& out)
   {
      int nx = domain().nx();         // number grid points
      int np = mixture.nPolymer();    // number of polymer species
      int nb_tot = mixture.nBlock();  // number of blocks in all polymers
      int ns = mixture.nSolvent();    // number of solvents

      out << "nx          "  <<  nx       << std::endl;
      out << "n_block     "  <<  nb_tot   << std::endl;
      out << "n_solvent   "  <<  ns       << std::endl;

      int nb;                   // number of blocks per polymer
      int i, j, k, l;           // loop indices
      double c;
      for (i = 0; i < nx; ++i) {
         out << Int(i, 5);
         for (j = 0; j < np; ++j) {
            nb = mixture.polymer(j).nBlock();
            for (k = 0; k < nb; ++k) {
               c = mixture.polymer(j).block(k).cField()[i];
               out << " " << Dbl(c, 15, 8);
            }
         }
         for (l = 0; l < ns; ++l) {
            c = mixture.solvent(l).cField()[i];
            out << " " << Dbl(c, 15, 8);
         }
         out << std::endl;
      }
   }

   /*
   * Write incoming q fields for a specified vertex.
   */
   void FieldIo::writeVertexQ(Mixture const & mixture,
                              int polymerId, int vertexId, 
                              std::string const & filename)
   {
      std::ofstream out;
      fileMaster().openOutputFile(filename, out);
      writeVertexQ(mixture, polymerId, vertexId, out);
      out.close();
   }

   /*
   * Write incoming q fields for a specified vertex.
   */
   void FieldIo::writeVertexQ(Mixture const & mixture,
                              int polymerId, int vertexId, 
                              std::ostream& out)
   {
      UTIL_CHECK(polymerId >= 0);
      UTIL_CHECK(polymerId < mixture.nPolymer());
      Polymer const & polymer = mixture.polymer(polymerId);
      UTIL_CHECK(vertexId >= 0);
      UTIL_CHECK(vertexId <= polymer.nBlock());
      Vertex const & vertex = polymer.vertex(vertexId);
      int nb = vertex.size();   // number of attached blocks
      int nx = domain().nx();   // number grid points

      Pair<int> pId;
      int bId;
      int dId;
      int i, j;
      double c, product;
      for (i = 0; i < nx; ++i) {
         out << Int(i, 5);
         product = 1.0;
         for (j = 0; j < nb; ++j) {
            pId = vertex.inPropagatorId(j);
            bId = pId[0];  // blockId
            UTIL_CHECK(bId >= 0);
            dId = pId[1];  // directionId
            UTIL_CHECK(dId >= 0);
            UTIL_CHECK(dId <= 1);
            c = polymer.propagator(bId, dId).tail()[i];
            out << " " << Dbl(c, 15, 8);
            product *= c;
         }
         out << " " << Dbl(product, 15, 8) << std::endl;
      }
     
   }

   /*
   * Interpolate fields onto new mesh, write result to output file
   */
   void FieldIo::remesh(DArray<Field> const & fields, int nx, 
                        std::string const & filename)
   {
      std::ofstream out;
      fileMaster().openOutputFile(filename, out);
      remesh(fields, nx, out);
      out.close();
   }

   /*
   * Interpolate fields onto new mesh, write result to outputstream
   */
   void 
   FieldIo::remesh(DArray<Field> const & fields, int nx, std::ostream& out)
   {
      // Query and check dimensions of fields array
      int nm = fields.capacity();
      UTIL_CHECK(nm > 0);
      for (int i = 0; i < nm; ++i) {
         UTIL_CHECK(fields[i].capacity() == domain().nx());
      }

      // Output new grid dimensions
      out << "nx     "  <<  nx              << std::endl;
      out << "nm     "  <<  nm              << std::endl;

      // Output first grid point
      int i, j;
      i = 0;
      out << Int(i, 5);
      for (j = 0; j < nm; ++j) {
         out << "  " << Dbl(fields[j][i]);
      }
      out << std::endl;

      // Spacing for new grid
      double dx = (domain().xMax() - domain().xMin())/double(nx-1);

      // Variables used for interpolation
      double y;  // Grid coordinate in old grid
      double fu; // fractional part of y
      double fl; // 1.0 - fl
      double w;  // interpolated field value
      int    yi; // Truncated integer part of y

      // Loop over intermediate points
      for (i = 1; i < nx -1; ++i) {
         y = dx*double(i)/domain().dx();
         yi = y;
         UTIL_CHECK(yi >= 0);
         UTIL_CHECK(yi + 1 < domain().nx());
         fu = y - double(yi);
         fl = 1.0 - fu;

         out << Int(i, 5);
         for (j = 0; j < nm; ++j) {
            w = fl*fields[j][yi] + fu*fields[j][yi+1];
            out << "  " << Dbl(w);
         }
         out << std::endl;
      }

      // Output last grid point
      out << Int(nx - 1, 5);
      i = domain().nx() - 1;
      for (j = 0; j < nm; ++j) {
         out << "  " << Dbl(fields[j][i]);
      }
      out << std::endl;

   }

   void 
   FieldIo::extend(DArray<Field> const & fields, int m, 
                   std::string const & filename)
   {
      std::ofstream out;
      fileMaster().openOutputFile(filename, out);
      extend(fields, m, out);
      out.close();
   }

   /*
   * Interpolate fields onto new mesh.
   */
   void 
   FieldIo::extend(DArray<Field> const & fields, int m, std::ostream& out)
   {
      // Query and check dimensions of fields array
      int nm = fields.capacity();
      UTIL_CHECK(nm > 0);
      int nx = fields[0].capacity();
      UTIL_CHECK(nx == domain().nx());
      if (nm > 1) {
         for (int i = 0; i < nm; ++i) {
            UTIL_CHECK(nx == fields[i].capacity());
         }
      }

      // Output new grid dimensions
      out << "nx     "  <<  nx + m << std::endl;
      out << "nm     "  <<  nm              << std::endl;

      // Loop over existing points;
      int i, j;
      for (i = 0; i < nx; ++i) {
         out << Int(i, 5);
         for (j = 0; j < nm; ++j) {
            out << "  " << Dbl(fields[j][i]);
         }
         out << std::endl;
      }

      // Add m new points, assign each the last value
      for (i = nx; i < nx + m; ++i) {
         out << Int(i, 5);
         for (j = 0; j < nm; ++j) {
            out << "  " << Dbl(fields[j][nx-1]);
         }
         out << std::endl;
      }

   }

} // namespace Fd1d
} // namespace Pscf
