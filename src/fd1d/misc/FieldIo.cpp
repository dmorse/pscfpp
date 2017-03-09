/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "FieldIo.h"

#include <pscf/inter/Interaction.h>
#include <pscf/inter/ChiInteraction.h>
#include <pscf/homogeneous/Clump.h>

#include <util/format/Str.h>
#include <util/format/Int.h>
#include <util/format/Dbl.h>

#include <string>

namespace Pscf {
namespace Fd1d
{

   using namespace Util;

   /*
   * Default constructor.
   */
   FieldIo::FieldIo()
    : SystemAccess()
   {}

   /*
   * Constructor.
   */
   FieldIo::FieldIo(System& system)
    : SystemAccess(system)
   {}

   /*
   * Destructor.
   */
   FieldIo::~FieldIo()
   {}

 
   void FieldIo::readFields(Array<Field> &  fields, 
                            std::string const & filename)
   {
      std::ifstream in;
      fileMaster().openInputFile(filename, in);
      readFields(fields, in);
      in.close();
   }

   void FieldIo::readFields(Array<Field>& fields, std::istream& in)
   {
      // Read grid dimensions
      std::string label;
      int nx, nm;
      in >> label;
      UTIL_CHECK(label == "nx");
      in >> nx;
      UTIL_CHECK(nx > 0);
      UTIL_CHECK(nx == domain().nx());
      in >> label;
      UTIL_CHECK (label == "nm");
      in >> nm;
      UTIL_CHECK(nm > 0);
      UTIL_CHECK(nm == mixture().nMonomer());

      // Read fields
      int i,j, idum;
      for (i = 0; i < nx; ++i) {
         in >> idum;
         UTIL_CHECK(idum == i);
         for (j = 0; j < nm; ++j) {
            in >> fields[j][i];
         }
      }
   }

   void FieldIo::writeFields(Array<Field> const &  fields, 
                             std::string const & filename)
   {
      std::ofstream out;
      fileMaster().openOutputFile(filename, out);
      writeFields(fields, out);
      out.close();
   }

   void FieldIo::writeFields(Array<Field> const & fields, std::ostream& out)
   {
      int nx = domain().nx();
      int nm = mixture().nMonomer();
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

   void FieldIo::writeBlockCFields(std::string const & filename)
   {
      std::ofstream out;
      fileMaster().openOutputFile(filename, out);
      writeBlockCFields(out);
      out.close();
   }

   /*
   * Write the concentrations associated with all blocks.
   */
   void FieldIo::writeBlockCFields(std::ostream& out)
   {
      int nx = domain().nx();          // number grid points
      int np = mixture().nPolymer();   // number of polymer species
      int nb;                          // number of blocks per polymer
      int i, j, k;
      double c;
      for (i = 0; i < nx; ++i) {
         out << Int(i, 5);
         for (j = 0; j < np; ++j) {
            nb = mixture().polymer(j).nBlock();
            for (k = 0; k < nb; ++k) {
               c = mixture().polymer(j).block(k).cField()[i];
               out << " " << Dbl(c, 15, 8);
            }
         }
         out << std::endl;
      }
   }

   void FieldIo::remesh(Array<Field> const &  fields, int nx, 
                        std::string const & filename)
   {
      std::ofstream out;
      fileMaster().openOutputFile(filename, out);
      remesh(fields, nx, out);
      out.close();
   }

   /*
   * Interpolate fields onto new mesh.
   */
   void 
   FieldIo::remesh(Array<Field> const & fields, int nx, std::ostream& out)
   {
      int nm = mixture().nMonomer();
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

      double y;  // Grid coordinate in old grid
      double fu; // fractional part of y
      double fl; // 1.0 - fl
      double w;  // interpolated field value
      int    yi; // Truncated integer part of y

      // Loop over intermediate points
      for (i = 1; i < nx -1; ++i) {
         y = dx*double(i)/domain().dx();
         yi = y;
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
   FieldIo::extend(Array<Field> const & fields, int m, 
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
   FieldIo::extend(Array<Field> const & fields, int m, std::ostream& out)
   {
      int nm = mixture().nMonomer();
      int nx = domain().nx();

      out << "nx     "  <<  nx + m << std::endl;
      out << "nm     "  <<  nm              << std::endl;
      int i, j;

      // Loop over existing points;
      for (i = 0; i < nx; ++i) {
         out << Int(i, 5);
         for (j = 0; j < nm; ++j) {
            out << "  " << Dbl(fields[j][i]);
         }
         out << std::endl;
      }

      // Add m new points
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
