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

 
   void FieldIo::readFields(std::string const & filename, 
                            Array<Field> &  fields)
   {
      std::ifstream in;
      fileMaster().openInputFile(filename, in);
      readFields(in, fields);
      in.close();
   }

   void FieldIo::readFields(std::istream &in, Array<Field> &  fields)
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

   void FieldIo::writeFields(std::string const & filename, 
                            Array<Field> const &  fields)
   {
      std::ofstream out;
      fileMaster().openOutputFile(filename, out);
      writeFields(out, fields);
      out.close();
   }

   void FieldIo::writeFields(std::ostream &out, 
                             Array<Field> const &  fields)
   {
      int i, j;
      int nx = domain().nx();
      int nm = mixture().nMonomer();
      out << "nx     "  <<  nx              << std::endl;
      out << "nm     "  <<  nm              << std::endl;

      // Write fields
      for (i = 0; i < nx; ++i) {
         out << Int(i, 5);
         for (j = 0; j < nm; ++j) {
            out << "  " << Dbl(fields[j][i], 18, 11);
         }
         out << std::endl;
      }
   }

   void FieldIo::remesh(std::string const & filename, 
                        Array<Field> const &  fields, int nx)
   {
      std::ofstream out;
      fileMaster().openOutputFile(filename, out);
      remesh(out, fields, nx);
      out.close();
   }

   /*
   * Interpolate fields onto new mesh.
   */
   void 
   FieldIo::remesh(std::ostream& out, Array<Field> const & fields, int nx)
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
   FieldIo::extend(std::string const & filename, 
                   Array<Field> const &  fields, int m)
   {
      std::ofstream out;
      fileMaster().openOutputFile(filename, out);
      extend(out, fields, m);
      out.close();
   }

   /*
   * Interpolate fields onto new mesh.
   */
   void 
   FieldIo::extend(std::ostream& out, Array<Field> const & fields, int m)
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
