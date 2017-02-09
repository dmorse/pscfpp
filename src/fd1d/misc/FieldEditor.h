#ifndef FD1D_FIELD_EDITOR_H
#define FD1D_FIELD_EDITOR_H

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
   class FieldEditor : public SystemAccess
   {

   public:

      /**
      * Default constructor.
      */
      FieldEditor();

      /**
      * Constructor.
      */
      FieldEditor(System& system);

      /**
      * Destructor.
      */
      ~FieldEditor();

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
