#ifndef PSPC_BASIS_FIELD_STATE_H
#define PSPC_BASIS_FIELD_STATE_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2021, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "FieldState.h"
#include <string>

namespace Pscf {
namespace Pspc
{

   using namespace Util;

   /**
   * FieldState for fields in symmetry-adapted basis format.
   */
   template <int D>
   class BasisFieldState : public FieldState<D, DArray<double> >
   {

      /**
      * Default constructor.
      */
      BasisFieldState();

      /**
      * Constructor, associate with system.
      */
      BasisFieldState(System<D>&);

      /**
      * Destructor.
      */
      ~BasisFieldState();

      /**
      * Read state from file.
      *
      * \param filename name of input w-field basis file
      */
      void read(std::string & filename);
   
      /**
      * Write state to file.
      *
      * \param filename name of input w-field r-grid file
      */
      void write(std::string & filename);

      /**
      * Get the current state of the associated system.
      */
      void getSystemState();

      /**
      * Set the state of the associated system to this state.
      */
      void setSystemState();

      using FieldState<D, DArray<double> >::fields;
      using FieldState<D, DArray<double> >::field;
      using FieldState<D, DArray<double> >::unitCell;
      using FieldState<D, DArray<double> >::fieldIo;
      using FieldState<D, DArray<double> >::system;

   private:

      void allocate();

   };

   #ifndef PSPC_FIELD_STATE_TPP
   // Suppress implicit instantiation
   extern template class BasisFieldState<1>;
   extern template class BasisFieldState<2>;
   extern template class BasisFieldState<3>;
   #endif

} // namespace Pspc
} // namespace Pscf
#endif
