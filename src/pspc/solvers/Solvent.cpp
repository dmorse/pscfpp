/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2019, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Solvent.tpp"

namespace Pscf {
namespace Pspc { 

   #if 0
   template class Solvent<1>;
   template class Solvent<2>;
   template class Solvent<3>;

   template std::istream& operator >> (std::istream& , Solvent<1> &);
   template std::ostream& operator << (std::ostream& , const Solvent<1>&);
   template std::istream& operator >> (std::istream& , Solvent<2>&);
   template std::ostream& operator << (std::ostream& , const Solvent<2>&);
   template std::istream& operator >> (std::istream& , Solvent<3>& );
   template std::ostream& operator << (std::ostream& , const Solvent<3>& );
   #endif

}
}
