#ifndef RPG_FILM_ENVIRONMENT_H
#define RPG_FILM_ENVIRONMENT_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/environment/FilmEnvironment.h> // base class template
#include <pscf/cuda/CudaTp.h>               // base class argument

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class FilmEnvironment<1, CudaTp<1> >;
      extern template class FilmEnvironment<2, CudaTp<2> >;
      extern template class FilmEnvironment<3, CudaTp<3> >;
   } 
} 
#endif
