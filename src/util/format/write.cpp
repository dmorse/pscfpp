/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/


#include "write.h"
#include "Dbl.h"
#include "Int.h"
#include "Lng.h"
#include "Bool.h"
//#include "Str.h"

namespace Util
{

   // Explicit specialization for double
   template <> void write(std::ostream& out, double data)
   { out << Dbl(data); }
   
   // Explicit specialization for double
   template <> void write(std::ostream& out, std::complex<double> data)
   { out << Dbl(data.real()) << Dbl(data.imag()); }

   // Explicit specialization for int
   template <> void write(std::ostream& out, int data)
   { out << Int(data); }
   
   // Explicit specialization for long
   template <> void write(std::ostream& out, long data)
   { out << Lng(data); }
   
   // Explicit specialization for bool
   template <> void write(std::ostream& out, bool data)
   { out << Bool(data); }

}
