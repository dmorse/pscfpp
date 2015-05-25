#include "PointGroup.h"

namespace Util {

   /**
   * Output stream inserter operator writing a PointGroup to stream
   */ 
   std::ostream& operator << (std::ostream& out, const PointGroup& g)
   {
      int i, size;
      size = g.size();
      out << "size = " << size << std::endl;
      for (i = 0; i < size; ++i) {
         out << std::endl;
         out << g[i];
      }
      return out;
   }

}
