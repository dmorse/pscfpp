#ifdef  UTIL_MPI

#include <string.h>
/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "MpiSendRecv.h"

namespace Util
{

   // bool explicit specializations
 
   template <>
   void send<bool>(MPI::Comm& comm, bool& data, int dest, int tag)
   { 
      int value = data ? 1 : 0;
      send<int>(comm, value, dest, tag);
   }

   template <>
   void recv<bool>(MPI::Comm& comm, bool& data, int source, int tag)
   { 
      int value;
      recv<int>(comm, value, source, tag);
      data = value ? true : false;
   }

   template <>
   void bcast<bool>(MPI::Intracomm& comm, bool& data, int root)
   { 
      int value;
      int rank = comm.Get_rank();
      if (rank == root) 
         value = data ? 1 : 0;
      bcast<int>(comm, value, root);
      if (rank != root) 
         data = value ? true : false;
   }

   // std::string explicit specializations
 
   template <>
   void send<std::string>(MPI::Comm& comm, std::string& data, int dest, int tag)
   { 

      // Send size of char C array
      int count = data.size() + 1;
      send<int>(comm, count, dest, tag);

      // Send string contents as char array
      char* cstr = new char[count];
      strcpy(cstr, data.c_str());
      send<char>(comm, cstr, count, dest, tag); 
      delete [] cstr;

   }

   template <>
   void recv<std::string>(MPI::Comm& comm, std::string& data, int source, int tag)
   { 

      // Receive size of char C array
      int  count;
      recv<int>(comm, count, source, tag);

      // Receive contents as char C array
      char* cstr = new char[count];
      recv<char>(comm, cstr, count, source, tag); 
      data = cstr;
      delete [] cstr;

   }

   template <>
   void bcast<std::string>(MPI::Intracomm& comm, std::string& data, int root)
   { 
      int rank = comm.Get_rank();
      int count;

      // Broadcast string count
      if (rank == root) 
         count = data.size() + 1;
      bcast<int>(comm, count, root);

      // Broadcast std::string contents as C string
      char* cstr = new char[count];
      if (rank == root) 
         strcpy(cstr, data.c_str());
      bcast<char>(comm, cstr, count, root); 
      if (rank != root) 
         data = cstr;
      delete [] cstr;

   }

}
#endif
