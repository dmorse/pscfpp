#ifndef UTIL_PAIR_H
#define UTIL_PAIR_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/containers/FArray.h>
#include <iostream>

namespace Util
{

   /**
   * An array of exactly 2 objects.
   *
   * \ingroup Array_Module
   */
   template <typename Data>
   class Pair : public FArray <Data, 2>
   {

      #ifdef UTIL_MPI
      public:

      /**
      * Commit associated MPI DataType.
      */
      static void commitMpiType();
      #endif

   };

   /** 
   * Input a Pair from an istream.
   *
   * \param in   istream from which to read
   * \param pair Pair to be read
   */
   template <typename Data>
   std::istream& operator>>(std::istream& in, Pair<Data> &pair)
   {
      in >> pair[0] >> pair[1];
      return in;
   }
   
   /**
   * Output a Pair to an ostream, without line breaks.
   *
   * \param out  ostream to which to write
   * \param pair Pair to be written
   */
   template <typename Data>
   std::ostream& operator<<(std::ostream& out, const Pair<Data> &pair) 
   {
      out << "    " << pair[0] << "    " << pair[1];
      return out;
   }

   #ifdef UTIL_MPI

   /*
   * Commit associated MPI Datatype.
   */
   template <typename Data>
   void Pair<Data>::commitMpiType() 
   {
      if (!MpiTraits< Pair<Data> >::hasType) {
         MpiStructBuilder builder;
         Pair<Data>       object;
         builder.setBase(&object);
         builder.addMember(&object[0], MpiTraits<Data>::type);
         builder.addMember(&object[1], MpiTraits<Data>::type);
         builder.commit(MpiTraits< Pair<Data> >::type);
         MpiTraits< Pair<Data> >::hasType = true;
      }
   }

   #endif 
} 
#endif
