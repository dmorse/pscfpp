#ifdef  UTIL_MPI
#ifndef UTIL_STRUCT_BUILDER_H
#define UTIL_STRUCT_BUILDER_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/global.h>

namespace Util
{

   /**
   * A MpiStructBuilder objects is used to create an MPI Struct datatype.
   *
   * This class provides methods to simplify construction of an MPI data 
   * type that can stores instances of a C struct or C++ class. 
   *
   * As an example, consider the creation of an MPI datatype MyClassMpi 
   * for class MyClass, with a class definition:
   * \code
   *
   * class MyClass 
   * {
   *    double x[3];
   *    int    i, j;
   * }
   *
   * \endcode
   * The code required to build and commit the MPI datatype MyClassMpi is:
   * \code
   *
   *    MyClass          object;
   *    MPI::Datatype  MyClassMpi;
   *    MpiStructBuilder builder;
   *
   *    builder.setBase(&object)
   *    builder.addMember(&object.x, MPI::DOUBLE, 3);
   *    builder.addMember(&object.i, MPI::INT, 1);
   *    builder.addMember(&object.j, MPI::INT, 1);
   *    builder.commit(&MyClassMpi);
   *
   * \endcode
   * The setBase and addMember classes require addresses of an instance 
   * of the class and of its members, respectively. These addresses must
   * all refer to same instance. The commit method calculates the offset
   * of each member by subtracting the address of the object from the
   * address of each of its members.
   */
   class MpiStructBuilder
   {
   
   public:
   
      /// Default constructor.
      MpiStructBuilder();
   
      /**
      * Set address of an class instance.
      */
      void setBase(void* objectAddress); 
   
      /**
      * Add a new member variable to the type map.
      *
      * This method must be called once for each member. The address parameter 
      * must be a pointer to a member variable of the object whose base address 
      * is passed to setBase(). 
      *
      * The count parameter is required only for array members: the default
      * value of count=1 may be used for scalar members. 
      *
      * \param memberAddress  displacement of variable, in bytes.
      * \param type           data type (MPI::INT, MPI::DOUBLE, etc.)
      * \param count          number of contiguous variables (array count)
      */
      void addMember(void* memberAddress, MPI::Datatype type, int count = 1);
   
      /**
      * Build and commit a user-defined MPI Struct datatype.
      *
      * The setBase() method must be called once and the addMember() method 
      * must be called once per member before calling this method.
      *
      * \param newType new MPI datatype (on output).
      */
      void commit(MPI::Datatype& newType);
   
   private:
   
      static const int MaxNBlock = 20;      // Maximum allowed number of members
   
      MPI::Aint     base_;                  // address of example object
      MPI::Datatype types_[MaxNBlock];      // datatypes of members
      MPI::Aint     addresses_[MaxNBlock];  // addresses of members
      int           counts_[MaxNBlock];     // counts of members
      int           nMember_;               // number of members added
   
   };
}
#endif
#endif
