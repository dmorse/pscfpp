#ifdef UTIL_MPI

#include "MpiStructBuilder.h"
#include <util/global.h>

namespace Util
{

   /// Default constructor.
   MpiStructBuilder::MpiStructBuilder() 
    : base_(0),
      nMember_(0)
   {}
   
   /* 
   * Record address of a pointer to an instance of the class.
   */
   void MpiStructBuilder::setBase(void* objectAddress) 
   { 
      base_ = MPI::Get_address(objectAddress);
   }

   /* 
   * Add a member variable to the struct definition
   */
   void MpiStructBuilder::addMember(void* memberAddress, MPI::Datatype type, int count)
   {
      addresses_[nMember_]  = MPI::Get_address(memberAddress);
      types_[nMember_]      = type;
      counts_[nMember_]     = count;
      ++nMember_;
   }

   /*
   * Build and commit a user-defined MPI Struct datatype.
   *
   * \param mpiType new MPI datatype (on output).
   */
   void MpiStructBuilder::commit(MPI::Datatype& mpiType) 
   {
      for (int i = 0; i < nMember_; ++i) {
         addresses_[i] = addresses_[i] - base_;
      }
      mpiType = 
            MPI::Datatype::Create_struct(nMember_, counts_, addresses_, types_);
      mpiType.Commit();
   }

}
#endif
