#ifdef  UTIL_MPI
/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "MpiLogger.h"
#include "MpiSendRecv.h"

namespace Util
{

   /*
   * Constructor.
   */
   MpiLogger::MpiLogger(MPI::Intracomm& comm)
    : communicatorPtr_(&comm),
      rank_(-1),
      size_(-1)
   {}

   /*
   */
   void MpiLogger::begin()
   {
      communicatorPtr_->Barrier();
      rank_ = communicatorPtr_->Get_rank();
      size_ = communicatorPtr_->Get_size();
      int data;
      if (rank_ > 0) {
         recv<int>(*communicatorPtr_, data, rank_ - 1, 0);
      } else {
         std::cout << std::endl;
         std::cout.flush();
      }
   }

   /*
   */
   void MpiLogger::end()
   {
      std::cout.flush();
      if (rank_ < size_ - 1) {
         send<int>(*communicatorPtr_, rank_, rank_ + 1, 0);
      }
      communicatorPtr_->Barrier();
   }

}
#endif
