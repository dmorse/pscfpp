/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "MemoryIArchive.h"
#include "MemoryOArchive.h"
#include <util/misc/Memory.h>

#include <sstream>

namespace Util
{

   /*
   * Constructor.
   */
   MemoryIArchive::MemoryIArchive()
    : buffer_(0),
      begin_(0),
      cursor_(0),
      end_(0),
      endAllocated_(0),
      oArchivePtr_(0),
      capacity_(0),
      version_(0),
      ownsData_(false)
   {}

   /*
   * Destructor.
   */
   MemoryIArchive::~MemoryIArchive()
   {  
      if (buffer_ && ownsData_) {
         //delete buffer_; 
         Memory::deallocate<Byte>(buffer_, capacity_ + sizeof(size_t));
      }
   }

   /*
   * Allocate a block of memory.
   */
   void MemoryIArchive::allocate(size_t capacity)
   {
      if (begin_) {
         UTIL_THROW("Re-allocation is prohibited");
      }
      // buffer_ = new Byte[capacity + sizeof(size_t)];
      Memory::allocate<Byte>(buffer_, capacity + sizeof(size_t));
      begin_  = buffer_ + sizeof(size_t);
      cursor_ = begin_;
      end_ = begin_;
      endAllocated_ = begin_ + capacity; 
      capacity_ = capacity;
      ownsData_ = true; 
   }

   /**
   * Assignment from MemoryOArchive.
   */
   MemoryIArchive& MemoryIArchive::operator = (MemoryOArchive& other)
   {
      if (isAllocated()) {
         UTIL_THROW("This MemoryIArchive is already allocated");
      }
      buffer_ = other.buffer_;
      begin_ = other.begin_;
      cursor_ = other.begin_;
      end_ = other.cursor_;
      endAllocated_ = other.endAllocated_;
      capacity_ = other.capacity_;
      ownsData_ = false;
      oArchivePtr_ = &other;
      other.isLocked_ = true;
      return *this;
   }

   /*
   * Reset cursor to beginning, prepare to reread.
   */
   void MemoryIArchive::reset()
   {
      if (!isAllocated()) {
         UTIL_THROW("Archive is not allocated");
      }  
      cursor_ = begin(); 
   }

   /*
   * Clear memory block (reset to empty state).
   */
   void MemoryIArchive::clear()
   {
      if (!isAllocated()) {
         UTIL_THROW("Archive is not allocated");
      }
      if (!ownsData_) {  
         UTIL_THROW("Archive does not own data");
      }
      cursor_ = begin(); 
      end_ = begin();
   }

   /*
   * Release associated data acquired by assignment.
   */
   void MemoryIArchive::release()
   {
      if (!begin_) {
         UTIL_THROW("This archive has no associated data");
      }
      if (ownsData_) {
         UTIL_THROW("This archive owns its data");
      }
      buffer_ = 0;
      begin_ = 0;
      cursor_ = 0;
      end_ = 0;
      endAllocated_ = 0;
      capacity_ = 0;
      ownsData_ = false;
      oArchivePtr_->isLocked_ = false;
      oArchivePtr_ = 0;
   }

   #ifdef UTIL_MPI
   /*
   * Receive a block.
   */
   void MemoryIArchive::recv(MPI::Intracomm& comm, int source)
   {
      int  myRank     = comm.Get_rank();
      int  comm_size  = comm.Get_size();

      // Preconditions
      if (source > comm_size - 1 || source < 0) {
         UTIL_THROW("Source rank out of bounds");
      }
      if (source == myRank) {
         UTIL_THROW("Source and desination identical");
      }

      size_t recvCapacity = capacity_ + sizeof(size_t);
      comm.Recv(buffer_, recvCapacity, MPI::UNSIGNED_CHAR, source, 5);

      begin_ = buffer_ + sizeof(size_t);
      cursor_ = begin_;

      size_t* sizePtr = (size_t*) buffer_;
      size_t  size = *sizePtr;
      end_  = buffer_ + size;
   }
   #endif

   /*
   * Load a std::string from MemoryIArchive.
   */
   template <>
   void serialize(MemoryIArchive& ar, std::string& data, 
                  const unsigned int version)
   {
      static std::vector<char> charvec;
      size_t size;
      ar.unpack(size);
      if (size > charvec.capacity()) {
         charvec.reserve(size + 8);
      }
      ar.unpack(&charvec[0], size);
      data = &charvec[0];
   }

}
