#ifndef PSSP_GPU_RING_PBUFFER_H
#define PSSP_GPU_RING_PBUFFER_H


#include <util/containers/RingBuffer.h>

namespace Pscf {
namespace Pssp_gpu 
{
	template< class Data >
	class RingPBuffer : public RingBuffer <Data>
	{
	public:
		/**
      * Homebrew version of the parent class
      * Add a new value to the buffer. This class expect data to be a pointer
      * Will delete the pointer to an old history before appending a new member
      * This method provides better runtime for large Data class where memory
      * transfer might be an issue
      * \param value new value to be added.
      */
      void append(Data* value)
      {
	      assert(last_ < capacity_);
	      ++last_;
	      if (last_ == capacity_) {
	         last_ = 0; // wrap around
	      };
	      
	      //if occupied. i.e. pointer is not null;
	      Data* ptr = data_ + last_;
	      if(!(ptr))
	      	delete (ptr);

	      ptr = data;
	   
	      // If buffer is not yet full, update size_
	      if (size_ < capacity_) {
	         ++size_;
      };
      }
	private:
	};
}
}


#endif