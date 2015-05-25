#ifndef UTIL_CARRAY_PARAM_H
#define UTIL_CARRAY_PARAM_H

#include <util/param/Parameter.h>
#include <util/global.h>

#ifdef UTIL_MPI
#include <util/mpi/MpiSendRecv.h>
#endif

#include <iomanip> 

namespace Util
{

   /**
   * A Parameter associated with a 1D C array. 
   *
   * \ingroup Param_Module
   */
   template <class Type>
   class CArrayParam : public Parameter
   {
   
   public:

      /**   
      * Constructor.
      */
      CArrayParam(const char *label, Type *value, int n, bool isRequired = true);

      /** 
      * Write parameter to stream.
      *
      * \param out output stream
      */
      void writeParam(std::ostream &out);

   protected:
      
      /**
      * Read parameter value from an input stream.
      * 
      * \param in input stream from which to read
      */
      virtual void readValue(std::istream& in);

      /**
      * Load bare parameter value from an archive.
      *
      * \param ar input archive from which to load
      */
      virtual void loadValue(Serializable::IArchive& ar);

      /**
      * Save parameter value to an archive.
      *
      * \param ar output archive to which to save
      */
      virtual void saveValue(Serializable::OArchive& ar);

      #ifdef UTIL_MPI
      /**
      * Broadcast parameter value within the ioCommunicator.
      */
      virtual void bcastValue();
      #endif

   private:
   
      /// Pointer to value.
      Type* value_;
   
      /// Array dimension
      int n_;    
   
   };

   /*
   * CArrayParam<Type> constructor.
   */
   template <class Type>
   CArrayParam<Type>::CArrayParam(const char *label, Type* value, int n, bool isRequired)
    : Parameter(label, isRequired),
      value_(value),
      n_(n)
   {}

   /*
   * Read C-array of n values from file.
   */
   template <class Type>
   void CArrayParam<Type>::readValue(std::istream &in)
   {  
      for (int i = 0; i < n_; ++i) {
         in >> value_[i];
      }
   }

   /*
   * Load C-array of n values from an input archive
   */
   template <class Type>
   void CArrayParam<Type>::loadValue(Serializable::IArchive& ar)
   {  ar.unpack(value_, n_); }

   /*
   * Save C-array of n values to an output archive
   */
   template <class Type>
   void CArrayParam<Type>::saveValue(Serializable::OArchive& ar)
   {  ar.pack(value_, n_); }

   #ifdef UTIL_MPI
   /*
   * Broadcast an array of n values
   */
   template <class Type>
   void CArrayParam<Type>::bcastValue()
   {  bcast<Type>(ioCommunicator(), value_, n_, 0); }
   #endif

   /*
   * Write a C array
   */
   template <class Type>
   void CArrayParam<Type>::writeParam(std::ostream &out) 
   {
      if (isActive()) {
         Label space("");
         for (int i = 0; i < n_; ++i) {
            if (i == 0) {
               out << indent() << label_;
            } else {
               out << indent() << space;
            }
            out << std::right << std::scientific 
                << std::setprecision(Parameter::Precision) 
                << std::setw(Parameter::Width)
                << value_[i] 
                << std::endl;
         }
      }
   }

} 
#endif
