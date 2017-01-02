#ifndef UTIL_D_ARRAY_PARAM_H
#define UTIL_D_ARRAY_PARAM_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/Parameter.h>    // base class
#include <util/containers/DArray.h>  // member
#include <util/global.h>

#include <iomanip> 

namespace Util
{

   /**
   * A Parameter associated with a DArray container.
   *
   * \ingroup Param_Module
   */
   template <class Type>
   class DArrayParam : public Parameter
   {

   public:

      /*   
      * Constructor.
      */
      DArrayParam(const char *label, DArray<Type>& array, int n, bool isRequired = true);
 
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
   
      /// Pointer to associated DArray.
      DArray<Type>* arrayPtr_;
   
      /// Logical array dimension
      int n_;

   };

   /*
   * DArrayParam<Type> constructor.
   */
   template <class Type>
   DArrayParam<Type>::DArrayParam(const char *label, DArray<Type>& array, int n, bool isRequired)
    : Parameter(label, isRequired),
      arrayPtr_(&array),
      n_(n)
   {}

   /*
   * Read array of values from isteam.
   */
   template <class Type>
   void DArrayParam<Type>::readValue(std::istream &in)
   {  
      if (!(arrayPtr_->isAllocated())) {
         UTIL_THROW("Cannot read unallocated DArray");
      }
      if (arrayPtr_->capacity() != n_) {
         UTIL_THROW("Error: DArray capacity < n");
      }
      for (int i = 0; i < n_; ++i) {
         in >> (*arrayPtr_)[i];
      }
   }

   /*
   * Load a DArray from input archive.
   */
   template <class Type>
   void DArrayParam<Type>::loadValue(Serializable::IArchive& ar)
   {  
      if (!(arrayPtr_->isAllocated())) {
         arrayPtr_->allocate(n_);
      }
      ar >> *arrayPtr_;
      if (arrayPtr_->capacity() < n_) {
         UTIL_THROW("Error: DArray capacity < n");
      }
   }

   /*
   * Save a DArray to an output archive.
   */
   template <class Type>
   void DArrayParam<Type>::saveValue(Serializable::OArchive& ar)
   {  
      if (!(arrayPtr_->isAllocated())) {
         UTIL_THROW("Cannot save unallocated DArray");
      }
      if (arrayPtr_->capacity() != n_) {
         UTIL_THROW("Error: DArray capacity < n");
      }
      ar << *arrayPtr_;
   }

   #ifdef UTIL_MPI
   /*
   * Broadcast a DArray.
   */
   template <class Type>
   void DArrayParam<Type>::bcastValue()
   {  bcast<Type>(ioCommunicator(), *arrayPtr_, n_, 0); }
   #endif

   /*
   * Write a DArray parameter.
   */
   template <class Type>
   void DArrayParam<Type>::writeParam(std::ostream &out) 
   {
      if (isActive()) {

         if (!(arrayPtr_->isAllocated())) {
            UTIL_THROW("Cannot write unallocated DArray");
         }
         if (arrayPtr_->capacity() != n_) {
            UTIL_THROW("Error: DArray capacity != n in writeParam");
         }
   
         Label space("");
         int i;
         for (i = 0; i < n_; ++i) {
            if (i == 0) {
               out << indent() << label_;
            } else {
               out << indent() << space;
            }
            out << std::right << std::scientific 
                << std::setprecision(Parameter::Precision) 
                << std::setw(Parameter::Width)
                << (*arrayPtr_)[i] 
                << std::endl;
         }

      } // if isActive
   }

} 
#endif
