#ifndef UTIL_CARRAY_2D_PARAM_H
#define UTIL_CARRAY_2D_PARAM_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/Parameter.h>
#include <util/global.h>

#ifdef UTIL_MPI
#include <util/mpi/MpiSendRecv.h>
#endif

#include <iomanip> 

namespace Util
{

   /**
   * A Parameter associated with a 2D built-in C array.
   *
   * \ingroup Param_Module
   */
   template <class Type>
   class CArray2DParam : public Parameter
   {
      
   public:
  
      /** 
      * Constructor.
      *
      * Example: A 2 X 2 matrix stored in an oversized 3 x 3 C array.
      * \code
      *    double                matrix[3][3];
      *    CArray2DParam<double> param("matrix", matrix[0], 2, 2, 3);
      * \endcode
      *
      * \param label  parameter label (usually a literal C-string)
      * \param ptr  pointer to first element of first row of 2D array
      * \param m  logical number of rows
      * \param n  logical number of columns
      * \param np  physical number of columns (allocated elements per row).
      * \param isRequired Is this a required parameter?
      */
      CArray2DParam(const char *label, Type* ptr, 
                    int m, int n, int np, bool isRequired = true);
 
      /**
      * Write 2D C array to file.
      */ 
      void writeParam(std::ostream &out);

   protected:
      
      /**
      * Read 2D array parameter from an input stream.
      * 
      * \param in input stream from which to read
      */
      virtual void readValue(std::istream& in);

      /**
      * Load 2D array from an archive.
      *
      * \param ar input archive from which to load
      */
      virtual void loadValue(Serializable::IArchive& ar);

      /**
      * Save 2D array to an archive.
      *
      * \param ar output archive to which to save
      */
      virtual void saveValue(Serializable::OArchive& ar);

      #ifdef UTIL_MPI
      /**
      * Broadcast 2D array within the ioCommunicator.
      */
      virtual void bcastValue();
      #endif

   private:
   
      /// Pointer to first element of first row in associated 2D C array
      Type* ptr_;
   
      /// Number of rows in array[][np]
      int m_; 

      /// Logical number of columns in array[][np]
      int n_; 

      /// Physical number of columns in array[][np]
      int np_; 
   
   };


   /*
   * CArray2D constructor.
   */
   template <class Type>
   CArray2DParam<Type>::CArray2DParam(const char* label, Type* ptr, int m, int n, int np, bool isRequired)
    : Parameter(label, isRequired),
      ptr_(ptr),
      m_(m),
      n_(n),
      np_(np)
   {}

   /*
   * Read a DArray from isteam.
   */
   template <class Type>
   void CArray2DParam<Type>::readValue(std::istream &in)
   {  
      int i, j;
      for (i = 0; i < m_; ++i) {
         for (j = 0; j < n_; ++j) {
            in >> ptr_[i*np_ + j];
         }
      }
   }

   /*
   * Load a DArray from input archive.
   */
   template <class Type>
   void CArray2DParam<Type>::loadValue(Serializable::IArchive& ar)
   {  ar.unpack(ptr_, m_, n_, np_); }

   /*
   *  Save a DArray to an output archive.
   */
   template <class Type>
   void CArray2DParam<Type>::saveValue(Serializable::OArchive& ar)
   {  ar.pack(ptr_, m_, n_, np_); }

   #ifdef UTIL_MPI
   /*
   * Broadcast a DArray.
   */
   template <class Type>
   void CArray2DParam<Type>::bcastValue()
   {  bcast<Type>(ioCommunicator(), ptr_, m_*np_, 0); }
   #endif

   /*
   * Write a CArray2DParam.
   */
   template <class Type>
   void CArray2DParam<Type>::writeParam(std::ostream &out)
   {
      if (isActive()) {
         Label space("");
         int i, j;
   
         for (i = 0; i < m_; ++i) {
            if (i == 0) {
               out << indent() << label_;
            } else {
               out << indent() << space;
            }
            for (j = 0; j < n_; ++j) {
               out << std::right << std::scientific 
                   << std::setprecision(Parameter::Precision) 
                   << std::setw(Parameter::Width)
                   << ptr_[i*np_ + j];
            }
            out << std::endl;
         }
      }
   }
} 
#endif
