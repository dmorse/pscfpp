#ifndef UTIL_DMATRIX_PARAM_H
#define UTIL_DMATRIX_PARAM_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/Parameter.h>
#include <util/containers/DMatrix.h>
#ifdef UTIL_MPI
#include <util/mpi/MpiSendRecv.h>
#endif
#include <util/global.h>

#include <iomanip> 

namespace Util
{

   /**
   * A Parameter associated with a 2D built-in C array.
   *
   * \ingroup Param_Module
   */
   template <class Type>
   class DMatrixParam : public Parameter
   {
      
   public:
  
      /** 
      * Constructor.
      *
      * \param label  parameter label (a literal C-string)
      * \param matrix  DMatrix<Type> object
      * \param m  number of rows
      * \param n  number of columns
      * \param isRequired  Is this a required parameter?
      */
      DMatrixParam(const char *label, DMatrix<Type>& matrix, int m, int n, bool isRequired = true);
 
      /**
      * Write DMatrix to file.
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
   
      /// Pointer to associated DMatrix.
      DMatrix<Type>* matrixPtr_;
   
      /// Number of rows in array[m][n]
      int m_; 

      /// Number of columns in array[m][n]
      int n_; 
   
   };

   /*
   * DMatrix constructor.
   */
   template <class Type>
   DMatrixParam<Type>::DMatrixParam(const char* label, DMatrix<Type>& matrix, int m, int n, bool isRequired)
    : Parameter(label, isRequired),
      matrixPtr_(&matrix),
      m_(m),
      n_(n)
   {}

   /*
   * Read a DMatrix from isteam.
   */
   template <class Type>
   void DMatrixParam<Type>::readValue(std::istream &in)
   {  
      // Preconditions
      if (!(matrixPtr_->isAllocated())) {
         UTIL_THROW("Cannot read unallocated DMatrix");
      }
      if (m_ != matrixPtr_->capacity1()) {
         UTIL_THROW("Error: Logical size m_ != DMatrix<Type>::capacity1()");
      }
      if (n_ != matrixPtr_->capacity2()) {
         UTIL_THROW("Error: Logical size n_ != DMatrix<Type>::capacity2()");
      }

      int i, j;
      for (i = 0; i < m_; ++i) {
         for (j = 0; j < n_; ++j) {
            in >> (*matrixPtr_)(i, j);
         }
      }
   }

   /*
   * Load a DMatrix from input archive.
   */
   template <class Type>
   void DMatrixParam<Type>::loadValue(Serializable::IArchive& ar)
   {  
      if (!(matrixPtr_->isAllocated())) {
         matrixPtr_->allocate(m_, n_);
      } else {
         if (m_ != matrixPtr_->capacity1()) {
            UTIL_THROW("Error: Logical size m_ != DMatrix<Type>::capacity1()");
         }
         if (n_ != matrixPtr_->capacity2()) {
            UTIL_THROW("Error: Logical size n_ != DMatrix<Type>::capacity2()");
         }
      }
      ar >> *matrixPtr_;
   }

   /*
   * Save a DMatrix to an output archive.
   */
   template <class Type>
   void DMatrixParam<Type>::saveValue(Serializable::OArchive& ar)
   {
      if (m_ != matrixPtr_->capacity1()) {
         UTIL_THROW("Error: Logical size m_ != DMatrix<Type>::capacity1()");
      }
      if (n_ != matrixPtr_->capacity2()) {
         UTIL_THROW("Error: Logical size n_ != DMatrix<Type>::capacity2()");
      }
      ar << *matrixPtr_; 
   }

   #ifdef UTIL_MPI
   /*
   * Broadcast a DMatrix.
   */
   template <class Type>
   void DMatrixParam<Type>::bcastValue()
   {  
      if (!(matrixPtr_->isAllocated())) {
         matrixPtr_->allocate(m_, n_);
      } else {
         if (m_ != matrixPtr_->capacity1()) {
            UTIL_THROW("Error: Logical size m_ > DMatrix<Type>::capacity1()");
         }
         if (n_ != matrixPtr_->capacity2()) {
            UTIL_THROW("Error: Logical size n_ > DMatrix<Type>::capacity2()");
         }
      }
      bcast<Type>(ioCommunicator(), *matrixPtr_, m_, n_, 0); 
   }
   #endif

   /*
   * Write a DMatrixParam.
   */
   template <class Type>
   void DMatrixParam<Type>::writeParam(std::ostream &out)
   {
      if (isActive()) {
         // Preconditions
         if (!(matrixPtr_->isAllocated())) {
            UTIL_THROW("Cannot read unallocated DMatrix");
         }
         if (m_ > matrixPtr_->capacity1()) {
            UTIL_THROW("Error: Logical size m_ > DMatrix<Type>::capacity1()");
         }
         if (n_ > matrixPtr_->capacity2()) {
            UTIL_THROW("Error: Logical size n_ > DMatrix<Type>::capacity2()");
         }
   
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
                   << (*matrixPtr_)(i, j);
            }
            out << std::endl;
         }
      }
   }

} 
#endif
