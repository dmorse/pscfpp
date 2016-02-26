#ifndef UTIL_SYMM_DMATRIX_PARAM_H
#define UTIL_SYMM_DMATRIX_PARAM_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/Parameter.h>
#include <util/containers/DMatrix.h>
#include <util/format/Int.h>
#ifdef UTIL_MPI
#include <util/mpi/MpiSendRecv.h>
#endif
#include <util/global.h>

#include <iomanip> 

namespace Util
{

   /**
   * A Parameter associated with a symmetric DMatrix.
   *
   * \ingroup Param_Module
   */
   template <class Type>
   class SymmDMatrixParam : public Parameter
   {
      
   public:
  
      /** 
      * Constructor.
      *
      * \param label  parameter label (a literal C-string)
      * \param matrix  DMatrix<Type> object
      * \param n number of rows or columns
      * \param isRequired  Is this a required parameter?
      */
      SymmDMatrixParam(const char *label, DMatrix<Type>& matrix, int n, bool isRequired = true);
 
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
   
      /// Number of rows or columns
      int n_; 
   
   };

   /*
   * DMatrix constructor.
   */
   template <class Type>
   SymmDMatrixParam<Type>::SymmDMatrixParam(const char* label, DMatrix<Type>& matrix, int n, bool isRequired)
    : Parameter(label, isRequired),
      matrixPtr_(&matrix),
      n_(n)
   {}

   /*
   * Read a DMatrix from isteam.
   */
   template <class Type>
   void SymmDMatrixParam<Type>::readValue(std::istream &in)
   {  
      // Preconditions
      if (!(matrixPtr_->isAllocated())) {
         UTIL_THROW("Cannot read unallocated DMatrix");
      }
      if (n_ != matrixPtr_->capacity1()) {
         UTIL_THROW("Error: Logical size n_ != DMatrix<Type>::capacity1()");
      }
      if (n_ != matrixPtr_->capacity2()) {
         UTIL_THROW("Error: Logical size n_ != DMatrix<Type>::capacity2()");
      }

      // Create matrix of int/bool flags for error checking
      DMatrix<int> flags;
      flags.allocate(n_, n_);
      int i, j;
      for (i = 0; i < n_; ++i) {
         for (j = 0; j < n_; ++j) {
            flags(i,j) = 0;
         }
      }

      double value;
      for (int k = 0; k < n_*(n_ + 1)/2; ++k) {
         in >> i >> j >> value;
         UTIL_CHECK(flags(i,j) == 0);
         (*matrixPtr_)(i, j) = value;
         flags(i, j) = 1;
         if (i != j) {
            UTIL_CHECK(flags(j,i) == 0);
            (*matrixPtr_)(j, i) = value;
            flags(j, i) = 1;
         }
      }

   }

   /*
   * Load a DMatrix from input archive.
   */
   template <class Type>
   void SymmDMatrixParam<Type>::loadValue(Serializable::IArchive& ar)
   {  
      if (!(matrixPtr_->isAllocated())) {
         matrixPtr_->allocate(n_, n_);
      } else {
         if (n_ != matrixPtr_->capacity1()) {
            UTIL_THROW("Error: Logical size n_ != DMatrix<Type>::capacity1()");
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
   void SymmDMatrixParam<Type>::saveValue(Serializable::OArchive& ar)
   {
      if (n_ != matrixPtr_->capacity1()) {
         UTIL_THROW("Error: Logical size n_ != DMatrix<Type>::capacity1()");
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
   void SymmDMatrixParam<Type>::bcastValue()
   {  
      if (!(matrixPtr_->isAllocated())) {
         matrixPtr_->allocate(n_, n_);
      } else {
         if (n_ != matrixPtr_->capacity1()) {
            UTIL_THROW("Error: Logical size n_ > DMatrix<Type>::capacity1()");
         }
         if (n_ != matrixPtr_->capacity2()) {
            UTIL_THROW("Error: Logical size n_ > DMatrix<Type>::capacity2()");
         }
      }
      bcast<Type>(ioCommunicator(), *matrixPtr_, n_, n_, 0); 
   }
   #endif

   /*
   * Write a SymmDMatrixParam.
   */
   template <class Type>
   void SymmDMatrixParam<Type>::writeParam(std::ostream &out)
   {
      if (isActive()) {
         // Preconditions
         if (!(matrixPtr_->isAllocated())) {
            UTIL_THROW("Cannot read unallocated DMatrix");
         }
         if (n_ > matrixPtr_->capacity1()) {
            UTIL_THROW("Error: Logical size n_ > DMatrix<Type>::capacity1()");
         }
         if (n_ > matrixPtr_->capacity2()) {
            UTIL_THROW("Error: Logical size n_ > DMatrix<Type>::capacity2()");
         }
   
         Label space("");
         int i, j;
         for (i = 0; i < n_; ++i) {
            for (j = 0; j <= i; ++j) {
               if (i == 0 && j == 0) {
                  out << indent() << label_;
               } else {
                  out << indent() << space;
               }
               out << Int(i) << "  " << Int(j) << "  "
                   << std::right << std::scientific 
                   << std::setprecision(Parameter::Precision) 
                   << std::setw(Parameter::Width)
                   << (*matrixPtr_)(i, j)
                   << std::endl;
            }
         }
      }
   }

} 
#endif
