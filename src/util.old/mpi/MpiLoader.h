#ifndef UTIL_MPI_LOADER_H
#define UTIL_MPI_LOADER_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/containers/DArray.h>   // template argument
#include <util/containers/FArray.h>   // template argument
#include <util/containers/DMatrix.h>  // template argument
#include <util/mpi/MpiFileIo.h>       // used in template implementation
#include <util/mpi/MpiTraits.h>       // used in template implementation
#include <util/mpi/MpiSendRecv.h>     // used in template implementation
#include <util/global.h>


namespace Util
{

   /**
   * Provides methods for MPI-aware loading of data from input archive.
   *
   * Each MpiLoader is associated with an IArchive input archive, and 
   * with a MpiFileIo, which are passed as arguments to the constructor. 
   * The MpiFileIo argument is often a ParamComposite, which is derived 
   * from MpiFileIo.
   *
   * The "load" function templates all load data from the archive and (if
   * appropriate) broadcast data among processors. If MPI is not enabled
   * (i.e., if UTIL_MPI is not defined), then the data is simply loaded
   * from the archive. If MPI is enabled and a parameter communicator is 
   * set, data is loaded from the archive by the ioProcessor and then 
   * broadcast to all other processors in the IO communicator. If MPI is 
   * enabled but no parameter communicator is set, every processor loads 
   * data independently.
   *
   * \ingroup Mpi_Module
   */
   template <class IArchive>
   class MpiLoader 
   {

   public:

      /**
      * Constructor
      *  
      * \param mpiFileIo associated MpiFileIo object
      * \param archive   input archive from which data will be loaded
      */
      MpiLoader(MpiFileIo& mpiFileIo, IArchive& archive);

      /**  
      * Load and broadcast a single Data value.
      *
      * \param value  reference to a Data
      */
      template <typename Data>
      void load(Data &value);
   
      /**  
      * Load and broadcast a C array.
      *
      * \param value  pointer to array
      * \param n      number of elements
      */
      template <typename Data>
      void load(Data *value, int n);
   
      /**  
      * Load and broadcast a DArray < Data > container.
      *
      * \param array  DArray object
      * \param n      number of elements
      */
      template <typename Data>
      void load(DArray<Data>& array, int n);
   
      /**  
      * Load and broadcast an FArray <Data , N > object.
      *
      * \param array  FArray object to be loaded
      */
      template <typename Data, int N>
      void load(FArray<Data, N >& array);
   
      /**  
      * Load and broadcast a 2D CArray of Data objects.
      *
      * Loads m rows of n elements into array declared as Data array[][np].
      *
      * \param value pointer to first element or row in array
      * \param m  logical number of rows (1st dimension)
      * \param n  logical number of columns (2nd dimension)
      * \param np  physcial number of columns (elements allocated per row)
      */
      template <typename Data> void 
      load(Data* value, int m, int n, int np);
  
      /**  
      * Load and broadcast a DMatrix<Data> object.
      *
      * \param matrix DMatrix object
      * \param m  number of rows (1st dimension)
      * \param n  number of columns (2nd dimension)
      */
      template <typename Data> 
      void load(DMatrix<Data>& matrix, int m, int n);

   private:

      // Pointer to associated MpiFileIo (passed to constructor).
      MpiFileIo*  mpiFileIoPtr_;

      // Pointer to associated input archive (passed to constructor).
      IArchive*  archivePtr_;
 
   };
 
   /*
   * Constructor
   */
   template <typename IArchive> 
   MpiLoader<IArchive>::MpiLoader(MpiFileIo& mpiFileIo, IArchive& archive)
    : mpiFileIoPtr_(&mpiFileIo),
      archivePtr_(&archive)
   {}

   /*  
   * Load and broadcast a single Data value.
   */
   template <typename IArchive> 
   template <typename Data> 
   void MpiLoader<IArchive>::load(Data &value)
   {
      if (mpiFileIoPtr_->isIoProcessor()) {
         *archivePtr_ >> value;
      }
      #ifdef UTIL_MPI
      if (mpiFileIoPtr_->hasIoCommunicator()) {
         bcast<Data>(mpiFileIoPtr_->ioCommunicator(), value, 0); 
      }
      #endif
   }
   
   /*  
   * Add a C array parameter, and load its elements.
   */
   template <typename IArchive> 
   template <typename Data>
   void MpiLoader<IArchive>::load(Data* value, int n)
   {
      if (mpiFileIoPtr_->isIoProcessor()) {
         for (int i = 0; i < n; ++i) {
            *archivePtr_ >> value[i];
         }
      }
      #ifdef UTIL_MPI
      if (mpiFileIoPtr_->hasIoCommunicator()) {
         bcast<Data>(mpiFileIoPtr_->ioCommunicator(), value, n, 0); 
      }
      #endif
   }

   /*
   * Load a DArray < Data > container.
   */
   template <typename IArchive> 
   template <typename Data>
   void 
   MpiLoader<IArchive>::load(DArray<Data>& array, int n)
   {
      if (mpiFileIoPtr_->isIoProcessor()) {
         *archivePtr_ >> array;
         if (array.capacity() < n) {
            UTIL_THROW("Error: DArray capacity < n");
         }
      }
      #ifdef UTIL_MPI
      if (mpiFileIoPtr_->hasIoCommunicator()) {
         bcast<Data>(mpiFileIoPtr_->ioCommunicator(), array, n, 0); 
      }
      #endif
   }

   /*
   * Load an FArray < Data, N > fixed-size 1D array container.
   */
   template <typename IArchive> 
   template <typename Data, int N>
   void MpiLoader<IArchive>::load(FArray<Data, N >& array)
   {
      if (mpiFileIoPtr_->isIoProcessor()) {
         for (int i = 0; i < N; ++i) {
            *archivePtr_ >> array[i];
         }
      }
      #ifdef UTIL_MPI
      if (mpiFileIoPtr_->hasIoCommunicator()) {
         bcast<Data>(mpiFileIoPtr_->ioCommunicator(), &(array[0]), N, 0); 
      }
      #endif
   }

   /*
   * Load a CArray2DParam < Data > C two-dimensional array parameter.
   */
   template <typename IArchive> 
   template <typename Data> 
   void 
   MpiLoader<IArchive>::load(Data *array, int m, int n, int np)
   {
      if (mpiFileIoPtr_->isIoProcessor()) {
         int i, j; 
         for (i = 0; i < m; ++i) {
            for (j = 0; j < n; ++j) {
               *archivePtr_ >> array[i*np + j];
            }
         }
      }
      #ifdef UTIL_MPI
      if (mpiFileIoPtr_->hasIoCommunicator()) {
         // Broadcast block of m rows of np elements each.
         bcast<Data>(mpiFileIoPtr_->ioCommunicator(), &(array[0]), m*np, 0); 
      }
      #endif
   }
  
   /*
   * Add and load a DMatrix < Data > C two-dimensional matrix parameter.
   */
   template <typename IArchive> 
   template <typename Data> 
   void MpiLoader<IArchive>::load(DMatrix<Data>& matrix, int m, int n)
   {
      if (mpiFileIoPtr_->isIoProcessor()) {
         *archivePtr_ >> matrix;
      }
      #ifdef UTIL_MPI
      if (mpiFileIoPtr_->hasIoCommunicator()) {
         bcast<Data>(mpiFileIoPtr_->ioCommunicator(), matrix, m, n, 0); 
      }
      #endif
   }
  
} 
#endif
