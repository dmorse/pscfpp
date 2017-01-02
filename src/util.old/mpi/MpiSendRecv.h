#ifdef UTIL_MPI
#ifndef UTIL_MPI_SEND_RECV_H
#define UTIL_MPI_SEND_RECV_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

/**
* \file MpiSendRecv.h
*
* This file contains templates for global functions send<T>, recv<T> 
* and bcast<T>. These are wrappers for the MPI send, recv (receive), and 
* bcast (broadcast) functions. Overloaded forms of these functions are
* provided to transmit single values and arrays. The main difference 
* between the wrappers and the underlying MPI functions is that the 
* wrapper functions do not require an MPI type handle as a parameter. 
* Instead, the MPI type associated with C++ type T (the template parameter) 
* is inferred by the function implementations, by methods that are 
* described below.  The most important advantage of this is that it 
* allows the wrapper functions to be used within other templates that 
* take type T as a template parameter. The corresponding MPI methods 
* cannot be used in generic templates because they require MPI type 
* handle parameters that have different values for different date types.
*
* The implementation of the templates send<T>, recv<T>, bcast<T> for
* single values of type T rely on the existence of an associated explicit 
* specialization of the class template MpiTraits<typename T>. If it 
* exists, the class MpiTraits<T> maps C++ type T onto an associated 
* MPI type.  Each specialization MpiTraits<T> has a static member 
* MpiTraits<T>::type that contains an opaque handle for the MPI type 
* associated with C++ type T. Explicit specializations for the most 
* common built-in C++ types are defined in MpiTraits.h. 
*
* The send<T>, recv<T>, and bcast<T> templates can also be used to 
* transmit instances of a user defined class T if an appropriate MPI 
* type exists.  To make this work, the user must define and commit an 
* associated user-defined MPI data type, and also define an explicit 
* specialization MpiTraits<T> to associate this MPI type with C++ type 
* T. Specialized MPI data types and MpiTraits classes for Util::Vector
* and Util::Vector are defined in the header and implementation files 
* for these classes. User defined MPI types must be committed before 
* they can be used. 
*
* Explicit specializations of send<T>, recv<T> and bcast<T> may also 
* be provided for some types for which the algorithm based on MpiTraits
* is awkward or unworkable. Explicit specializations are declared in
* this file for bool and std::string. The implementations of send<T>
* recv<T>, and bcast<T> for T=bool transmit boolean values as integers. 
* The implementations for T = std::string transmit strings as character 
* arrays. No MpiTraits classes are needed or provided for bool or 
* std::string, because the compiler will always use these explicit 
* specializations, which do not rely on MpiTraits classes, rather than
* the main function templates. It may also be more convenient for some
* user-defined classes to provide explicit specializations of these
* three functions, rather than defining an associated MPI type and 
* MpiTraits specialization.  The templates defined here can be used to 
* transmit instances of type T either if: (i) Explicit specializations 
* are defined for these three functions, or (ii) an associated MpiTraits 
* class and MPI data type are defined.
*
* Overloaded forms of send<T>, recv<T>, and bcast<T> are provided to
* transmit 1D and 2D C arrays of data and DArray<T> and DMatrix<T> 
* containers. These functions send the data in one transmission, as 
* a contiguous buffer, if an MPI type is available, but send each 
* element in a separate transmission if no MpiType exists but an 
* explicit specialization exists for the required scalar form of 
* send, recv, or bcast.
*/

#include <util/global.h>

#include <util/mpi/MpiTraits.h>
#include <util/containers/DArray.h>
#include <util/containers/DMatrix.h>

namespace Util
{

   // Scalar parameters

   /**
   * Send a single T value.
   *
   * Throws an Exception if no associated MPI data type is available, 
   * i.e., if MpiTraits<T>::hasType is false.
   * 
   * \param comm MPI communicator
   * \param data value
   * \param dest MPI rank of receiving processor in comm
   * \param tag  user-defined integer identifier for message
   */
   template <typename T>
   void send(MPI::Comm& comm, T& data, int dest, int tag)
   {
      if (!MpiTraits<T>::hasType)
         UTIL_THROW("No committed MPI type in send<T>");
      comm.Send(&data, 1, MpiTraits<T>::type, dest, tag); 
   }
  
   /**
   * Receive a single T value.
   *
   * Throws an Exception if no associated MPI data type is available, 
   * i.e., if MpiTraits<T>::hasType is false.
   * 
   * \param comm   MPI communicator
   * \param data   value
   * \param source MPI rank of sending processor in comm
   * \param tag    user-defined integer identifier for message
   */
   template <typename T>
   void recv(MPI::Comm& comm, T& data, int source, int tag)
   {  
      if (!MpiTraits<T>::hasType)
         UTIL_THROW("No committed MPI type in recv<T>");
      comm.Recv(&data, 1, MpiTraits<T>::type, source, tag); 
   }
  
   /**
   * Broadcast a single T value.
   *
   * Throws an Exception if no associated MPI data type is available, 
   * i.e., if MpiTraits<T>::hasType is false.
   * 
   * \param comm   MPI communicator
   * \param data   value
   * \param root   MPI rank of root (sending) processor in comm
   */
   template <typename T>
   void bcast(MPI::Intracomm& comm, T& data, int root)
   {  
      if (!MpiTraits<T>::hasType)
         UTIL_THROW("No committed MPI type in bcast<T>");
      comm.Bcast(&data, 1, MpiTraits<T>::type, root); 
   }

   // C Array partial specializations

   /**
   * Send a C-array of T values
   *
   * Throws an exception if their exists neither an associated MPI
   * data type nor an explicit specialization of the scalar send<T>.
   *
   * \param comm   MPI communicator
   * \param array  address of first element in array
   * \param count  number of elements in array
   * \param dest   MPI rank of destination (receiving) processor in comm
   * \param tag    user-defined integer identifier for this message
   */
   template <typename T>
   void send(MPI::Comm& comm, T* array, int count, int dest, int tag)
   {  
      if (MpiTraits<T>::hasType) {
         comm.Send(array, count, MpiTraits<T>::type, dest, tag); 
      } else { 
         // Try send<T> by element, in case of explicit specialization.
         // If there is no specialization or type, send<T> throws.
         for (int i = 0; i < count; ++i) {
            send<T>(comm, array[i], dest, tag); 
         }
      }
   } 

   /**
   * Receive a C-array of T objects.
   *
   * Throws an exception if their exists neither an associated MPI
   * data type nor an explicit specialization of the scalar recv<T>.
   *
   * \param comm   MPI communicator
   * \param array  address of first element in array
   * \param count  number of elements in array
   * \param source MPI rank of source (sending) processor in comm
   * \param tag    user-defined integer identifier for this message
   */
   template <typename T>
   void recv(MPI::Comm& comm, T* array, int count, int source, int tag)
   {  
      if (MpiTraits<T>::hasType) {
         comm.Recv(array, count, MpiTraits<T>::type, source, tag); 
      } else {
         // Try recv<T> by element, in case of explicit specialization.
         // If there is no specialization or type, recv<T> throws.
         for (int i = 0; i < count; ++i) {
            recv<T>(comm, array[i], source, tag); 
         }
      }
   }

   /**
   * Broadcast a C-array of T objects.
   *
   * Throws an exception if their exists neither an associated MPI
   * data type nor an explicit specialization of the scalar bcast<T>.
   *
   * \param comm   MPI communicator
   * \param array  address of first element in array
   * \param count  number of elements in array
   * \param root   MPI rank of root (sending) processor in comm
   */
   template <typename T>
   void bcast(MPI::Intracomm& comm, T* array, int count, int root)
   {  
      if (MpiTraits<T>::hasType) {
         comm.Bcast(array, count, MpiTraits<T>::type, root); 
      } else {
         // Try bcast<T> by element, in case of explicit specialization.
         // If there is no specialization or type, bcast<T> throws.
         for (int i = 0; i < count; ++i) {
            bcast<T>(comm, array[i], root); 
         }
      }
   }

   // DArray container partial specializations

   /**
   * Send a DArray<T> container.
   *
   * Throws an exception if their exists neither an associated MPI data 
   * type nor an explicit specialization of the scalar send<T> method.
   *
   * \param comm   MPI communicator
   * \param array  DArray object
   * \param count  logical number of elements in array
   * \param dest   MPI rank of destination (receiving) processor in comm
   * \param tag    user-defined integer identifier for this message
   */
   template <typename T>
   void send(MPI::Comm& comm, DArray<T>& array, int count, int dest, int tag)
   {  
      // Preconditions
      if (!(array.isAllocated())) {
         UTIL_THROW("Cannot read unallocated DArray");
      }
      if (count > array.capacity()) {
         UTIL_THROW("Error: Logical size count > DArray capacity");
      }

      if (MpiTraits<T>::hasType) {
         comm.Send(&array[0], count, MpiTraits<T>::type, dest, tag); 
      } else {
         // Try send<T> by element, in case of explicit specialization.
         // If there is no specialization or type, send<T> throws.
         for (int i = 0; i < count; ++i) {
            send<T>(comm, array[i], dest, tag); 
         }
      }
   }
  
   /**
   * Receive a DArray<T> container.
   *
   * Throws an exception if their exists neither an associated MPI data
   * type nor an explicit specialization of the scalar recv<T> method.
   *
   * \param comm   MPI communicator
   * \param array  DArray object
   * \param count  logical number of elements in array
   * \param source MPI rank of source (sending) processor in comm
   * \param tag    user-defined integer identifier for this message
   */
   template <typename T>
   void recv(MPI::Comm& comm, DArray<T>& array, int count, int source, int tag)
   {  
      // Preconditions
      if (!(array.isAllocated())) {
         UTIL_THROW("Cannot read unallocated DArray");
      }
      if (count > array.capacity()) {
         UTIL_THROW("Error: Logical size count > DArray capacity");
      }

      if (MpiTraits<T>::hasType) {
         comm.Recv(&array[0], count, MpiTraits<T>::type, source, tag); 
      } else {
         // Try recv<T> by element, in case of explicit specialization.
         // If there is no specialization or type, recv<T> throws.
         for (int i = 0; i < count; ++i) {
            recv<T>(comm, array[i], source, tag); 
         }
      }
   }
  
   /**
   * Broadcast a DArray<T> container.
   *
   * Throws an exception if their exists neither an associated MPI
   * data type nor an explicit specialization of the scalar bcast<T>.
   *
   * \param comm   MPI communicator
   * \param array  address of first element in array
   * \param count  number of elements in array
   * \param root   MPI rank of root (sending) processor in comm
   */
   template <typename T>
   void bcast(MPI::Intracomm& comm, DArray<T>& array, int count, int root)
   {  
      // Preconditions
      if (!(array.isAllocated())) {
         UTIL_THROW("Cannot read unallocated DArray");
      }
      if (count > array.capacity()) {
         UTIL_THROW("Error: Logical size count > DArray capacity");
      }

      if (MpiTraits<T>::hasType) {
         comm.Bcast(&array[0], count, MpiTraits<T>::type, root); 
      } else {
         // try bcast<T> by element, in case of explicit specialization.
         // If there is no specialization or type, bcast<T> throws.
         for (int i = 0; i < count; ++i) {
            bcast<T>(comm, array[i], root); 
         }
      }
   }
  
   // DMatrix container partial specializations

   /**
   * Send a DMatrix<T> container.
   *
   * Throws an exception if their exists neither an associated MPI
   * data type nor an explicit specialization of the scalar send<T>.
   *
   * \param comm    MPI communicator
   * \param matrix  DMatrix object to send
   * \param m       logical number of rows in matrix
   * \param n       logical number of columns in matrix
   * \param dest    MPI rank of destination (receiving) processor in comm
   * \param tag     user-defined integer identifier for this message
   */
   template <typename T>
   void send(MPI::Comm& comm, DMatrix<T>& matrix, int m, int n, int dest, int tag)
   {  
      // Preconditions
      if (!(matrix.isAllocated())) {
         UTIL_THROW("Cannot read unallocated DMatrix");
      }
      if (m > matrix.capacity1()) {
         UTIL_THROW("Error: Logical size m > DMatrix<T>::capacity1()");
      }
      if (n > matrix.capacity2()) {
         UTIL_THROW("Error: Logical size n > DMatrix<T>::capacity2()");
      }

      if (MpiTraits<T>::hasType) {
         int mp = matrix.capacity1();
         int np = matrix.capacity2();
         comm.Send(&matrix(0, 0), mp*np, MpiTraits<T>::type, dest, tag); 
         // Note: This method sends the entire physical memory block.
      } else {
         // try send<T> by element, in case of explicit specialization.
         // If there is no specialization or type, send<T> throws.
         int i, j;
         for (i = 0; i < m; ++i) {
            for (j = 0; j < n; ++j) {
               send<T>(comm, matrix(i, j), dest, tag); 
            }
         }
      }
   }
  
   /**
   * Receive a DMatrix<T> container.
   *
   * Throws an exception if their exists neither an associated MPI
   * data type nor an explicit specialization of the scalar recv<T>.
   *
   * \param comm    MPI communicator
   * \param matrix  DMatrix object to receive
   * \param m       logical number of rows in matrix
   * \param n       logical number of columns in matrix
   * \param source  MPI rank of source (sending) processor in comm
   * \param tag     user-defined integer identifier for this message
   */
   template <typename T>
   void recv(MPI::Comm& comm, DMatrix<T>& matrix, int m, int n, 
             int source, int tag)
   {  
      // Preconditions
      if (!(matrix.isAllocated())) {
         UTIL_THROW("Cannot recv unallocated DMatrix");
      }
      if (m > matrix.capacity1()) {
         UTIL_THROW("Error: Logical size m > DMatrix<T>::capacity1()");
      }
      if (n > matrix.capacity2()) {
         UTIL_THROW("Error: Logical size n > DMatrix<T>::capacity2()");
      }

      if (MpiTraits<T>::hasType) {
         int mp = matrix.capacity1();
         int np = matrix.capacity2();
         comm.Recv(&matrix(0, 0), mp*np, MpiTraits<T>::type, source, tag); 
         // Note: This method receives the entire physical memory block.
      } else {
         // try recv<T> by element, in case of explicit specialization.
         // If there is no specialization or type, recv<T> throws.
         int i, j;
         for (i = 0; i < m; ++i) {
            for (j = 0; j < n; ++j) {
               recv<T>(comm, matrix(i, j), source, tag); 
            }
         }
      }
   }
  
   /**
   * Broadcast a DMatrix<T> container.
   *
   * Throws an exception if their exists neither an associated MPI
   * data type nor an explicit specialization of the scalar bcast<T>.
   *
   * \param comm    MPI communicator
   * \param matrix  DMatrix object
   * \param m       logical number of rows in matrix
   * \param n       logical number of columns in matrix
   * \param root    MPI rank of root (sending) processor in comm
   */
   template <typename T>
   void bcast(MPI::Intracomm& comm, DMatrix<T>& matrix, int m, int n, int root)
   {  
      // Preconditions
      if (!(matrix.isAllocated())) {
         UTIL_THROW("Cannot bcast unallocated DMatrix");
      }
      if (m > matrix.capacity1()) {
         UTIL_THROW("Error: Logical size m > DMatrix<T>::capacity1()");
      }
      if (n > matrix.capacity2()) {
         UTIL_THROW("Error: Logical size n > DMatrix<T>::capacity2()");
      }

      if (MpiTraits<T>::hasType) {
         int mp = matrix.capacity1();
         int np = matrix.capacity2();
         comm.Bcast(&matrix(0, 0), mp*np, MpiTraits<T>::type, root); 
         // Note: This method receives the entire physical memory block.
      } else {
         // Try bcast<T> by element, in case of explicit specialization.
         // If there is no specialization or type, bcast<T> throws.
         int i, j;
         for (i = 0; i < m; ++i) {
            for (j = 0; j < n; ++j) {
               bcast<T>(comm, matrix(i, j), root); 
            }
         }
      }
   }

   // bool (explicit specializations)

   /**
   * Explicit specialization of send for bool data.
   */ 
   template <>
   void send<bool>(MPI::Comm& comm, bool& data, int dest, int tag);

   /**
   * Explicit specialization of recv for bool data.
   */ 
   template <>
   void recv<bool>(MPI::Comm& comm, bool& data, int source, int tag);

   /**
   * Explicit specialization of bcast for bool data.
   */ 
   template <>
   void bcast<bool>(MPI::Intracomm& comm, bool& data, int root);

   // std::string (explicit specializations)
 
   /**
   * Explicit specialization of send for std::string data.
   */ 
   template <> void 
   send<std::string>(MPI::Comm& comm, std::string& data, int dest, int tag);

   /**
   * Explicit specialization of recv for std::string data.
   */ 
   template <> void 
   recv<std::string>(MPI::Comm& comm, std::string& data, int source, int tag);

   /**
   * Explicit specialization of bcast for std::string data.
   */ 
   template <>
   void bcast<std::string>(MPI::Intracomm& comm, std::string& data, int root);

}
#endif
#endif
