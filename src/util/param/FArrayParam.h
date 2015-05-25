#ifndef UTIL_F_ARRAY_PARAM_H
#define UTIL_F_ARRAY_PARAM_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/Parameter.h>
#include <util/containers/FArray.h>
#include <util/global.h>

#include <iomanip> 

namespace Util
{

   /**
   * A Parameter associated with a FArray container.
   *
   * \ingroup Param_Module
   */
   template <class Type, int N>
   class FArrayParam : public Parameter
   {
   
   public:

      /**
      * Constructor.
      *
      * \param label  label string for parameter file
      * \param array  associated FArray variable
      * \param isRequired  Is this a required parameter?
      */
      FArrayParam(const char *label, FArray<Type, N>& array, bool isRequired = true);
 
      /** 
      * Write FArray parameter to stream.
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
   
      /// Pointer to associated FArray.
      FArray<Type, N>* arrayPtr_;
   
   };

   /*
   * FArrayParam<Type, N> constructor.
   */
   template <class Type, int N>
   FArrayParam<Type, N>::FArrayParam(const char *label, FArray<Type, N>& array, bool isRequired)
    : Parameter(label, isRequired),
      arrayPtr_(&array)
   {}

   /*
   * Read a FArray from isteam.
   */
   template <class Type, int N>
   void FArrayParam<Type, N>::readValue(std::istream &in)
   {  
      for (int i = 0; i < N; ++i) {
         in >> (*arrayPtr_)[i];
      }
   }

   /*
   * Load a FArray from input archive.
   */
   template <class Type, int N>
   void FArrayParam<Type, N>::loadValue(Serializable::IArchive& ar)
   {  ar >> *arrayPtr_; }  

   /*
   * Save a FArray to an output archive.
   */
   template <class Type, int N>
   void FArrayParam<Type, N>::saveValue(Serializable::OArchive& ar)
   { ar << *arrayPtr_; }

   #ifdef UTIL_MPI
   /*
   * Broadcast a FArray.
   */
   template <class Type, int N>
   void FArrayParam<Type, N>::bcastValue()
   {  bcast<Type>(ioCommunicator(), &((*arrayPtr_)[0]), N, 0); }
   #endif

   /*
   * Write a FArray parameter.
   */
   template <class Type, int N>
   void FArrayParam<Type, N>::writeParam(std::ostream &out) 
   {
      if (isActive()) {
         Label space("");
         for (int i = 0; i < N; ++i) {
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
      }
   }

} 
#endif
