#ifndef UTIL_SCALAR_PARAM_H
#define UTIL_SCALAR_PARAM_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/Parameter.h>
#include <util/archives/Serializable_includes.h>
#include <util/global.h>

#ifdef UTIL_MPI
#include <util/mpi/MpiSendRecv.h>
#endif

#include <iomanip> 

namespace Util
{

   /** 
   * Template for a Parameter object associated with a scalar variable.
   *
   * This template can be used to define a Parameter subclass for any 
   * data type for which there exist inserter (<<) and extractor (>>) 
   * operators for stream io. 
   *
   * \ingroup Param_Module
   */
   template <class Type>
   class ScalarParam : public Parameter
   {
   
   public:
  
      /**
      * Constructor.
      *
      * \param label  label string const.
      * \param value  reference to parameter value.
      * \param isRequired  Is this a required parameter?
      */
      ScalarParam(const char *label, Type& value, bool isRequired = true); 

      /** 
      * Write parameter to stream.
      *
      * \param out output stream
      */
      void writeParam(std::ostream& out);

      /**
      * Set the pointer to point a specific variable.
      *
      * \param value variable that holds the parameter value. 
      */
      void setValue(Type& value);

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
      Type* valuePtr_;

      /// Private and not implemented to prevent copying.
      ScalarParam(const ScalarParam<Type>& other);

      /// Private and not implemented to prevent assignment.
      ScalarParam<Type> operator = (const ScalarParam<Type>& other);
   
   };

   // Member Function Definitions

   /*
   * ScalarParam<Type> constructor.
   */
   template <class Type>
   ScalarParam<Type>::ScalarParam(const char *label, Type &value, 
                                  bool isRequired)
    : Parameter(label, isRequired),
      valuePtr_(&value)
   {}

   template <class Type>
   void ScalarParam<Type>::readValue(std::istream &in)
   {  in >> *valuePtr_; }

   template <class Type>
   void ScalarParam<Type>::loadValue(Serializable::IArchive& ar)
   {  ar & *valuePtr_; }

   template <class Type>
   void ScalarParam<Type>::saveValue(Serializable::OArchive& ar)
   {  ar & *valuePtr_; }

   #ifdef UTIL_MPI
   template <class Type>
   void ScalarParam<Type>::bcastValue()
   {  bcast<Type>(ioCommunicator(), *valuePtr_, 0); }
   #endif

   /*
   * Write a parameter.
   */
   template <class Type>
   void ScalarParam<Type>::writeParam(std::ostream& out)
   {
      if (isActive()) {
         out << indent();
         out << label_;
         out << std::right << std::scientific 
             << std::setprecision(Parameter::Precision) 
             << std::setw(Parameter::Width) << *valuePtr_ 
             << std::endl;
      }
   }

   /*
   * Set the pointer to the parameter value.
   */
   template <class Type>
   void ScalarParam<Type>::setValue(Type& value)
   {  valuePtr_ = &value; }

} 
#endif
