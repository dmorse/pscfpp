#ifndef UTIL_PARAMETER_H
#define UTIL_PARAMETER_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/ParamComponent.h>
#include <util/archives/Serializable_includes.h>
#include <util/param/Label.h>
#include <util/global.h>

namespace Util
{

   /**
   * A single variable in a parameter file.
   *
   * Parameter is a base class for objects that read and write the value
   * of a single C++ variable from or to a parameter file. The parameter
   * file format for a parameter contains a string label followed by a 
   * value for the variable. Different subclasses of parameter are 
   * specialized for different variable types, which can include primitive
   * C/C++ variables, user defined types that overload the << and >> 
   * operators, or any of several different types of container.
   *
   * A Parameter may be required or optional element in a parameter file,
   * depending on the value of the bool isRequired parameter of the
   * constructor. An optional element becomes "active" when an entry
   * with the correct label is read from a parameter file, or when an
   * active value is loaded from an archive. By convention, a required 
   * Parameter is always active, even before its value is read or loaded. 
   * The bool functions isRequired() and isActive() can be used to query
   * the state of a Parameter.
   *
   * The overloaded saveOptional() static member functions can be used
   * to save optional parameters to an archive in a form that records
   * whether or not they are active.
   *
   * \ingroup Param_Module
   */
   class Parameter : public ParamComponent
   {

   public:

      // Static members

      /// Width of output field for a scalar variable.
      static const int Width = 20;

      /// Precision for io of floating point data field.
      static const int Precision = 12;

      /**
      * Save an optional parameter value to an output archive
      *
      * \param ar  output archive to which to save
      * \param value  reference to value of optional parameter
      * \param isActive  Is this parameter present in the parameter file?
      */
      template <class Type>
      static void
      saveOptional(Serializable::OArchive& ar, Type& value, bool isActive);

      /**
      * Save an optional C-array of n values to an output archive
      *
      * \param ar  output archive to which to save
      * \param ptr  pointer to first element of optional C-array parameter
      * \param n  number of elements in array
      * \param isActive  Is this parameter present in the parameter file?
      */
      template <class Type>
      static void
      saveOptionalCArray(Serializable::OArchive& ar, Type* ptr, int n,
                         bool isActive);

      /**
      * Save an optional two-dimensional C array to an output archive.
      *
      * \param ar  output archive to which to save
      * \param ptr  pointer to first element optional 2D C-array parameter
      * \param m  logical number of rows in array
      * \param n  logical number of columns in array
      * \param np  logical number of columns in array
      * \param isActive  Is this parameter present in the parameter file?
      */
      template <class Type>
      static void
      saveOptionalCArray2D(Serializable::OArchive& ar, Type* ptr, 
                           int m, int n, int np, bool isActive);

      // Non-static member functions

      /**
      * Constructor.
      *
      * \param label  label string preceding value in file format
      * \param isRequired  Is this a required parameter?
      */
      Parameter(const char *label, bool isRequired = true);

      /**
      * Destructor.
      */
      virtual ~Parameter();

      /**
      * Read a label and (if the label matches) a parameter value.
      *
      * The parameter file format for a Parameter consists of a label string
      * followed by value. The value is read if and only if the label matches
      * the expected value for this Parameter. If this Parameter is required
      * and the input label not match, an error message is printed to the log
      * file and Exception is thrown. If the Parameter is not required and the
      * input label does not match, the label string is retained in an buffer
      * for later processing by the readParam method of other ParamComponent
      * objects.
      *
      * Upon entry to this function, a label string is read into a label buffer
      * if and only if the buffer is empty. This buffer is a static member of
      * the Label class, which can retain a label between invocations of the
      * readParameter method of different ParamComponent objects.  Once a
      * label string is read from file, it remains in the label buffer until
      * until it is matched, at which point the buffer is cleared to allow
      * processing of the next label.
      *
      * \param in  input stream from which to read
      */
      virtual void readParam(std::istream &in);

      /**
      * Load from an archive.
      *
      * An optional Parameter loads the value of an isActive flag, and then
      * loads the parameter value only if the isActive is true. A required
      * Parameter simply loads the parameter value. The variable associated
      * with an optional Parameter must be set to its default value before
      * attempting to load the parameter. Optional parameters should be
      * saved either using the save() method of an associated Parameter
      * object or using the appropriate overloaded Parameter::saveOptional() 
      * static member function, which both use the required format.
      *
      * \param ar input archive from which to load
      */
      virtual void load(Serializable::IArchive& ar);

      /**
      * Save to an archive.
      *
      * An optional Parameter saves the value of the isActive flag, and then
      * saves a parameter value only if the isActive is true. A required
      * Parameter simply saves its value. The label string is not saved to
      * the archive.
      *
      * The overloaded static saveOptional functions can also be used to 
      * save optional parameter values in this format.
      *
      * \param ar output archive to which to save
      */
      virtual void save(Serializable::OArchive& ar);

      /**
      * Return label string.
      */
      std::string label() const;

      /**
      * Is this an optional parameter?
      */
      bool isRequired() const;

      /**
      * Is this parameter active?
      */
      bool isActive() const;

   protected:

      /// Label object that contains parameter label string.
      Label label_;

      /// Is this parameter active (always true if isRequired).
      bool isActive_;

      /**
      * Read parameter value from an input stream.
      *
      * \param in input stream from which to read
      */
      virtual void readValue(std::istream& in){}

      /**
      * Load bare parameter value from an archive.
      *
      * \param ar input archive from which to load
      */
      virtual void loadValue(Serializable::IArchive& ar){}

      /**
      * Save parameter value to an archive.
      *
      * \param ar output archive to which to save
      */
      virtual void saveValue(Serializable::OArchive& ar){}

      #ifdef UTIL_MPI
      /**
      * Broadcast parameter value within the ioCommunicator.
      */
      virtual void bcastValue(){}
      #endif

   };

   /*
   * Save a parameter value to an output archive
   */
   template <class Type>
   void Parameter::saveOptional(Serializable::OArchive& ar, 
                                Type& value, bool isActive)
   {
      ar << isActive;
      if (isActive) {
         ar & value;
      }
   }

   /*
   * Save C-array of n ptrs to an output archive
   */
   template <class Type>
   void Parameter::saveOptionalCArray(Serializable::OArchive& ar, 
                                      Type* ptr, int n, bool isActive)
   {
      ar << isActive;
      if (isActive) {
         ar.pack(ptr, n);
      }
   }

   /*
   * Save a two-dimensional C array to an output archive.
   */
   template <class Type> 
   void Parameter::saveOptionalCArray2D(Serializable::OArchive& ar, 
                                        Type* ptr, int m, int n, int np, 
                                        bool isActive)
   {
      ar << isActive;
      if (isActive) {
         ar.pack(ptr, m, n, np);
      }
   }

}
#endif
