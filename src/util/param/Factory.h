#ifndef UTIL_FACTORY_H
#define UTIL_FACTORY_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/mpi/MpiFileIo.h>         // member
#include <util/param/ParamComposite.h>  // used in implementation
#include <util/param/Begin.h>           // used in implementation
#include <util/param/End.h>             // used in implementation
#include <util/global.h>

#ifdef UTIL_MPI
#include <util/mpi/MpiSendRecv.h>
#endif

#include <string>

namespace Util
{

   /**
   * Factory template.
   *
   * \ingroup Manager_Module
   */
   template <typename Data>
   class Factory
   {

   public:

      /**
      * Constructor.
      */
      Factory();

      /**
      * Destructor.
      */
      virtual ~Factory();

      /**
      * Add a new subfactory to the list.
      *
      * \param subfactory New subfactory to be added
      */
      void addSubfactory(Factory<Data>& subfactory);

      /**
      * Returns a pointer to a new instance of specified subclass.
      *
      * This method takes the name className of a subclass of Data as a
      * parameter, and attempts to instantiate an object of that class.
      * If it recognizes the className, it creates an instance of that
      * class and returns a Data* base class pointer to the new object.
      * If it does not recognize the className, it returns a null pointer.
      *
      * An implementation should first call "trySubfactories(className)"
      * and immediately return if this returns a non-null pointer, before
      * attempting to match the className against specific strings.
      *
      * \param  className name of subclass
      * \return base class pointer to new object, or a null pointer.
      */
      virtual Data* factory(const std::string &className) const = 0;

      /**
      * Read a class name, instantiate an object, and read its parameters.
      *
      * This method:
      *  - reads a comment line of the form className + {
      *  - invokes the factory method to create an instance of className
      *  - invokes the readParam() method of the new object
      *
      * When compiled with MPI, if the parent ParamComposite has a param
      * communicator, this method reads the comment line on the Io processor,
      * broadcasts it to all others, and then lets each processor 
      * independently match this string. 
      *
      * \throws Exception if className is not recognized.
      *
      * \param  in        input stream
      * \param  parent    parent ParamComposite object
      * \param  className (output) name of subclass of Data
      * \param  isEnd     (output) is the input a closing bracket "}" ?
      * \return pointer to new instance of className
      */
      Data* readObject(std::istream &in, ParamComposite& parent,
                       std::string& className, bool& isEnd);

      /**
      * Load a class name, instantiate an object, and load the object.
      *
      * This method:
      *  - loads a className from an input archive
      *  - invokes the factory method to create an instance of className
      *  - invokes the load() method of the new object
      *
      * When compiled with MPI, if the parent ParamComposite has a param
      * communicator, this method loads the comment line on the Io processor,
      * broadcasts it to all others, and then lets each processor 
      * independently match this string. 
      *
      * \throws Exception if className is not recognized.
      *
      * \param  ar        input/loading archive
      * \param  parent    parent ParamComposite object
      * \param  className (output) name of subclass of Data
      * \return pointer to new instance of className
      */
      Data* loadObject(Serializable::IArchive& ar, ParamComposite& parent,
                       std::string& className);

   protected:

      /**
      * Search through subfactories for match.
      *
      * This method iterates through all registered subfactories, calls
      * the factory(const std::string& ) method of each, and immediately
      * returns a pointer to a new object if any of them returns a non-null
      * pointer. If all of them return a null pointer, this method also
      * returns a null pointer.
      *
      * \param  className name of subclass
      * \return base class pointer to new object, or a null pointer.
      */
      Data* trySubfactories(const std::string& className) const;

      #ifdef UTIL_MPI
      /**
      * Set associated Mpi communicator.
      *
      * Is not recursive (is not applied to subfactories).
      *
      * \param communicator MPI Intra-communicator to use for input
      */
      void setIoCommunicator(MPI::Intracomm& communicator);

      /**
      * Does this factory have a param communicator?
      */
      bool hasIoCommunicator() const;
      #endif

   private:

      /// Vector of pointers to child factories.
      std::vector< Factory<Data>* > subfactories_;

      /// Object to identify if this processor can do file Io.
      MpiFileIo   paramFileIo_;

   };

   /*
   * Constructor.
   */
   template <typename Data>
   Factory<Data>::Factory()
    : paramFileIo_()
   {}

   /*
   * Destructor.
   */
   template <typename Data>
   Factory<Data>::~Factory()
   {}

   #ifdef UTIL_MPI
   /*
   * Set the param communicator.
   */
   template <typename Data>
   void Factory<Data>::setIoCommunicator(MPI::Intracomm& communicator)
   {
      if (paramFileIo_.hasIoCommunicator()) {
         if (&paramFileIo_.ioCommunicator() != &communicator) {
            UTIL_THROW("Attempt to modify Factory param communicator");
         }
      } else {
         paramFileIo_.setIoCommunicator(communicator);
      }
   }

   /*
   * Does thus factory have a param communicator?
   */
   template <typename Data>
   bool Factory<Data>::hasIoCommunicator() const
   {  return paramFileIo_.hasIoCommunicator(); }
   #endif

   /*
   * Add a subfactory to the list of children.
   */
   template <typename Data>
   void Factory<Data>::addSubfactory(Factory<Data>& subfactory)
   {  subfactories_.push_back(&subfactory); }




   /*
   * Read subclass name, create object, and read its parameters.
   */
   template <typename Data>
   Data* Factory<Data>::readObject(std::istream &in, ParamComposite& parent,
                                   std::string& className, bool& isEnd)
   {
      std::string  commentString;
      Data*        typePtr = 0;
      int          length;
      bool         hasData = false; // initialized to avoid compiler warning

      #ifdef UTIL_MPI
      // Set ioCommunicator to that of parent, if any.
      if (parent.hasIoCommunicator()) {
         setIoCommunicator(parent.ioCommunicator());
      }
      #endif

      // Read a first line of the form "ClassName{"
      if (paramFileIo_.isIoProcessor()) {
         in >> commentString;
      }
      #ifdef UTIL_MPI
      // Broadcast the full string to all processors.
      if (paramFileIo_.hasIoCommunicator()) {
         bcast<std::string>(paramFileIo_.ioCommunicator(), commentString, 0);
      }
      // Hereafter, each processor independently processes the same string.
      #endif
      length = commentString.size();

      // If commentString = '}', set isEnd=true and return null ptr.
      if (length == 1 && commentString[0] == '}') {
         className = std::string();
         isEnd = true;
         return 0; 
      } else {
         isEnd = false;
      }

      // Isolate className by stripping the trailing "{" bracket
      if (commentString[length-1] == '{') {
         className = commentString.substr(0, commentString.size() - 1);
         hasData = true;
      } else
      if (commentString[length-1] == '}') {
         className = commentString.substr(0, commentString.size() - 2);
         hasData = false;
      } else {
         if (paramFileIo_.isIoProcessor()) {
            Log::file() << "commentString = " << commentString << std::endl;
            Log::file() << "className     = " << className     << std::endl;
            UTIL_THROW("Invalid first line\n");
         }
      }

      // Create new object of the specified class.
      typePtr = factory(className);

      // If the subclass name was recognized:
      if (typePtr) {

         // Add Begin object to new child ParamComposite, indented for child.
         Begin* beginPtr;
         beginPtr = &typePtr->addBegin(className.c_str());
         beginPtr->setIndent(parent);

         // Echo Begin object, if appropriate
         if (ParamComponent::echo() && paramFileIo_.isIoProcessor()) {
            beginPtr->writeParam(Log::file());
         }

         // Read parameters for the new child object, if any
         if (hasData) {
            parent.addParamComposite(*typePtr);
            typePtr->readParameters(in);
         }

         // Read closing bracket, set indentation as for child.
         typePtr->readEnd(in).setIndent(parent);

         // Note: The readParameters() methods for managed objects should 
         // not read begin and end lines, which read here. 

      }
      return typePtr;
   }

   /*
   * Load subclass name, create object, and load object.
   */
   template <typename Data>
   Data* Factory<Data>::loadObject(Serializable::IArchive& ar, ParamComposite& parent,
                                   std::string& className)
   {
      #ifdef UTIL_MPI
      // Set ioCommunicator to that of parent, if any.
      if (parent.hasIoCommunicator()) {
         setIoCommunicator(parent.ioCommunicator());
      }
      #endif

      // Read the class name.
      if (paramFileIo_.isIoProcessor()) {
         ar & className;
      }

      #ifdef UTIL_MPI
      // Broadcast the full string to all processors.
      if (paramFileIo_.hasIoCommunicator()) {
         bcast<std::string>(paramFileIo_.ioCommunicator(), className, 0);
      }
      #endif

      // Create and load a new object of the specified class.
      Data* typePtr = factory(className);
      if (typePtr) {
         parent.loadParamComposite(ar, *typePtr);
      } else {
         Log::file() << "Failed attempt to create instance of " 
                     << className << std::endl;
      }
      return typePtr;
   }

   /*
   * Try all subfactories in sequence searching for a match.
   */
   template <typename Data>
   Data* Factory<Data>::trySubfactories(const std::string& className) const
   {
      Data* typePtr = 0;
      int n = subfactories_.size();
      for (int i = 0; i < n && typePtr == 0; ++i) {
         typePtr = subfactories_[i]->factory(className);
      }
      return typePtr;
   }

}
#endif
