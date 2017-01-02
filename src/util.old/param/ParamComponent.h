#ifndef UTIL_PARAM_COMPONENT_H
#define UTIL_PARAM_COMPONENT_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/mpi/MpiFileIo.h>          // base class
#include <util/archives/Serializable.h>  // base class (typedef)
#include <util/global.h>

#include <iostream>
#include <sstream>
#include <string>

namespace Util
{

   /**
   *  Abstract base class for classes that input and ouput parameters to file.
   *
   *  The readParam  method reads a parameter or parameter list from iostream. 
   *  The writeParam method writes a parameter or parameter list to an ostream. 
   *  The same io format should be used by write and read methods. 
   *
   *  \ingroup Param_Module
   */
   class ParamComponent : public Serializable, public MpiFileIo
   {

   public:

      /**
      * Destructor.
      */
      virtual ~ParamComponent();

      /**
      * Read parameter(s) from file.
      *
      * \param in input stream
      */
      virtual void readParam(std::istream& in) = 0;

      /**
      * Read parameter(s) to file.
      *
      * \param out output stream
      */
      virtual void writeParam(std::ostream& out) = 0;

      /**
      * Load internal state from an archive.
      *
      * The default implementation is empty. This default is used by
      * the Begin, End, and Blank subclasses.
      */ 
      virtual void load(Serializable::IArchive &ar)
      {}

      /**
      * Save internal state to an archive.
      *
      * The default implementation is empty. This default is used by
      * the Begin, End, and Blank subclasses.
      */ 
      virtual void save(Serializable::OArchive &ar)
      {}

      /**
      * Nontrivial implementation provided by ParamComposite subclass.
      *
      * The default implementation is empty. This default is used by
      * all leaf nodes (all other than ParamComposite and subclasses).
      */
      virtual void resetParam()
      {}

      /**
      * Set indent level.
      *
      * If next=true (default) set indent level one higher than that
      * of parent. If next=false, set indent level the same as parent.
      *
      * \param parent parent ParamComponent object
      * \param next   If true, set level one higher than for parent.
      */
      void setIndent(const ParamComponent& parent, bool next = true);

      /**
      * Return indent string for this object (string of spaces).
      */
      std::string indent() const;

      /**
      * Serialize this ParamComponent as a string.
      *
      * \param ar      saving or loading archive
      * \param version version id for archive
      */
      template <class Archive>
      void serialize(Archive& ar, const unsigned int version);

      // Public static methods

      /**
      * Initialize static echo member to false.
      */
      static void initStatic();

      /**
      * Enable or disable echoing for all subclasses of ParamComponent.
      *
      * When echoing is enabled, all parameters are echoed to a log
      * file immediately after being read. This is useful as an aid
      * to debugging the parameter file, by showing where the error
      * occurred.
      *
      * \param echo set true to enable echoing, false to disable. 
      */
      static void setEcho(bool echo = true);

      /**
      * Get echo parameter. 
      *
      * \return true if echoing is enabled, false if disabled.
      */
      static bool echo();

   protected:

      /**
      * Constructor.
      *
      * Protected to prevent instantiation of a conceptually abstract 
      * class. 
      *
      * On return the indent string is empty. If UTIL_MPI is defined,
      * no communicator is set upon construction.
      */
      ParamComponent();

      /**
      * Copy constructor.
      */
      ParamComponent(const ParamComponent& other);

   private:

      /// Indentation string, a string of spaces.
      std::string indent_;

      /// Parameter to enable (true) or disable (false) echoing.
      static bool echo_;

   // friend:

      #ifdef UTIL_MPI
      template <class T> friend class Factory;

      /*
      * Note: This friend declaration allows the Factory::readObject() 
      * method to set the Factory ioCommunicator to be that of the 
      * parent ParamComposite. 
      */
      #endif

   };

   /*
   * Serialize a ParamComponent as a string.
   */
   template <class Archive>
   void ParamComponent::serialize(Archive& ar, const unsigned int version)
   {
      std::string str;
      if (Archive::is_saving()) {
         std::ostringstream buffer;
         writeParam(buffer);
         str = buffer.str();
      } 
      ar & str;
      if (Archive::is_loading()) {
         std::istringstream buffer(str);
         readParam(buffer);
      } 
   }

} 
#endif
