/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "ParamComposite.h"   // class header
#include "Begin.h"
#include "End.h"
#include "Blank.h"
#include <util/global.h>

#include <cstdio>
#include <cstring>

namespace Util
{

   /*
   * Default Constructor.
   */
   ParamComposite::ParamComposite()
    : ParamComponent(),
      list_(),
      isLeaf_(),
      size_(0),
      className_("ParamComposite"),
      isRequired_(true),
      isActive_(true)
   {}

   /*
   * Constructor.
   */
   ParamComposite::ParamComposite(int capacity)
    : ParamComponent(),
      list_(),
      isLeaf_(),
      size_(0),
      className_("ParamComposite"),
      isRequired_(true),
      isActive_(true)
   {
      if (capacity <= 0 ) {
         UTIL_THROW("Attempt to reserve capacity <= 0");
      }
      list_.reserve(capacity);
      isLeaf_.reserve(capacity);
   }

   /*
   * Copy constructor.
   */
   ParamComposite::ParamComposite(const ParamComposite& other)
    : ParamComponent(other),
      list_(),
      isLeaf_(),
      size_(0),
      className_(other.className_),
      isRequired_(true),
      isActive_(true)
   {}

   /*
   * Destructor.
   */
   ParamComposite::~ParamComposite()
   {
      if (size_ > 0) {
         for (int i=0; i < size_; ++i) {
            if (isLeaf_[i]) {
               delete list_[i];
            }
            /* Only delete Parameter, Begin, End & Blank leaf objects.
            * Do NOT delete node objects here. These are instances of
            * of subclasses of ParamComposite that are never created
            * by this object, and so should not be destroyed by this
            * object.
            */
         }
      }
   }

   /*
   * Read required parameter block, including begin and end.
   */
   void ParamComposite::readParam(std::istream &in)
   {
      assert(className_.size() > 0);
      isRequired_ = true;
      isActive_ = true;
      readBegin(in, className_.c_str());
      readParameters(in);
      readEnd(in);
   }

   /*
   * Read optional parameter block, including begin and end.
   */
   void ParamComposite::readParamOptional(std::istream &in)
   {
      assert(className_.size() > 0);
      isRequired_ = false;
      isActive_ = false;
      Begin* beginPtr = &readBegin(in, className_.c_str(), isRequired_);
      if (beginPtr->isActive()) {
         readParameters(in);
         readEnd(in);
         isActive_ = true;
      } else {
         delete beginPtr;
      }
   }

   /*
   * Default writeParam implementation.
   */
   void ParamComposite::writeParam(std::ostream &out)
   {
      if (isActive_) {
         for (int i=0; i < size_; ++i) {
            list_[i]->writeParam(out);
         }
      }
   }

   /*
   * Default load implementation, adds begin and end.
   */
   void ParamComposite::load(Serializable::IArchive& ar)
   {
      assert(className_.size() > 0);
      Begin* beginPtr = &addBegin(className_.c_str());
      if (ParamComponent::echo()) {
         if (isIoProcessor()) {
            beginPtr->writeParam(Log::file());
         }
      }
      loadParameters(ar);
      End* endPtr = &addEnd();
      if (ParamComponent::echo()) {
         if (isIoProcessor()) {
            endPtr->writeParam(Log::file());
         }
      }
   }

   /*
   * Optionally load this object.
   */
   void ParamComposite::loadOptional(Serializable::IArchive& ar)
   {
      // Load and broadcast isActive_ flag
      if (isIoProcessor()) {
         ar & isActive_;
         if (!isActive_) {
            if (ParamComponent::echo()) {
               Log::file() << indent() 
                           << className() << "{ [absent] }"
                           << std::endl; 
            }
         }
      } else {
         #ifdef UTIL_MPI
         if (!hasIoCommunicator()) {
            UTIL_THROW("Error: not isIoProcessor and not hasIoCommunicator");
         }
         #else
         UTIL_THROW("Error: not isIoProcessor and no MPI");
         #endif
      }
      #ifdef UTIL_MPI
      if (hasIoCommunicator()) {
         bcast<bool>(ioCommunicator(), isActive_, 0); 
      }
      #endif

      // Load object data iff isActive_
      if (isActive_) {
         load(ar);
      }
   }

   /*
   * Default save implementation.
   */
   void ParamComposite::save(Serializable::OArchive& ar)
   {
      for (int i=0; i < size_; ++i) {
         list_[i]->save(ar);
      }
   }

   /*
   * Save this optional ParamComposite.
   */
   void ParamComposite::saveOptional(Serializable::OArchive& ar)
   {
      ar & isActive_;
      if (isActive_) {
         save(ar);
      }
   }

   /*
   * Reset list to empty state.
   */
   void ParamComposite::resetParam()
   {
      for (int i=0; i < size_; ++i) {
         if (isLeaf_[i]) {
            delete list_[i];
         } else {
            list_[i]->resetParam();
         }
      }
      size_ = 0;
   }

   /*
   * Set this to be the parent of a child component.
   */
   void ParamComposite::setParent(ParamComponent& param, bool next)
   {
      param.setIndent(*this, next);
      #ifdef UTIL_MPI
      if (hasIoCommunicator()) {
         param.setIoCommunicator(ioCommunicator());
      }
      #endif 
   }

   /*
   * Add a leaf ParamComponent to the format array.
   */
   void ParamComposite::addComponent(ParamComponent& param, bool isLeaf)
   {
      list_.push_back(&param);
      isLeaf_.push_back(isLeaf);
      ++size_;
   }

   // ParamComposite object

   /*
   * Add a required ParamComposite node.
   */
   void ParamComposite::addParamComposite(ParamComposite &child, bool next)
   {
      bool isLeaf = false;
      setParent(child, next);
      addComponent(child, isLeaf);
   }

   /*
   * Add a required ParamComposite node, and call its readParam() method.
   */
   void
   ParamComposite::readParamComposite(std::istream &in, ParamComposite &child, 
                                      bool next)
   {
      addParamComposite(child, next);
      child.readParam(in);
   }

   /*
   * Add an optional ParamComposite, and call its readParamOptional method.
   */
   void
   ParamComposite::readParamCompositeOptional(std::istream &in, 
                                              ParamComposite &child, bool next)
   {
      addParamComposite(child, next);
      child.readParamOptional(in);
   }

   /*
   * Add a required ParamComposite node and load its data from archive ar.
   */
   void
   ParamComposite::loadParamComposite(Serializable::IArchive &ar, 
                                      ParamComposite &child, bool next)
   {
      addParamComposite(child, next);
      child.load(ar);
   }

   /*
   * Add an optional ParamComposite node, and load data if isActive.
   */
   void
   ParamComposite::loadParamCompositeOptional(Serializable::IArchive &ar, 
                                      ParamComposite &child, bool next)
   {
      addParamComposite(child, next);
      child.loadOptional(ar);
   }

   // Begin

   /*
   * Create and add a new Begin object.
   */
   Begin& ParamComposite::addBegin(const char *label)
   {
      Begin* ptr = new Begin(label);
      setParent(*ptr, false);
      addComponent(*ptr);
      return *ptr;
   }

   /*
   * Read the opening line of a ParamComposite.
   */
   Begin& ParamComposite::readBegin(std::istream &in, const char *label,
                                    bool isRequired)
   {
      Begin* ptr = new Begin(label, isRequired);
      setParent(*ptr, false);
      ptr->readParam(in);
      if (ptr->isActive()) {
         addComponent(*ptr);
      }
      return *ptr;
   }

   // End

   /*
   * Create and add a new End object.
   */
   End& ParamComposite::addEnd()
   {
      End* ptr = new End();
      setParent(*ptr, false);
      addComponent(*ptr);
      return *ptr;
   }

   /*
   * Read the closing bracket of a ParamComposite.
   */
   End& ParamComposite::readEnd(std::istream &in)
   {
      End* ptr = &addEnd();
      ptr->readParam(in);
      return *ptr;
   }

   // Blank

   /*
   * Create and add a new Blank object (a blank line).
   */
   Blank& ParamComposite::addBlank()
   {
      Blank* ptr = new Blank();
      setParent(*ptr);
      addComponent(*ptr);
      return *ptr;
   }

   /*
   * Read a blank line.
   */
   Blank& ParamComposite::readBlank(std::istream &in)
   {
      Blank* ptr = &addBlank();
      ptr->readParam(in);
      return *ptr;
   }

   /*
   * Set the class name string.
   */
   void ParamComposite::setClassName(const char * className)
   {  className_ = className; }

   /*
   * Set or unset the isActive flag.
   */
   void ParamComposite::setIsRequired(bool isRequired)
   {  
      isRequired_ = isRequired; 
      if (isRequired_) {
         isActive_ = true;
      }
   }

   /*
   * Set or unset the isActive flag.
   */
   void ParamComposite::setIsActive(bool isActive)
   {
      if (isRequired_ && !isActive) {
         UTIL_THROW("Error: cannot be required but not active");
      }
      isActive_ = isActive; 
   }

}
