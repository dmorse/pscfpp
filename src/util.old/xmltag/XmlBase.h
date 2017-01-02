#ifndef UTIL_XML_BASE_H
#define UTIL_XML_BASE_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/global.h>
#include <string>

namespace Util
{

   /**
   * Base class for classes that parse XML markup tags.
   * 
   * \ingroup XmlTag_Module
   */
   class XmlBase 
   {
   public:

      /**
      * Constructor.
      */
      XmlBase();

      /**
      * Destructor.
      */
      ~XmlBase();

      /**
      * Initialize string and cursor.
      */
      void setString(const std::string& string, int cursor = 0);

      /**
      * Set cursor. String must already be set.
      */
      void setCursor(int cursor);

      /**
      * Skip leading white space, if any.
      */
      void skip();
   
      /**
      * Advance to the next character.
      */
      void next();

      /**
      * Return the associated string.
      */
      const std::string& string() const;

      /**
      * Return the index of the current character.
      */
      int cursor() const;

      /**
      * Return the current character.
      */
      int c() const;

      /**
      * Has the cursor reached the end of the string?
      */
      bool isEnd() const;
   
   private:
   
      std::string const * stringPtr_;
      int end_;
      int cursor_; 
      char c_;

   };

   /*
   * Skip leading white space, if any.
   */
   inline void XmlBase::skip()
   {
      if (cursor_ == end_) return;
      while (c_ == ' ') {
         ++cursor_; 
         if (cursor_ == end_) {
            c_ = '\0';
            return;
         }
         c_ = (*stringPtr_)[cursor_];
      }
   }

   /*
   * Advance to the next character.
   */
   inline void XmlBase::next()
   {
      if (cursor_ == end_) return;
      ++cursor_; 
      if (cursor_ != end_) {
         c_ = (*stringPtr_)[cursor_];
      } else {
         c_ = '\0';
      }
   }

   /*
   * Return the associated string.
   */
   inline const std::string& XmlBase::string() const
   {  return (*stringPtr_); }

   /**
   * Return the index of the current character.
   */
   inline int XmlBase::cursor() const
   {  return cursor_; }

   /**
   * Return the current character.
   */
   inline int XmlBase::c() const
   {  return c_; }

   /**
   * Has the cursor reached the end of the string?
   */
   inline bool XmlBase::isEnd() const
   {  return (cursor_ == end_); }
   
}
#endif
