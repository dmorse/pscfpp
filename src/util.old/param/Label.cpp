/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Label.h"

#include <iomanip>

namespace Util
{

   /// Define static variables.
   bool  Label::isClear_ = true;
   std::string  Label::input_;

   // Static member functions

   /*
   * Clear input buffer (static member function).
   */
   void Label::clear()
   {
      input_.clear(); 
      isClear_ = true;
   }

   /*
   * Is the input buffer clear? (static member function).
   */
   bool Label::isClear() 
   {  return isClear_; }

   // Non-static member functions

   /*
   * Constructor.
   */
   Label::Label(bool isRequired)
    : isRequired_(isRequired),
      string_()
   {}

   /*
   * Constructor.
   */
   Label::Label(const char* string, bool isRequired)
    : isRequired_(isRequired),
      string_(string)
   {}

   /*
   * Copy constructor.
   */
   Label::Label(const Label& other)
    : isRequired_(other.isRequired_),
      string_(other.string_)
   {}

   /*
   * Destructor.
   */
   Label::~Label()
   {}

   /*
   * Set the label string.
   */
   void Label::setString(std::string string) 
   {  string_ = string; }

   /*
   * Return label string.
   */
   std::string Label::string() const
   {  return string_; }

   /*
   * Extract a label from an input stream.
   */
   std::istream& operator>>(std::istream& in, Label label)
   {
      // If previous input value matched, read a new one.
      if (label.isClear_) {
         in >> label.input_;
         label.isClear_ = false;
      }
      if (label.input_ == label.string_) {
         label.clear(); // Clear label input buffer
      } else {
         if (label.isRequired_) {
            Log::file() << "Error reading label"        << std::endl;
            Log::file() << "Expected: " << label.string_ << std::endl;
            Log::file() << "Scanned:  " << label.input_ << std::endl;
            UTIL_THROW("Incorrect label");
         }
      };
      return in;
   }

   /*
   * Insert a Label into an output stream.
   */
   std::ostream& operator<<(std::ostream& out, Label label)
   {
      out << std::left << std::setw(Label::LabelWidth) 
          << label.string_; 
      out << std::right;
      return out;
   }

} 
