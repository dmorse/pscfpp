#ifndef CHEM_MONOMER_H
#define CHEM_MONOMER_H

/*
* PFTS - Polymer Field Theory Simulator
*
* Copyright 2013, David Morse (morse012@.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <string>
#include <iostream>

namespace Chem{

   /**
   * Descriptor for a monomer or particle type.
   */
   class Monomer
   {
   public:

      /**
      * Constructor.
      */
      Monomer();

      /**
      * Unique integer index for monomer type.
      */
      int id() const;

      /**
      * Statistical segment length.
      */
      double step() const;

      /**
      * Monomer name string.
      */
      std::string name() const;

      /**
      * Serialize to or from an archive.
      *
      * \param ar Archive object 
      * \param version archive format version index
      */
      template <class Archive>
      void serialize(Archive ar, const unsigned int version);

   private:

      int  id_;
      double  step_;
      std::string  name_;

   //friends

      friend 
      std::istream& operator >> (std::istream& in, Monomer& monomer);

      friend 
      std::ostream& operator << (std::ostream& out, const Monomer& monomer);

   };

   /**
   * istream extractor for a Monomer.
   *
   * \param in  input stream
   * \param monomer  Monomer to be read from stream
   * \return modified input stream
   */
   std::istream& operator >> (std::istream& in, Monomer& monomer);

   /**
   * ostream inserter for a Monomer.
   *
   * \param out  output stream
   * \param monomer  Monomer to be written to stream
   * \return modified output stream
   */
   std::ostream& operator << (std::ostream& out, const Monomer& monomer);

   // inline member functions

   /*
   * Get monomer type index.
   */
   inline int Monomer::id() const
   {  return id_; }

   /*
   * Statistical segment length.
   */
   inline double Monomer::step() const
   {  return step_; }

   /*
   * Monomer name string.
   */
   inline std::string Monomer::name() const
   {  return name_; }

   /*
   * Serialize to or from an archive.
   */
   template <class Archive>
   void Monomer::serialize(Archive ar, const unsigned int version)
   {
      ar & id_;
      ar & step_;
      ar & name_;
   }

}
#endif
