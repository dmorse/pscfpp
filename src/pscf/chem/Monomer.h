#ifndef PSCF_MONOMER_H
#define PSCF_MONOMER_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <string>
#include <iostream>

namespace Pscf
{

   /**
   * Descriptor for a monomer type.
   *
   * A Monomer has:
   *
   *  - a unique integer monomer id
   *  - a statistical segment length
   *
   * Iostream extractor (>>) and inserter (<<) operators are defined for 
   * a Monomer, allowing the description of a monomer to be read from or
   * written to file like a primitive variable. The text representation 
   * contains only the value of the kuhn (statistical segment) length,
   * as described \ref pscf_Monomer_page "here".
   *
   * Data for all monomers in a system is normally read from a parameter
   * file into an array-valued parameter named "monomers" in which each
   * element is a Monomer object. The id of each Monomer should be set 
   * to its element index within this array.
   *
   * \ingroup Pscf_Chem_Module
   */
   class Monomer
   {
   public:

      /**
      * Constructor.
      */
      Monomer();

      /**
      * Set the integer index for this monomer type.
      *
      * \param id  new value for index 
      */
      void setId(int id);

      /**
      * Set statistical segment length. 
      * 
      * \param kuhn  value of statistical segment length
      */
      void setKuhn(double kuhn);

      /**
      * Unique integer index for monomer type.
      */
      int id() const;

      /**
      * Statistical segment length (random walk step size).
      */
      double kuhn() const;

      /**
      * Serialize to or from an archive.
      *
      * \param ar Archive object 
      * \param version archive format version index
      */
      template <class Archive>
      void serialize(Archive ar, const unsigned int version);

   private:

      // Integer identifier
      int  id_;

      // Statistical segment length / kuhn length
      double  kuhn_;

   //friends

      friend 
      std::istream& operator >> (std::istream& in, Monomer& monomer);

      friend 
      std::ostream& operator << (std::ostream& out, const Monomer& monomer);

   };

   /**
   * Stream extractor (>>) for a Monomer.
   *
   * The text representation is given by the value of the kuhn data member
   * (i.e., the monomer statistical segment length). The type id is thus 
   * not read from a stream, and so must be set explicitly with setId.
   *
   * \param in  input stream
   * \param monomer  Monomer to be read from stream
   * \return modified input stream
   */
   std::istream& operator >> (std::istream& in, Monomer& monomer);

   /**
   * Stream inserter (<<) for a Monomer.
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
   inline double Monomer::kuhn() const
   {  return kuhn_; }

   /*
   * Set statistical segment length.
   *  
   *  \param kuhn  new value for statistical segement length
   */
   inline void Monomer::setKuhn(double kuhn)
   {  kuhn_ = kuhn; }

   /*
   * Serialize to or from an archive.
   */
   template <class Archive>
   void Monomer::serialize(Archive ar, const unsigned int version)
   {
      ar & id_;
      ar & kuhn_;
   }

}
#endif
