#ifndef PSCF_MONOMER_H
#define PSCF_MONOMER_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2019, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <string>
#include <iostream>

namespace Pscf
{

   /**
   * Descriptor for a monomer or particle type.
   *
   * Iostream extractor (>>) and inserter (<<) operators are defined for 
   * a Monomer, allowing the description of a monomer to be read from or
   * written to file like a primitive variable. The text representation 
   * contains a monomer name string and the value of the kuhn (statistical 
   * segment) length, as described \ref pscf_Monomer_page "here".
   *
   * Data for all monomers in a system is normally read from a parameter
   * file into an array-valued parameter named "monomers". 
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
      * Unique integer index for monomer type.
      */
      int id() const;

      /**
      * Statistical segment length (random walk step size).
      */
      double kuhn() const;

      /**
      * Set statistical segment length. 
      */
      void setKuhn(double kuhn);

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

      // Integer identifier
      int  id_;

      // Statistical segment length / kuhn length
      double  kuhn_;

      // Species name string
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
      ar & kuhn_;
      ar & name_;
   }

}
#endif
