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
      double step() const;

      /**
      * Set statistical segment length. 
      */
      void setStep(double step);

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
   * Set statistical segment length.
   */
   inline void Monomer::setStep(double step)
   {  step_ = step; }

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
