#ifndef PSCF_SPECIES_H
#define PSCF_SPECIES_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/ParamComposite.h>   // base class
#include <util/global.h>

namespace Pscf {

   using namespace Util;

   /**
   * Base class for a molecular species (polymer or solvent).
   *
   * Species is a base class for both polymeric and solvent species 
   * descriptors and solvers.  Each Species has values of phi, mu and q, 
   * and an ensemble enumeration value.  The value of either phi or mu 
   * must be provided as a parameter in the parameter file, thereby also
   * determining the choice of ensemble.  A subclass of Species that 
   * actually solves the single single-molecule statistical mechanics 
   * problem must use the function setQ(double) to set set the value 
   * of the molecular partition function q and to compute a new value 
   * for either mu or phi, depending on the ensemble. 
   *
   * \ingroup Pscf_Chem_Module
   */
   class Species : public ParamComposite
   {
   public:

      /**
      * Statistical ensemble for number of molecules.
      */
      enum Ensemble {Unknown, Closed, Open};

      /**
      * Default constructor.
      */
      Species();

      /**
      * Read phi or mu (but not both) and set ensemble accordingly.
      * 
      * This function either reads a value for either phi (species volume 
      * fraction) and sets the ensemble to Closed or reads mu (species 
      * chemical potential) and sets the ensemble to Open. An Exception is
      * thrown if neither a phi or mu parameter appears in the input stream.
      * 
      * \param in  input stream (parameter file)
      */
      virtual void readParameters(std::istream& in);

      /**
      * Get the overall volume fraction for this species.
      */
      double phi() const;

      /**
      * Get the chemical potential for this species (units kT=1).
      */
      double mu() const;

      /**
      * Get the molecular partition function for this species.
      */
      double q() const;

      /**
      * Get the statistical ensemble for this species (open or closed).
      */
      Ensemble ensemble() const;

      /**
      * Set value of phi (volume fraction), if ensemble is closed.
      *
      * An initial value for phi or mu is normally read from a parameter
      * file. This function is provided for use by a sweep or other
      * procedure in which phi for a species with a closed enesmble is
      * modified after initialization. It is an error to call setPhi for
      * a polymer species with an open ensemble.
      *
      * \throw Exception if ensemble is open
      * \param phi  new volume fraction value for this species
      */
      void setPhi(double phi);

      /**
      * Set value of mu (chemical potential), if ensemble is closed.
      *
      * An initial value for phi or mu is normally read from a parameter
      * file. This function is provided for use in a sweep or other
      * procedure in which mu for a species with an open enesmble is
      * modified after initialization. It is an error to call setMu for
      * a polymer species with a closed ensemble.
      *
      * \throw Exception if ensemble is closed
      * \param mu  new chemical potential value for this species
      */
      void setMu(double mu);

   protected:

      /**
      * Set q and compute phi or mu (depending on the ensemble).
      *
      * This function alllows a subclass to set the value of the 
      * molecular partition function q, and use this value to compute a
      * corresponding value for either mu (for a closed ensemble) or phi 
      * (for an open ensemble). Upon return these variables are related
      * by the equation: \f$ \phi = \exp(\mu) q \f$. 
      *
      * \param q  new value of molecular partition function q
      */ 
      void setQ(double q);

   private:

      /**
      * Species volume fraction.
      */
      double phi_;

      /**
      * Species chemical potential.
      */
      double mu_;

      /**
      * Molecular partition function.
      */
      double q_;

      /**
      * Statistical ensemble for this species (open or closed).
      */
      Ensemble ensemble_;

   };

   /*
   * Get species volume fraction.
   */
   inline double Species::phi() const
   {  return phi_; }

   /*
   * Get species chemical potential.
   */
   inline double Species::mu() const
   {  return mu_; }

   /*
   * Get species partition function q.
   */
   inline double Species::q() const
   {  return q_; }

   /*
   * Get statistical ensemble for this species (open or closed).
   */
   inline Species::Ensemble Species::ensemble() const
   {  return ensemble_; }

   /**
   * istream extractor for a Species::Ensemble enumeration.
   *
   * \param  in       input stream
   * \param  policy   Species::Ensemble to be read
   * \return modified input stream
   */
   std::istream& operator >> (std::istream& in, Species::Ensemble& policy);

   /**
   * ostream inserter for an Species::Ensemble enumeration.
   *
   * Text representations of allowed values are "Open" and "Closed".
   *
   * \param  out      output stream
   * \param  policy   Species::Ensemble to be written
   * \return modified output stream
   */
   std::ostream& operator << (std::ostream& out, Species::Ensemble policy);

   /**
   * Serialize a Species::Ensemble
   *
   * \param ar      archive object
   * \param policy  object to be serialized
   * \param version archive version id
   */
   template <class Archive>
   void serialize(Archive& ar, Species::Ensemble& policy,
                  const unsigned int version)
   {  serializeEnum(ar, policy, version); }

}
#endif
