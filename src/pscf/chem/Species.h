#ifndef PSCF_SPECIES_H
#define PSCF_SPECIES_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/ParamComposite.h>   // base class
#include <pscf/chem/Ensemble.h>          // member

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
   * problem must use the function setQ(WT) to set set the value of 
   * the molecular partition function q and to compute a new value 
   * for either mu or phi, depending on the ensemble. 
   *
   * \ingroup Pscf_Chem_Module
   */
   template <typename WT = double>
   class Species : public ParamComposite
   {
   public:

      /**
      * Constructor.
      */
      Species();

      /**
      * Destructor.
      */
      virtual ~Species() = default;

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
      WT phi() const;

      /**
      * Get the chemical potential for this species (units kT=1).
      */
      WT mu() const;

      /**
      * Get the molecular partition function for this species.
      */
      WT q() const;

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
      void setQ(WT q);

   private:

      /**
      * Species volume fraction.
      */
      WT phi_;

      /**
      * Species chemical potential.
      */
      WT mu_;

      /**
      * Molecular partition function.
      */
      WT q_;

      /**
      * Input value of phi (closed ensemble) or mu (open).
      * 
      * The input value of this control parameter is required
      * to be a real number.
      */
      double phiMu_;

      /**
      * Statistical ensemble for this species (open or closed).
      */
      Ensemble ensemble_;

   };

   /*
   * Get species volume fraction.
   */
   template <typename WT>
   inline WT Species<WT>::phi() const
   {  return phi_; }

   /*
   * Get species chemical potential.
   */
   template <typename WT>
   inline WT Species<WT>::mu() const
   {  return mu_; }

   /*
   * Get species partition function q.
   */
   template <typename WT>
   inline WT Species<WT>::q() const
   {  return q_; }

   /*
   * Get statistical ensemble for this species (open or closed).
   */
   template <typename WT>
   inline Ensemble Species<WT>::ensemble() const
   {  return ensemble_; }

   // Explicit instantiation declaration
   extern template class Species<double>;

}
#endif
