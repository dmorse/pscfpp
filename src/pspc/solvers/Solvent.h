#ifndef PSPC_SOLVENT_H
#define PSPC_SOLVENT_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2019, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pspc/solvers/Propagator.h>
#include <pscf/mesh/Mesh.h>
#include <pscf/chem/Species.h>             // base class
#include <util/param/ParamComposite.h>     // base class

#include <iostream>

namespace Pscf { 
namespace Pspc { 

   using namespace Util;

   // Forward declaration of class template
   template <int D> class Solvent;

   /**
   * istream extractor for a Solvent.
   *
   * \param in  input stream
   * \param solvent  Solvent to be read from stream
   * \return modified input stream
   */
   template <int D>
   std::istream& operator>> (std::istream& in, Solvent<D>& solvent);

   /**
   * ostream inserter for a Solvent.
   *
   * \param out  output stream
   * \param solvent  Solvent to be written to stream
   * \return modified output stream
   */
   template <int D>
   std::ostream& operator<< (std::ostream& out, const Solvent<D> & solvent);

   /**
   * Solver and descriptor for a solvent species.
   *
   * \ingroup Pspc_Solver_Module
   */
   template <int D>
   class Solvent : public Species, public ParamComposite
   {

   public:

      /**
      * Monomer concentration field type.
      */
      typedef typename Propagator<D>::CField CField;

      /** 
      * Monomer chemical potential field type.
      */
      typedef typename Propagator<D>::WField WField;

      /**
      * Constructor.
      */
      Solvent();
   
      /**
      * Constructor.
      */
      ~Solvent();
   
      /**
      * Read and initialize.
      *
      * \param in input parameter stream
      */
      virtual void readParameters(std::istream& in);

      /**
      * Set association with Mesh and allocate concentration field array.
      *
      * \param mesh associated Mesh<D> object
      */
      void setMesh(Mesh<D> const & mesh);

      /**
      * Compute monomer concentration field, q and phi and/or mu.
      *
      * Upon return, cField, phi, mu, and q are all set.
      *
      * \param wField monomer chemical potential field of relevant type.
      */
      void compute(WField const & wField );
  
      /**
      * Serialization function.
      */ 
      template <class Archive>
      void serialize(Archive& ar, unsigned int);

      /// \name Setters (set member data)
      //@{

      /**
      * Set value of phi (volume fraction), if ensemble is closed.
      *
      * \throw Exception if ensemble is open
      * \param phi desired volume fraction for this species
      */
      void setPhi(double phi);

      /**
      * Set value of mu (chemical potential), if ensemble is closed.
      *
      * \throw Exception if ensemble is open
      * \param phi desired chemical potential for this species
      */
      void setMu(double mu);

      /**
      * Set the monomer id for this solvent.
      *
      * \param monomerId integer id of monomer type (>=0)
      */ 
      void setMonomerId(int monomerId);
  
      /**
      * Set the size or volume of this solvent species.
      *
      * The ``size" is steric volume / reference volume.
      *
      * \param size volume of solvent
      */ 
      void setSize(double size);

      //@}
      /// \name Accessors (getters)
      //@{
 
      /**
      * Get the monomer type id.
      */ 
      int monomerId() const;
  
      /**
      * Get the size (number of monomers) in this solvent.
      */
      double size() const;

      /**
      * Get monomer concentration field for this solvent.
      */
      const CField& cField() const
      {  return cField_;  }
  
      //@}

      // Inherited accessor functions 
      using Pscf::Species::phi;
      using Pscf::Species::mu;
      using Pscf::Species::q;
      using Pscf::Species::ensemble;
      using Pscf::ParamComposite::read;

   protected:

      // Inherited data members
      using Pscf::Species::phi_;
      using Pscf::Species::mu_;
      using Pscf::Species::q_;
      using Pscf::Species::ensemble_;

   private:

      // Concentration field for this solvent
      CField cField_;
  
      /// Identifier for the associated monomer type.
      int monomerId_;

      /// Size of this block = volume / monomer reference volume. 
      double size_;

      // Pointer to associated mesh
      Mesh<D> const *  meshPtr_;

   // friend:
 
      friend 
      std::istream& operator>> (std::istream& in, Solvent& solvent);

      friend 
      std::ostream& operator<< (std::ostream& out, const Solvent& solvent);

   };
   
   // Inline member functions

   /*
   * Get the monomer type id.
   */ 
   template <int D>
   inline int Solvent<D>::monomerId() const
   {  return monomerId_; }

   /*
   * Get the size (number of monomers) in this block.
   */
   template <int D>
   inline double Solvent<D>::size() const
   {  return size_; }
    
   // Non-inline functions

   #if 0
   /*
   * Extract a Solvent<D> from an istream.
   */
   template <int D>
   std::istream& operator>> (std::istream& in, Solvent<D> & solvent)
   {
      in >> solvent.monomerId_;
      in >> solvent.size_;
      return in;
   }

   /*
   * Output a Solvent<D> to an ostream, without line breaks.
   */
   template <int D>
   std::ostream& operator << (std::ostream& out, const Solvent<D> & solvent)
   {
      out << solvent.monomerId_;
      out << "  ";
      out.setf(std::ios::scientific);
      out.width(16);
      out.precision(8);
      out << solvent.size_;
      return out;
   }
   #endif

   /*
   * Serialize to/from an archive.
   */
   template <int D>
   template <class Archive>
   void Solvent<D>::serialize(Archive& ar, unsigned int)
   {
      ar & phi_;
      ar & mu_;
      ar & q_;
      ar & ensemble_;
      ar & monomerId_;
      ar & size_;
   }

   #if 0
   #ifndef PSPC_SOLVENT_TPP
   // Supress implicit instantiation
   //extern template class Solvent<1>;
   //extern template class Solvent<2>;
   //extern template class Solvent<3>;
   //extern template std::istream& operator>> (std::istream&, Solvent<1> &);
   //extern template std::ostream& operator<< (std::ostream&, const Solvent<1>&);
   //extern template std::istream& operator>> (std::istream&, Solvent<2>&);
   //extern template std::ostream& operator<< (std::ostream&, const Solvent<2>&);
   //extern template std::istream& operator>> (std::istream&, Solvent<3>& );
   //extern template std::ostream& operator<< (std::ostream&, const Solvent<3>& );
   #endif
   #endif

}
}
#include "Solvent.tpp" 
#endif 
