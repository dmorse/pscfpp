#ifndef PSCF_POLYMER_MODEL_H
#define PSCF_POLYMER_MODEL_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/archives/serialize.h>
#include <iostream>


namespace Pscf{ 
   using namespace Util;

   /**
   * Enumeration and functions to specify a model for polymer chains.
   *
   * PSCF allows the use of two types of model for polymer chains, 
   * wich are referred to as the "thread" and "bead" models: The 
   * PolymerModel::Type enumeration type has two allowed values of
   * PolymerModel::Thread and PolymerModel::Bead that each identify 
   * one of these two models. 
   *
   * The "thread" model treats each polymer block as a continuous random
   * walk, and uses numerical approximations that attempt to closely 
   * approximate the behavior of that physical model. Specifically, it 
   * uses a pseudo-spectral scheme with step-size extrapolation to solve 
   * the modified diffusion equation and an Simpson's rule algorithm for 
   * integration with respect to contour coordinates that both converge 
   * rapidly to the behavior of the continuous random walk model with 
   * decreasing contour length step size ds.  The thread model describes 
   * the length of each block by a real value referred to as the "length" 
   * of a block, which can take on any real positive value.
   *
   * The "Bead" model is instead based directly on a model of each block
   * of a block polymer as a string containing an integer number of 
   * discrete beads connected by Gaussian springs.  In the bead model, 
   * the w fields act on the discrete beads, rather * than acting 
   * continuously along the chain. The bead model specifies the length of 
   * a block by an integer value for the number of beads within the block,
   * denoted by "nBead". 
   *
   * Functions defined in the PolymerModel namespace can be used to 
   * set and query the value of a globally accessible pseudo-private 
   * PolymerModel::Type variable that specifies which of these two models 
   * should be used throughout a program.  The model type is initialized 
   * to PolymerModel::Thread by default. In PSCF programs that can support
   * either thread or bead models, a choice of model should be made made 
   * in the early stages of program execution and not changed again 
   * before the program terminates.
   */
   namespace PolymerModel{ 
   
      /**
      * Scoped enumeration of polymer model types.
      *
      * Allowed values:
      *
      *   - Thread: continuous thread model (4th order accurate)
      *
      *   - Bead: discrete bead-spring model
      *
      * \ingroup Pscf_Chem_Module
      */
      enum Type {Thread, Bead};
   
      /**
      * Set the global polymer model enumeration value.
      */
      void setModel(PolymerModel::Type model);
   
      /**
      * Get the global polymer model type enumeration value.
      */
      PolymerModel::Type model();
   
      /**
      * Is the global polymer model a thread model ?
      */
      bool isThread();

      /**
      * Is the global polymer model a discrete bead-spring model ?
      */
      bool isBead();

      /**
      * How many times has the choice of model been set?
      *
      * This returns a counter that is initialized to zero, and is 
      * incremented every times setModel is invoked.
      */
      int nSet();

      /**
      * Make the model immutable hereafter (this is irreversible).
      * 
      * Advice: Lock the model early in execution in the main program for
      * production code, but do not lock the model in unit tests.
      */
      void lock();

      /**
      * Has the model type been locked (i.e., made immutable) ?
      */
      bool isLocked();

      /**
      * Input stream extractor for a PolymerModel::Type enumeration.
      *
      * \ingroup Pscf_Chem_Module
      *
      * \param in  input stream
      * \param model  value of PolymerModel::Type to be read from in
      */ 
      std::istream& operator >> (std::istream& in, 
                                 PolymerModel::Type& model); 
   
      /**
      * Output stream inserter for a PolymerModel::Type enumeration.
      *
      * \ingroup Pscf_Chem_Module
      *
      * \param out  output stream
      * \param model  value of PolymerModel::Type to be written to out
      */ 
      std::ostream& operator << (std::ostream& out, 
                                 PolymerModel::Type const & model); 
   
      /**
      * Serialize a PolymerModel::Type enumeration.
      *
      * \ingroup Pscf_Chem_Module
      *
      * \param ar  archive
      * \param data  enumeration data to be serialized
      * \param version  version id
      */ 
      template <class Archive>
      inline void 
      serialize(Archive& ar, PolymerModel::Type& data, 
                const unsigned int version)
      {  serializeEnum(ar, data, version); }

   }
}
#endif
