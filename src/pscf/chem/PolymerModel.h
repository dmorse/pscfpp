#ifndef PSCF_POLYMER_MODEL_H
#define PSCF_POLYMER_MODEL_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/archives/serialize.h>
#include <iostream>


namespace Pscf{ 
   using namespace Util;

   /**
   * Enumeration and functions to specify a model for polymer chains.
   *
   * PSCF allows the use of two types of model for a polymer chain,
   * which we refer to as the "Thread" and "Bead" models:
   *
   * The "Thread" model treats each polymer block as a continuous 
   * random walk, and uses numerical approximations that attempt 
   * to closely approximate the behavior of that physical model. 
   * Specifically, it uses a pseudo-spectral scheme for solving the
   * modified diffusion equation and algorithms for integration with
   * respect to contour coordinates that converges rapidly to the 
   * behavior of that model with decreasing contour length step size.
   * The thread model describes the length of each block by a real
   * value denoted by the "length" of a block, which can take on 
   * any value.
   *
   * The "Bead" model is instead based directly on a model of each
   * block of the polymer as a string of an integer number of beads
   * connected by Gaussian springs, and attempts to approximate the
   * behavior of such a discrete chain as closely as possible. In
   * this model, the w fields act on the discrete beads, rather 
   * than acting continuously along the chain. The bead model 
   * specifies the length of a block by an integer value for the
   * number of beads within the block, denoted by "nBead". 
   * 
   * The two allowed PolymerModel::Thread and PolymerModel::Bead of 
   * the PolymerModel::Type each identify one of these two models. 
   *
   * Functions defined in the PolymerModel namespace can be used
   * to set and query the value of a pseudo-private 
   * PolymerModel::Type that specifies which model should be used 
   * throughout a program.  The model type is initialized to
   * PolymerModel::Thread by default.  In PSCF programs that can 
   * support either thread or bead models, a choice of model should 
   * be made made in the early stages of program execution and not
   * not changed again before the program terminates.
   */
   namespace PolymerModel{ 
   
      /**
      * Scoped enumeration of polymer models.
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
      * Make the model immutable hereafter (this is irreversible).
      * 
      * Advice: Lock the model early in execution in production code, but 
      * do not lock the model in unit tests.
      */
      void lock();

      /**
      * Get the global polymer model enumeration value.
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
      * Has the model type been locked (i.e., made immutable) ?
      */
      bool isLocked();

      /**
      * Input stream extractor for a PolymerModel::Type enumeration.
      *
      * \ingroup Pscf_Chem_Module
      *
      * \param in input stream
      * \param model value of PolymerModel to be read from file
      */ 
      std::istream& operator >> (std::istream& in, PolymerModel::Type& model); 
   
      /**
      * Output stream inserter for a PolymerModel::Type enumeration.
      *
      * \ingroup Pscf_Chem_Module
      *
      * \param out  output stream
      * \param model  value of PolymerModel to be written 
      */ 
      std::ostream& operator << (std::ostream& out, PolymerModel::Type const & model); 
   
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
