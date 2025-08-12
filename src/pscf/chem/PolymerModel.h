#ifndef PSCF_POLYMER_MODEL_H
#define PSCF_POLYMER_MODEL_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/archives/serialize.h>
#include <iostream>


namespace Pscf{

   /**
   * Enumeration and functions to specify a model for polymer chains.
   *
   * PSCF allows the use of either a continuous "thread" or a discrete
   * "bead" model for chain conformations. The PolymerModel::Type
   * enumeration type has two allowed values of PolymerModel::Thread
   * and PolymerModel::Bead that identify these two models.
   *
   * Functions defined in the PolymerModel namespace can be used to set
   * and query a globally accessible pseudo private variable of type
   * PolymerModel::Type that specifies which of these two models is being
   * used throughout a program. This model type variable is initialized
   * to PolymerModel::Thread by default.  In PSCF programs that can
   * support either thread or bead models, a choice of model should be
   * set in the early stages of program execution and not changed again
   * before the program terminates.
   */
   namespace PolymerModel{

      using namespace Util;

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
      *
      * \param model  polymer model type enumeration value
      */
      void setModel(PolymerModel::Type model);

      /**
      * How many times has setModel been called?
      *
      * This returns a counter that is initialized to zero, and is
      * incremented every times setModel is invoked. It can be used
      * confirm that polymer model was not changed within a section
      * of code..
      */
      int nSet();

      /**
      * Get the global polymer model type enumeration value.
      */
      PolymerModel::Type model();

      /**
      * Is the thread model in use ?
      */
      bool isThread();

      /**
      * Is the bead model in use ?
      */
      bool isBead();

      /**
      * Make the polymer model immutable.
      *
      * Locking the polymer model is irreversible - there is no unlock
      * function.
      *
      * Advice: Lock the model shortly after it is set in the execution
      * path  of the main program for PSCF executable programs. Do not
      * lock the model in unit tests, so that a unit test runner can run
      * different tests that test the thread and bead model.
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
