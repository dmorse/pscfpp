#ifndef UTIL_AUTOCORRELATION_H
#define UTIL_AUTOCORRELATION_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/global.h>
#include "AutoCorrStage.h"               // base class
#include <util/param/ParamComposite.h>   // base class
#include <util/containers/GArray.h>      // member

namespace Util
{

   /**
   * Auto-correlation function, using hierarchical algorithm.
   *
   * This class represents the primary stage of a linked list of
   * AutoCorrStage objects that implement a hierarchical blocking
   * algorithm for an auto-correlation function.
   *
   * \ingroup Accumulators_Module
   */
   template <typename Data, typename Product>
   class AutoCorrelation : public AutoCorrStage<Data, Product>,
                           public ParamComposite
   {

   public:

      /**
      * Constructor
      */
      AutoCorrelation();

      /**
      * Read parameters from file and initialize.
      *
      * \param in input parameter file
      */
      virtual void readParameters(std::istream& in);

      /**
      * Load internal state from an archive.
      *
      * \param ar input/loading archive
      */
      virtual void load(Serializable::IArchive &ar);

      /**
      * Save internal state to an archive.
      *
      * \param ar output/saving archive
      */
      virtual void save(Serializable::OArchive &ar);

      /**
      * Return maximum delay, in primary samples.
      */
      int maxDelay() const;

      using AutoCorrStage<Data, Product>::setParam;
      using AutoCorrStage<Data, Product>::sample;
      using AutoCorrStage<Data, Product>::clear;
      using AutoCorrStage<Data, Product>::serialize;
      using AutoCorrStage<Data, Product>::output;
      using AutoCorrStage<Data, Product>::bufferCapacity;
      using AutoCorrStage<Data, Product>::nSample;
      using AutoCorrStage<Data, Product>::stageInterval;
      using AutoCorrStage<Data, Product>::autoCorrelation;
      using AutoCorrStage<Data, Product>::corrTime;

   protected:

      using AutoCorrStage<Data, Product>::bufferCapacity_;
      using AutoCorrStage<Data, Product>::maxStageId_;
      using AutoCorrStage<Data, Product>::blockFactor_;
      using AutoCorrStage<Data, Product>::allocate;
      using AutoCorrStage<Data, Product>::serializePrivate;

      /**
      * Register a descendant stage.
      *
      * This should be called only by a root stage.
      *
      * \param ptr pointer to a descendant AutoCorrelation.
      */
      virtual void registerDescendant(AutoCorrStage<Data, Product>* ptr);

   private:

      // Pointers to descendant AutoCorrStage objects
      GArray< AutoCorrStage<Data, Product>* > descendants_;

   };

}
#endif
