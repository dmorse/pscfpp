#ifndef PSPC_SYSTEM_TPP
#define PSPC_SYSTEM_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "System.h"

#include <pspc/sweep/Sweep.h>
#include <pspc/sweep/SweepFactory.h>

#include <pspc/iterator/IteratorFactory.h>

#include <pspc/solvers/Mixture.h>
#include <pspc/solvers/Polymer.h>
#include <pspc/solvers/Solvent.h>

#include <pscf/mesh/MeshIterator.h>
#include <pscf/crystal/shiftToMinimum.h>
#include <pscf/inter/Interaction.h>
#include <pscf/inter/Interaction.h>
#include <pscf/homogeneous/Clump.h>

#include <pspc/field/BFieldComparison.h>
#include <pspc/field/RFieldComparison.h>

#include <util/param/BracketPolicy.h>
#include <util/misc/ioUtil.h>
#include <util/format/Str.h>
#include <util/format/Int.h>
#include <util/format/Dbl.h>

#include <ctime>
#include <iomanip>
#include <sstream>
#include <string>
#include <unistd.h>

namespace Pscf {
namespace Pspc
{

   using namespace Util;

   /*
   * Constructor.
   */
   template <int D>
   System<D>::System()
    : mixture_(),
      domain_(),
      fileMaster_(),
      homogeneous_(),
      interactionPtr_(0),
      iteratorPtr_(0),
      iteratorFactoryPtr_(0),
      sweepPtr_(0),
      sweepFactoryPtr_(0),
      w_(),
      c_(),
      f_(),
      e_(),
      fHelmholtz_(0.0),
      pressure_(0.0),
      hasMixture_(false),
      isAllocated_(false)
      // hasIterator_(true)
   {  
      setClassName("System"); 
      domain_.setFileMaster(fileMaster_);
      w_.setFieldIo(domain_.fieldIo());
      interactionPtr_ = new Interaction(); 
      iteratorFactoryPtr_ = new IteratorFactory<D>(*this); 
      sweepFactoryPtr_ = new SweepFactory<D>(*this);
      BracketPolicy::set(BracketPolicy::Optional);
   }

   /*
   * Destructor.
   */
   template <int D>
   System<D>::~System()
   {
      if (interactionPtr_) {
         delete interactionPtr_;
      }
      if (iteratorPtr_) {
         delete iteratorPtr_;
      }
      if (iteratorFactoryPtr_) {
         delete iteratorFactoryPtr_;
      }
      if (sweepPtr_) {
         delete sweepPtr_;
      }
      if (sweepFactoryPtr_) {
         delete sweepFactoryPtr_;
      }
   }

   /*
   * Process command line options.
   */
   template <int D>
   void System<D>::setOptions(int argc, char **argv)
   {
      bool eflag = false;  // echo
      bool pFlag = false;  // param file 
      bool cFlag = false;  // command file 
      bool iFlag = false;  // input prefix
      bool oFlag = false;  // output prefix
      char* pArg = 0;
      char* cArg = 0;
      char* iArg = 0;
      char* oArg = 0;
   
      // Read program arguments
      int c;
      opterr = 0;
      while ((c = getopt(argc, argv, "er:p:c:i:o:f")) != -1) {
         switch (c) {
         case 'e':
            eflag = true;
            break;
         case 'p': // parameter file
            pFlag = true;
            pArg  = optarg;
            break;
         case 'c': // command file
            cFlag = true;
            cArg  = optarg;
            break;
         case 'i': // input prefix
            iFlag = true;
            iArg  = optarg;
            break;
         case 'o': // output prefix
            oFlag = true;
            oArg  = optarg;
            break;
         case '?':
           Log::file() << "Unknown option -" << optopt << std::endl;
           UTIL_THROW("Invalid command line option");
         }
      }
   
      // Set flag to echo parameters as they are read.
      if (eflag) {
         Util::ParamComponent::setEcho(true);
      }

      // If option -p, set parameter file name
      if (pFlag) {
         fileMaster_.setParamFileName(std::string(pArg));
      }

      // If option -c, set command file name
      if (cFlag) {
         fileMaster_.setCommandFileName(std::string(cArg));
      }

      // If option -i, set path prefix for input files
      if (iFlag) {
         fileMaster_.setInputPrefix(std::string(iArg));
      }

      // If option -o, set path prefix for output files
      if (oFlag) {
         fileMaster_.setOutputPrefix(std::string(oArg));
      }

   }

   /*
   * Read parameters and initialize.
   */
   template <int D>
   void System<D>::readParameters(std::istream& in)
   {
      readParamComposite(in, mixture_);
      hasMixture_ = true;

      int nm = mixture_.nMonomer(); 
      int np = mixture_.nPolymer(); 
      int ns = mixture_.nSolvent(); 

      // Initialize homogeneous object NOTE: THIS OBJECT IS NOT USED AT ALL.
      homogeneous_.setNMolecule(np+ns);
      homogeneous_.setNMonomer(nm);
      initHomogeneous();

      interaction().setNMonomer(mixture_.nMonomer());
      readParamComposite(in, interaction());

      readParamComposite(in, domain_);

      mixture_.setMesh(mesh());
      mixture_.setupUnitCell(unitCell());

      allocate();

      // Initialize iterator through the factory and mediator
      std::string className;
      bool isEnd;
      iteratorPtr_ = iteratorFactoryPtr_->readObject(in, *this, 
                                                     className, isEnd);
      if (!iteratorPtr_) {
         std::string msg = "Unrecognized Iterator subclass name ";
         msg += className;
         UTIL_THROW(msg.c_str());
      }
      iterator().setup();

      // Optionally instantiate a Sweep object
      sweepPtr_ = 
         sweepFactoryPtr_->readObjectOptional(in, *this, className, isEnd);
      if (sweepPtr_) {
         sweepPtr_->setSystem(*this);
      }
   }

   /*
   * Read default parameter file.
   */
   template <int D>
   void System<D>::readParam(std::istream& in)
   {
      readBegin(in, className().c_str());
      readParameters(in);  
      readEnd(in);
   }

   /*
   * Read default parameter file.
   */
   template <int D>
   void System<D>::readParam()
   {  readParam(fileMaster_.paramFile()); }

   /*
   * Write parameter file, omitting any sweep block.
   */
   template <int D>
   void System<D>::writeParamNoSweep(std::ostream& out) const
   {
      out << "System{" << std::endl;
      mixture().writeParam(out);
      interaction().writeParam(out);
      domain().writeParam(out);
      iterator().writeParam(out);
      out << "}" << std::endl;
   }

   /*
   * Allocate memory for fields.
   */
   template <int D>
   void System<D>::allocate()
   {
      // Preconditions
      UTIL_CHECK(hasMixture_);

      int nMonomer = mixture_.nMonomer();

      w_.allocate(nMonomer, basis().nBasis(), mesh().dimensions());
      c_.allocate(nMonomer, basis().nBasis(), mesh().dimensions());

      // Allocate temporary work fields 
      tmpFieldsBasis_.allocate(nMonomer);
      tmpFieldsRGrid_.allocate(nMonomer);
      tmpFieldsKGrid_.allocate(nMonomer);
      for (int i = 0; i < nMonomer; ++i) {
         tmpFieldsBasis_[i].allocate(basis().nBasis());
         tmpFieldsRGrid_[i].allocate(mesh().dimensions());
         tmpFieldsKGrid_[i].allocate(mesh().dimensions());
      }
      isAllocated_ = true;
   }

   /*
   * Read a filename string and echo to log file (used in readCommands).
   */
   template <int D>
   void System<D>::readEcho(std::istream& in, std::string& string) const
   {
      in >> string;
      Log::file() << " " << Str(string, 20) << std::endl;
   }

   
   /*
   * Read and execute commands from a specified command file.
   */
   template <int D>
   void System<D>::readCommands(std::istream &in) 
   {
      UTIL_CHECK(isAllocated_);
      std::string command, filename, inFileName, outFileName;

      bool readNext = true;
      while (readNext) {

         in >> command;
         Log::file() << command <<std::endl;

         if (command == "FINISH") {
            Log::file() << std::endl;
            readNext = false;
         } else
         if (command == "READ_W_BASIS") {
            readEcho(in, filename);
            readWBasis(filename);
         } else
         if (command == "READ_W_RGRID") {
            readEcho(in, filename);
            readWRGrid(filename);
         } else
         if (command == "READ_H_BASIS") {
            readEcho(in, filename);
            if (!h_.isAllocated()) {
               h_.allocate(mixture_.nMonomer(), basis().nBasis(), 
                           mesh().dimensions());
            }
            h_.setFieldIo(domain_.fieldIo());
            h_.readBasis(filename, domain_.unitCell());
         } else
         if (command == "READ_H_RGRID") {
            readEcho(in, filename);
            if (!h_.isAllocated()) {
               h_.allocate(mixture_.nMonomer(), basis().nBasis(), 
                           mesh().dimensions());
            }
            h_.setFieldIo(domain_.fieldIo());
            h_.readRGrid(filename, domain_.unitCell());
         } else
         if (command == "READ_MASK_BASIS") {
            readEcho(in, filename);
            if (!mask_.isAllocated()) {
               mask_.allocate(basis().nBasis(), mesh().dimensions());
            }
            mask_.setFieldIo(domain_.fieldIo());
            mask_.readBasis(filename, domain_.unitCell());
         } else
         if (command == "READ_MASK_RGRID") {
            readEcho(in, filename);
            if (!mask_.isAllocated()) {
               mask_.allocate(basis().nBasis(), mesh().dimensions());
            }
            mask_.setFieldIo(domain_.fieldIo());
            mask_.readBasis(filename, domain_.unitCell());
         } else
         if (command == "COMPUTE") {
            // Solve the modified diffusion equation, without iteration
            compute();
         } else
         if (command == "ITERATE") {
            // Attempt iteration to convergence
            bool isContinuation = false;
            int fail = iterate(isContinuation);
            if (fail) {
               readNext = false;
            }
         } else
         if (command == "SWEEP") {
            // Do a series of iterations.
            sweep();
         } else
         if (command == "WRITE_W_BASIS") {
            readEcho(in, filename);
            UTIL_CHECK(w_.hasData());
            UTIL_CHECK(w_.isSymmetric());
            fieldIo().writeFieldsBasis(filename, w_.basis(), unitCell());
         } else 
         if (command == "WRITE_W_RGRID") {
            readEcho(in, filename);
            UTIL_CHECK(w_.hasData());
            fieldIo().writeFieldsRGrid(filename, w_.rgrid(), unitCell());
         } else 
         if (command == "WRITE_C_BASIS") {
            readEcho(in, filename);
            UTIL_CHECK(hasCFields_);
            fieldIo().writeFieldsBasis(filename, c_.basis(), unitCell());
         } else
         if (command == "WRITE_C_RGRID") {
            readEcho(in, filename);
            UTIL_CHECK(hasCFields_);
            fieldIo().writeFieldsRGrid(filename, c_.rgrid(), unitCell());
         } else
         if (command == "WRITE_C_BLOCK_RGRID") {
            readEcho(in, filename);
            writeBlockCRGrid(filename);
         } else
         if (command == "WRITE_Q_SLICE") {
            int polymerId, blockId, directionId, segmentId;
            readEcho(in, filename);
            in >> polymerId;
            in >> blockId;
            in >> directionId;
            in >> segmentId;
            Log::file() << Str("polymer ID  ", 21) << polymerId << "\n"
                        << Str("block ID  ", 21) << blockId << "\n"
                        << Str("direction ID  ", 21) << directionId << "\n"
                        << Str("segment ID  ", 21) << segmentId << std::endl;
            writeQSlice(filename, polymerId, blockId, directionId, 
                                  segmentId);
         } else
         if (command == "WRITE_Q_TAIL") {
            int polymerId, blockId, directionId;
            readEcho(in, filename);
            in >> polymerId;
            in >> blockId;
            in >> directionId;
            Log::file() << Str("polymer ID  ", 21) << polymerId << "\n"
                        << Str("block ID  ", 21) << blockId << "\n"
                        << Str("direction ID  ", 21) << directionId << "\n";
            writeQTail(filename, polymerId, blockId, directionId);
         } else
         if (command == "WRITE_Q") {
            int polymerId, blockId, directionId;
            readEcho(in, filename);
            in >> polymerId;
            in >> blockId;
            in >> directionId;
            Log::file() << Str("polymer ID  ", 21) << polymerId << "\n"
                        << Str("block ID  ", 21) << blockId << "\n"
                        << Str("direction ID  ", 21) << directionId << "\n";
            writeQ(filename, polymerId, blockId, directionId);
         } else
         if (command == "WRITE_Q_ALL") {
            readEcho(in, filename);
            writeQAll(filename);
         } else
         if (command == "WRITE_PARAM") {
            readEcho(in, filename);
            std::ofstream file;
            fileMaster().openOutputFile(filename, file);
            writeParamNoSweep(file);
            file.close();
         } else
         if (command == "WRITE_THERMO") {
            readEcho(in, filename);
            std::ofstream file;
            fileMaster().openOutputFile(filename, file, 
                                        std::ios_base::app);
            writeThermo(file);
            file.close();
         } else
         if (command == "WRITE_STARS") {
            readEcho(in, filename);
            writeStars(filename);
         } else
         if (command == "WRITE_WAVES") {
            readEcho(in, filename);
            writeWaves(filename);
         } else 
         if (command == "WRITE_MASK_BASIS") {
            readEcho(in, filename);
            UTIL_CHECK(mask_.hasData());
            UTIL_CHECK(mask_.isSymmetric());
            fieldIo().writeFieldBasis(filename, mask_.basis(), unitCell());
         } else 
         if (command == "WRITE_MASK_RGRID") {
            readEcho(in, filename);
            UTIL_CHECK(mask_.hasData());
            fieldIo().writeFieldRGrid(filename, mask_.rgrid(), unitCell());
         } else 
         if (command == "WRITE_H_BASIS") {
            readEcho(in, filename);
            UTIL_CHECK(h_.hasData());
            UTIL_CHECK(h_.isSymmetric());
            fieldIo().writeFieldsBasis(filename, h_.basis(), unitCell());
         } else 
         if (command == "WRITE_H_RGRID") {
            readEcho(in, filename);
            UTIL_CHECK(h_.hasData());
            fieldIo().writeFieldsRGrid(filename, h_.rgrid(), unitCell());
         } else 
         if (command == "BASIS_TO_RGRID") {
            readEcho(in, inFileName);
            readEcho(in, outFileName);
            basisToRGrid(inFileName, outFileName);
         } else 
         if (command == "RGRID_TO_BASIS") {
            readEcho(in, inFileName);
            readEcho(in, outFileName);
            rGridToBasis(inFileName, outFileName);
         } else
         if (command == "KGRID_TO_RGRID") {
            readEcho(in, inFileName);
            readEcho(in, outFileName);
            kGridToRGrid(inFileName, outFileName);
         } else
         if (command == "RGRID_TO_KGRID") {
            readEcho(in, inFileName);
            readEcho(in, outFileName);
            rGridToKGrid(inFileName, outFileName);
         } else
         if (command == "BASIS_TO_KGRID") {
            readEcho(in, inFileName);
            readEcho(in, outFileName);
            basisToKGrid(inFileName, outFileName);
         } else
         if (command == "KGRID_TO_BASIS") {
            readEcho(in, inFileName);
            readEcho(in, outFileName);
            kGridToBasis(inFileName, outFileName);
         } else
         if (command == "CHECK_RGRID_SYMMETRY") {
            readEcho(in, inFileName);
            bool hasSymmetry;
            hasSymmetry = checkRGridFieldSymmetry(inFileName);
            if (hasSymmetry) {
               Log::file() 
                   << "Symmetry of r-grid file matches this space group." 
                   << std::endl;
            } else {
               Log::file() 
                   << "Symmetry of r-grid file does not match this space group" 
                   << std::endl
                   << "to within our error threshold of 1E-8."
                   << std::endl;
            }
         } else
         if (command == "GUESS_W_FROM_C") {
            readEcho(in, inFileName);
            readEcho(in, outFileName);
            guessWfromC(inFileName, outFileName);
         } else
         if (command == "COMPARE_BASIS") {

            // Get two filenames for comparison
            std::string filecompare1, filecompare2;
            readEcho(in, filecompare1);
            readEcho(in, filecompare2);
            
            DArray< DArray<double> > Bfield1, Bfield2;
            fieldIo().readFieldsBasis(filecompare1, Bfield1, domain_.unitCell());
            fieldIo().readFieldsBasis(filecompare2, Bfield2, domain_.unitCell());
            // Note: Bfield1 and Bfield2 will be allocated by readFieldsBasis

            // Compare and output report
            compare(Bfield1, Bfield2);

         } else
         if (command == "COMPARE_RGRID") {
            // Get two filenames for comparison
            std::string filecompare1, filecompare2;
            readEcho(in, filecompare1);
            readEcho(in, filecompare2);
            
            DArray< RField<D> > Rfield1, Rfield2;
            fieldIo().readFieldsRGrid(filecompare1, Rfield1, domain_.unitCell());
            fieldIo().readFieldsRGrid(filecompare2, Rfield2, domain_.unitCell());
            // Note: Rfield1 and Rfield2 will be allocated by readFieldsRGrid

            // Compare and output report
            compare(Rfield1, Rfield2);

         } else
         {
            Log::file() << "Error: Unknown command  " 
                        << command << std::endl;
            readNext = false;
         }
      }
   }

   /*
   * Read and execute commands from the default command file.
   */
   template <int D>
   void System<D>::readCommands()
   {  
      if (fileMaster_.commandFileName().empty()) {
         UTIL_THROW("Empty command file name");
      }
      readCommands(fileMaster_.commandFile()); 
   }

   // Chemical Potential Field Modifier Functions

   /*
   * Read w-field in symmetry adapted basis format.
   */
   template <int D>
   void System<D>::readWBasis(const std::string & filename)
   {
      w_.readBasis(filename, domain_.unitCell());
      hasCFields_ = false;
   }

   /*
   * Read w-fields in real-space grid (r-grid) format.
   */
   template <int D>
   void System<D>::readWRGrid(const std::string & filename)
   {
      w_.readRGrid(filename, domain_.unitCell());
      hasCFields_ = false;
   }

   /*
   * Set new w-field values.
   */
   template <int D>
   void System<D>::setWBasis(DArray< DArray<double> > const & fields)
   {
      w_.setBasis(fields);
      hasCFields_ = false;
   }

   /*
   * Set new w-field values, using r-grid fields as inputs.
   */
   template <int D>
   void System<D>::setWRGrid(DArray<Field> const & fields)
   {
      w_.setRGrid(fields);
      hasCFields_ = false;
   }

   // Unit Cell Modifier / Setter 

   /*
   * Set parameters of the associated unit cell.
   */
   template <int D>
   void System<D>::setUnitCell(UnitCell<D> const & unitCell)
   {
      UTIL_CHECK(domain_.unitCell().lattice() == unitCell.lattice());
      domain_.unitCell() = unitCell;
      mixture_.setupUnitCell(unitCell);
   }

   /*
   * Set parameters of the associated unit cell.
   */
   template <int D>
   void System<D>::setUnitCell(FSArray<double, 6> const & parameters)
   {
      UTIL_CHECK(domain_.unitCell().nParameter() == parameters.size());
      domain_.unitCell().setParameters(parameters);
      mixture_.setupUnitCell(domain_.unitCell());
   }

   // Primary SCFT Computations

   /*
   * Solve MDE for current w-fields, without iteration.
   */
   template <int D>
   void System<D>::compute(bool needStress)
   {
      UTIL_CHECK(w_.hasData());

      // Solve the modified diffusion equation (without iteration)
      mixture_.compute(w_.rgrid(), c_.rgrid(), mask_.phiTot());
      hasCFields_ = true;

      if (w_.isSymmetric()) {
         domain_.fieldIo().convertRGridToBasis(c_.rgrid(), c_.basis());
         if (needStress) {
            mixture_.computeStress();
         }
      }

   }

   /*
   * Iteratively solve a SCFT problem for specified parameters.
   */
   template <int D>
   int System<D>::iterate(bool isContinuation)
   {
      UTIL_CHECK(w_.hasData());
      UTIL_CHECK(w_.isSymmetric());
      hasCFields_ = false;

      Log::file() << std::endl;
      Log::file() << std::endl;

      // Call iterator (return 0 for convergence, 1 for failure)
      int error = iterator().solve(isContinuation);
      hasCFields_ = true;

      // If converged, compute related properties
      if (!error) {   
         if (!iterator().isFlexible()) {
            mixture().computeStress();
         }
         computeFreeEnergy();
         writeThermo(Log::file());
      }
      return error;
   }

   /*
   * Perform sweep along a line in parameter space.
   */
   template <int D>
   void System<D>::sweep()
   {
      UTIL_CHECK(w_.hasData());
      UTIL_CHECK(w_.isSymmetric());
      UTIL_CHECK(hasSweep());
      Log::file() << std::endl;
      Log::file() << std::endl;

      // Perform sweep
      sweepPtr_->sweep();
   }
   
   // Thermodynamic Properties

   /*
   * Compute Helmoltz free energy and pressure
   */
   template <int D>
   void System<D>::computeFreeEnergy()
   {
      UTIL_CHECK(hasCFields_);
      UTIL_CHECK(w_.hasData());
      UTIL_CHECK(w_.isSymmetric());

      // Initialize to zero
      fHelmholtz_ = 0.0;
      fIdeal_ = 0.0;
      fInter_ = 0.0;
 
      double phi, mu;
      int np = mixture_.nPolymer();
      int ns = mixture_.nSolvent();

      // Compute polymer ideal gas contributions to fHelhmoltz_
      if (np > 0) {
         Polymer<D>* polymerPtr;
         double length;
         for (int i = 0; i < np; ++i) {
            polymerPtr = &mixture_.polymer(i);
            phi = polymerPtr->phi();
            mu = polymerPtr->mu();
            length = polymerPtr->length();
            // Recall: mu = ln(phi/q)
            if (phi > 1.0E-08) {
               fIdeal_ += phi*( mu - 1.0 )/length;
            }
         }
      }

      // Compute solvent ideal gas contributions to fHelhmoltz_
      if (ns > 0) {
         Solvent<D>* solventPtr;
         double size;
         for (int i = 0; i < ns; ++i) {
            solventPtr = &mixture_.solvent(i);
            phi = solventPtr->phi();
            mu = solventPtr->mu();
            size = solventPtr->size();
            if (phi > 1.0E-08) {
               fIdeal_ += phi*( mu - 1.0 )/size;
            }
         }
      }

      int nm  = mixture_.nMonomer();
      int nBasis = basis().nBasis();

      double temp(0.0);
      // Compute Legendre transform subtraction
      // Use expansion in symmetry-adapted orthonormal basis
      for (int i = 0; i < nm; ++i) {
         for (int k = 0; k < nBasis; ++k) {
            temp -= w_.basis(i)[k] * c_.basis(i)[k];
         }
      }

      // If the system has a mask, then the volume that should be used
      // in calculating free energy/pressure is the volume available to
      // the polymers, not the total unit cell volume. We thus divide
      // all terms that involve integrating over the unit cell volume by
      // mask().phiTot(), the volume fraction of the unit cell that is 
      // occupied by the polymers. This properly scales them to the 
      // correct value. fExt_, fInter_, and the Legendre transform 
      // component of fIdeal_ all require this scaling. If no mask is 
      // present, mask.phiTot() = 1 and no scaling occurs.
      temp /= mask().phiTot(); 
      fIdeal_ += temp;
      fHelmholtz_ += fIdeal_;

      // Compute contribution from external fields, if fields exist
      if (hasExternalFields()) {
         fExt_ = 0.0;
         for (int i = 0; i < nm; ++i) {
            for (int k = 0; k < nBasis; ++k) {
               fExt_ += h_.basis(i)[k] * c_.basis(i)[k];
            }
         }
         fExt_ /= mask().phiTot();
         fHelmholtz_ += fExt_;
      }

      // Compute excess interaction free energy [ phi^{T}*chi*phi ]
      double chi;
      for (int i = 0; i < nm; ++i) {
         for (int j = i + 1; j < nm; ++j) {
            chi = interaction().chi(i,j);
            for (int k = 0; k < nBasis; ++k) {
               fInter_ += chi * c_.basis(i)[k] * c_.basis(j)[k];
            }
         }
      }
      fInter_ /= mask().phiTot();
      fHelmholtz_ += fInter_;

      // Initialize pressure
      pressure_ = -fHelmholtz_;

      // Polymer corrections to pressure
      if (np > 0) {
         Polymer<D>* polymerPtr;
         double length;
         for (int i = 0; i < np; ++i) {
            polymerPtr = &mixture_.polymer(i);
            phi = polymerPtr->phi();
            mu = polymerPtr->mu();
            length = polymerPtr->length();
            if (phi > 1E-08) {
               pressure_ += mu * phi /length;
            }
         }
      }

      // Solvent corrections to pressure
      if (ns > 0) {
         Solvent<D>* solventPtr;
         double size;
         for (int i = 0; i < ns; ++i) {
            solventPtr = &mixture_.solvent(i);
            phi = solventPtr->phi();
            mu = solventPtr->mu();
            size = solventPtr->size();
            if (phi > 1E-08) {
               pressure_ += mu * phi /size;
            }
         }
      }

   }

   /*
   * Write thermodynamic properties to file.
   */
   template <int D>
   void System<D>::writeThermo(std::ostream& out) const
   {
      out << std::endl;
      out << "fHelmholtz    " << Dbl(fHelmholtz(), 18, 11) << std::endl;
      out << "pressure      " << Dbl(pressure(), 18, 11) << std::endl;
      out << std::endl;
      
      out << "Free energy components:" << std::endl;
      out << "fIdeal        " << Dbl(fIdeal_, 18, 11) << std::endl;
      out << "fInter        " << Dbl(fInter_, 18, 11) << std::endl;
      if (hasExternalFields()) {
         out << "fExt          " << Dbl(fExt_, 18, 11) << std::endl;
      }
      out << std::endl;

      int np = mixture_.nPolymer();
      int ns = mixture_.nSolvent();

      if (np > 0) {
         out << "Polymers:" << std::endl;
         out << "     "
             << "        phi         "
             << "        mu          " 
             << std::endl;
         for (int i = 0; i < np; ++i) {
            out << Int(i, 5) 
                << "  " << Dbl(mixture_.polymer(i).phi(),18, 11)
                << "  " << Dbl(mixture_.polymer(i).mu(), 18, 11)  
                << std::endl;
         }
         out << std::endl;
      }

      if (ns > 0) {
         out << "Solvents:" << std::endl;
         out << "     "
             << "        phi         "
             << "        mu          " 
             << std::endl;
         for (int i = 0; i < ns; ++i) {
            out << Int(i, 5) 
                << "  " << Dbl(mixture_.solvent(i).phi(),18, 11)
                << "  " << Dbl(mixture_.solvent(i).mu(), 18, 11)  
                << std::endl;
         }
         out << std::endl;
      }

   }

   template <int D>
   void System<D>::initHomogeneous()
   {

      // Set number of molecular species and monomers
      int nm = mixture_.nMonomer(); 
      int np = mixture_.nPolymer(); 
      int ns = mixture_.nSolvent(); 
      UTIL_CHECK(homogeneous_.nMolecule() == np + ns);
      UTIL_CHECK(homogeneous_.nMonomer() == nm);

      // Allocate e_ work array, if necessary
      if (e_.isAllocated()) {
         UTIL_CHECK(e_.capacity() == nm);
      } else {
         e_.allocate(nm);
      }

      int i;   // molecule index
      int j;   // monomer index
      int k;   // block or clump index
      int nb;  // number of blocks
      int nc;  // number of clumps
 
      // Loop over polymer molecule species
      if (np > 0) {
         for (i = 0; i < np; ++i) {
   
            // Initial array of clump sizes 
            for (j = 0; j < nm; ++j) {
               e_[j] = 0.0;
            }
   
            // Compute clump sizes for all monomer types.
            nb = mixture_.polymer(i).nBlock(); 
            for (k = 0; k < nb; ++k) {
               Block<D>& block = mixture_.polymer(i).block(k);
               j = block.monomerId();
               e_[j] += block.length();
            }
    
            // Count the number of clumps of nonzero size
            nc = 0;
            for (j = 0; j < nm; ++j) {
               if (e_[j] > 1.0E-8) {
                  ++nc;
               }
            }
            homogeneous_.molecule(i).setNClump(nc);
    
            // Set clump properties for this Homogeneous::Molecule
            k = 0; // Clump index
            for (j = 0; j < nm; ++j) {
               if (e_[j] > 1.0E-8) {
                  homogeneous_.molecule(i).clump(k).setMonomerId(j);
                  homogeneous_.molecule(i).clump(k).setSize(e_[j]);
                  ++k;
               }
            }
            homogeneous_.molecule(i).computeSize();
   
         }
      }

      // Add solvent contributions
      if (ns > 0) {
         double size;
         int monomerId;
         for (int is = 0; is < ns; ++is) {
            i = is + np;
            monomerId = mixture_.solvent(is).monomerId();
            size = mixture_.solvent(is).size();
            homogeneous_.molecule(i).setNClump(1);
            homogeneous_.molecule(i).clump(0).setMonomerId(monomerId);
            homogeneous_.molecule(i).clump(0).setSize(size);
            homogeneous_.molecule(i).computeSize();
         }
      }

   }

   // Command functions

   /*
   * Write all concentration fields in real space (r-grid) format, for 
   * each block (or solvent) individually rather than for each species.
   */
   template <int D>
   void System<D>::writeBlockCRGrid(const std::string & filename) const
   {
      UTIL_CHECK(hasCFields_);

      // Create and allocate the DArray of fields to be written
      DArray<Field> blockCFields;
      blockCFields.allocate(mixture_.nSolvent() + mixture_.nBlock());
      int n = blockCFields.capacity();
      for (int i = 0; i < n; i++) {
         blockCFields[i].allocate(mesh().dimensions());
      }

      // Get data from Mixture and write to file
      mixture_.createBlockCRGrid(blockCFields);
      fieldIo().writeFieldsRGrid(filename, blockCFields, unitCell());
   }

   /*
   * Write the last time slice of the propagator in r-grid format.
   */
   template <int D>
   void System<D>::writeQSlice(const std::string & filename, 
                               int polymerId, int blockId, 
                               int directionId, int segmentId) 
   const
   {
      UTIL_CHECK(polymerId >= 0);
      UTIL_CHECK(polymerId < mixture_.nPolymer());
      Polymer<D> const& polymer = mixture_.polymer(polymerId);
      UTIL_CHECK(blockId >= 0);
      UTIL_CHECK(blockId < polymer.nBlock());
      UTIL_CHECK(directionId >= 0);
      UTIL_CHECK(directionId <= 1);
      Propagator<D> const& propagator 
                               = polymer.propagator(blockId, directionId);
      RField<D> const& field = propagator.q(segmentId);
      fieldIo().writeFieldRGrid(filename, field, unitCell());
   }

   /*
   * Write the last time slice of the propagator in r-grid format.
   */
   template <int D>
   void System<D>::writeQTail(const std::string & filename, 
                              int polymerId, int blockId, int directionId)
   const
   {
      UTIL_CHECK(polymerId >= 0);
      UTIL_CHECK(polymerId < mixture_.nPolymer());
      Polymer<D> const& polymer = mixture_.polymer(polymerId);
      UTIL_CHECK(blockId >= 0);
      UTIL_CHECK(blockId < polymer.nBlock());
      UTIL_CHECK(directionId >= 0);
      UTIL_CHECK(directionId <= 1);
      RField<D> const& 
            field = polymer.propagator(blockId, directionId).tail();
      fieldIo().writeFieldRGrid(filename, field, unitCell());
   }

   /*
   * Write the propagator for a block and direction.
   */
   template <int D>
   void System<D>::writeQ(const std::string & filename, 
                          int polymerId, int blockId, int directionId)
   const
   {
      UTIL_CHECK(polymerId >= 0);
      UTIL_CHECK(polymerId < mixture_.nPolymer());
      Polymer<D> const& polymer = mixture_.polymer(polymerId);
      UTIL_CHECK(blockId >= 0);
      UTIL_CHECK(blockId < polymer.nBlock());
      UTIL_CHECK(directionId >= 0);
      UTIL_CHECK(directionId <= 1);
      Propagator<D> const& propagator 
                              = polymer.propagator(blockId, directionId);
      int ns = propagator.ns();

      // Open file
      std::ofstream file;
      fileMaster_.openOutputFile(filename, file);

      // Write header
      fieldIo().writeFieldHeader(file, 1, unitCell());
      file << "ngrid" << std::endl
           << "          " << mesh().dimensions() << std::endl
           << "nslice"    << std::endl
           << "          " << ns << std::endl;

      // Write data
      bool hasHeader = false;
      for (int i = 0; i < ns; ++i) {
          file << "slice " << i << std::endl;
          fieldIo().writeFieldRGrid(file, propagator.q(i), unitCell(), 
                                    hasHeader);
      }
   }

   /*
   * Write propagators for all blocks of all polymers to files.
   */
   template <int D>
   void System<D>::writeQAll(std::string const & basename)
   {
      std::string filename;
      int np, nb, ip, ib, id;
      np = mixture_.nPolymer();
      for (ip = 0; ip < np; ++ip) {
         //Polymer<D> const * polymerPtr = &mixture_.polymer(ip);
         //nb = polymerPtr->nBlock();
         nb = mixture_.polymer(ip).nBlock();
         for (ib = 0; ib < nb; ++ib) {
            for (id = 0; id < 2; ++id) {
               filename = basename;
               filename += "_";
               filename += toString(ip);
               filename += "_";
               filename += toString(ib);
               filename += "_";
               filename += toString(id);
               filename += ".rf";
               writeQ(filename, ip, ib, id);
            }
         }
      }
   }

   /*
   * Write description of symmetry-adapted stars and basis to file.
   */
   template <int D>
   void System<D>::writeStars(const std::string & outFileName) const
   {
      std::ofstream outFile;
      fileMaster_.openOutputFile(outFileName, outFile);
      fieldIo().writeFieldHeader(outFile, mixture_.nMonomer(),
                                 unitCell());
      basis().outputStars(outFile);
   }

   /*
   * Write a list of waves and associated stars to file.
   */
   template <int D>
   void System<D>::writeWaves(const std::string & outFileName) const
   {
      std::ofstream outFile;
      fileMaster_.openOutputFile(outFileName, outFile);
      fieldIo().writeFieldHeader(outFile, mixture_.nMonomer(), 
                                 unitCell());
      basis().outputWaves(outFile);
   }

   // Field conversion command functions

   /*
   * Convert fields from symmetry-adpated basis to real-space grid format.
   */
   template <int D>
   void System<D>::basisToRGrid(const std::string & inFileName,
                                const std::string & outFileName) const
   {
      UnitCell<D> tmpUnitCell;
      fieldIo().readFieldsBasis(inFileName, tmpFieldsBasis_, tmpUnitCell);
      fieldIo().convertBasisToRGrid(tmpFieldsBasis_, tmpFieldsRGrid_);
      fieldIo().writeFieldsRGrid(outFileName, tmpFieldsRGrid_, 
                                 tmpUnitCell);
   }

   /*
   * Convert fields from real-space grid to symmetry-adapted basis format.
   */
   template <int D>
   void System<D>::rGridToBasis(const std::string & inFileName,
                                const std::string & outFileName) const
   {
      UnitCell<D> tmpUnitCell;
      fieldIo().readFieldsRGrid(inFileName, tmpFieldsRGrid_, tmpUnitCell);
      fieldIo().convertRGridToBasis(tmpFieldsRGrid_, tmpFieldsBasis_);
      fieldIo().writeFieldsBasis(outFileName, tmpFieldsBasis_, tmpUnitCell);
   }

   /*
   * Convert fields from Fourier (k-grid) to real-space (r-grid) format.
   */
   template <int D>
   void System<D>::kGridToRGrid(const std::string & inFileName,
                                const std::string& outFileName) const
   {
      UnitCell<D> tmpUnitCell;
      fieldIo().readFieldsKGrid(inFileName, tmpFieldsKGrid_, tmpUnitCell);
      for (int i = 0; i < mixture_.nMonomer(); ++i) {
         fft().inverseTransform(tmpFieldsKGrid_[i], tmpFieldsRGrid_[i]);
      }
      fieldIo().writeFieldsRGrid(outFileName, tmpFieldsRGrid_, 
                                 tmpUnitCell);
   }

   /*
   * Convert fields from real-space (r-grid) to Fourier (k-grid) format.
   */
   template <int D>
   void System<D>::rGridToKGrid(const std::string & inFileName,
                                const std::string & outFileName) const
   {
      UnitCell<D> tmpUnitCell;
      fieldIo().readFieldsRGrid(inFileName, tmpFieldsRGrid_, 
                                tmpUnitCell);
      for (int i = 0; i < mixture_.nMonomer(); ++i) {
         fft().forwardTransform(tmpFieldsRGrid_[i], tmpFieldsKGrid_[i]);
      }
      fieldIo().writeFieldsKGrid(outFileName, tmpFieldsKGrid_, 
                                 tmpUnitCell);
   }

   /*
   * Convert fields from Fourier (k-grid) to symmetry-adapted basis format.
   */
   template <int D>
   void System<D>::kGridToBasis(const std::string & inFileName,
                                const std::string& outFileName) const
   {
      UnitCell<D> tmpUnitCell;
      fieldIo().readFieldsKGrid(inFileName, tmpFieldsKGrid_, tmpUnitCell);
      fieldIo().convertKGridToBasis(tmpFieldsKGrid_, tmpFieldsBasis_);
      fieldIo().writeFieldsBasis(outFileName, tmpFieldsBasis_, tmpUnitCell);
   }

   /*
   * Convert fields from symmetry-adapted basis to Fourier (k-grid) format.
   */
   template <int D>
   void System<D>::basisToKGrid(const std::string & inFileName,
                                const std::string & outFileName) const
   {
      UnitCell<D> tmpUnitCell;
      fieldIo().readFieldsBasis(inFileName, tmpFieldsBasis_, tmpUnitCell);
      fieldIo().convertBasisToKGrid(tmpFieldsBasis_, tmpFieldsKGrid_);
      fieldIo().writeFieldsKGrid(outFileName, tmpFieldsKGrid_, tmpUnitCell);
   }

   /*
   * Convert fields from real-space grid to symmetry-adapted basis format.
   */
   template <int D>
   bool System<D>::checkRGridFieldSymmetry(const std::string & inFileName) 
   const
   {
      UnitCell<D> tmpUnitCell;
      fieldIo().readFieldsRGrid(inFileName, tmpFieldsRGrid_, tmpUnitCell);
      for (int i = 0; i < mixture_.nMonomer(); ++i) {
         bool symmetric = fieldIo().hasSymmetry(tmpFieldsRGrid_[i],true);
         if (!symmetric) {
            return false;
         }
      }
      return true;
   }

   /*
   * Compare two fields in basis format.
   */ 
   template <int D>
   void System<D>::compare(const DArray< DArray<double> > field1, 
                           const DArray< DArray<double> > field2)
   {
      BFieldComparison comparison(1);
      comparison.compare(field1,field2);

      Log::file() << "\n Basis expansion field comparison results" << std::endl;
      Log::file() << "     Maximum Absolute Difference:   " 
                  << comparison.maxDiff() << std::endl;
      Log::file() << "     Root-Mean-Square Difference:   " 
                  << comparison.rmsDiff() << "\n" << std::endl;
   }

   /*
   * Compare two fields in coordinate grid format.
   */ 
   template <int D>
   void System<D>::compare(const DArray< RField<D> > field1, 
                           const DArray< RField<D> > field2)
   {
      RFieldComparison<D> comparison;
      comparison.compare(field1, field2);

      Log::file() << "\n Real-space field comparison results" << std::endl;
      Log::file() << "     Maximum Absolute Difference:   " 
                  << comparison.maxDiff() << std::endl;
      Log::file() << "     Root-Mean-Square Difference:   " 
                  << comparison.rmsDiff() << "\n" << std::endl;
   }

   /*
   * Construct guess for omega (w-field) from rho (c-field).
   *
   * Modifies wFields and wFieldsRGrid and outputs wFields.
   */
   template <int D>
   void System<D>::guessWfromC(std::string const & inFileName, 
                               std::string const & outFileName)
   {
      int nm = mixture_.nMonomer();
      DArray< DArray<double> > tmpCFieldsBasis;
      tmpCFieldsBasis.allocate(nm);
      for (int i = 0; i < nm; ++i) {
         tmpCFieldsBasis[i].allocate(basis().nBasis());
      }

      fieldIo().readFieldsBasis(inFileName, tmpCFieldsBasis, 
                                domain_.unitCell());

      // Compute w fields from c fields
      for (int i = 0; i < basis().nBasis(); ++i) {
         for (int j = 0; j < nm; ++j) {
            tmpFieldsBasis_[j][i] = 0.0;
            for (int k = 0; k < nm; ++k) {
               tmpFieldsBasis_[j][i] += interaction().chi(j,k) 
                                      * tmpCFieldsBasis[k][i];
            }
         }
      }
      w_.setBasis(tmpFieldsBasis_);

      fieldIo().writeFieldsBasis(outFileName, w_.basis(), unitCell());

      hasCFields_ = false;
   }

} // namespace Pspc
} // namespace Pscf
#endif
