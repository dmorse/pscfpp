#ifndef PSPC_SYSTEM_TPP
#define PSPC_SYSTEM_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2019, The Regents of the University of Minnesota
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
#include <pscf/inter/ChiInteraction.h>
#include <pscf/homogeneous/Clump.h>

#include <pspc/field/BFieldComparison.h>
#include <pspc/field/RFieldComparison.h>

#include <util/param/BracketPolicy.h>
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
      wFieldsBasis_(),
      wFieldsRGrid_(),
      cFieldsBasis_(),
      cFieldsRGrid_(),
      f_(),
      c_(),
      fHelmholtz_(0.0),
      pressure_(0.0),
      hasMixture_(false),
      isAllocated_(false),
      hasWFields_(false),
      hasCFields_(false),
      hasSymmetricFields_(false),
      hasSweep_(false)
      // hasIterator_(true)
   {  
      setClassName("System"); 
      domain_.setFileMaster(fileMaster_);
      interactionPtr_ = new ChiInteraction(); 
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
      iteratorPtr_ = iteratorFactoryPtr_->readObject(in, *this, className, isEnd);
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
         hasSweep_ = true;
         sweepPtr_->setSystem(*this);
      } else {
         hasSweep_ = false;
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
   * Allocate memory for fields.
   */
   template <int D>
   void System<D>::allocate()
   {
      // Preconditions
      UTIL_CHECK(hasMixture_);

      // Allocate wFields and cFields
      int nMonomer = mixture_.nMonomer();
      wFieldsBasis_.allocate(nMonomer);
      wFieldsRGrid_.allocate(nMonomer);

      cFieldsBasis_.allocate(nMonomer);
      cFieldsRGrid_.allocate(nMonomer);
      
      tmpFields_.allocate(nMonomer);
      tmpFieldsRGrid_.allocate(nMonomer);
      tmpFieldsKGrid_.allocate(nMonomer);

      for (int i = 0; i < nMonomer; ++i) {
         wFieldsBasis_[i].allocate(basis().nBasis());
         wFieldsRGrid_[i].allocate(mesh().dimensions());

         cFieldsBasis_[i].allocate(basis().nBasis());
         cFieldsRGrid_[i].allocate(mesh().dimensions());

         tmpFields_[i].allocate(basis().nBasis());
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
         if (command == "COMPUTE") {
            // Read w (chemical potential fields) if not done previously 
            if (!hasWFields_) {
               readEcho(in, filename);
               readWBasis(filename);
            }
            // Solve the modified diffusion equation, without iteration
            compute();
         } else
         if (command == "ITERATE") {
            // Read w (chemical potential) fields if not done previously 
            if (!hasWFields_) {
               readEcho(in, filename);
               readWBasis(filename);
            }
            // Attempt iteration to convergence
            int fail = iterate();
            if (fail) {
               readNext = false;
            }
         } else
         if (command == "SWEEP") {
            // Do a series of iterations.
            sweep();
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
         if (command == "WRITE_W_BASIS") {
            readEcho(in, filename);
            writeWBasis(filename);
         } else 
         if (command == "WRITE_W_RGRID") {
            readEcho(in, filename);
            writeWRGrid(filename);
         } else 
         if (command == "WRITE_C_BASIS") {
            readEcho(in, filename);
            writeCBasis(filename);
         } else
         if (command == "WRITE_C_RGRID") {
            readEcho(in, filename);
            writeCRGrid(filename);
         } else
         if (command == "WRITE_C_BLOCK_RGRID") {
            readEcho(in, filename);
            writeBlockCRGrid(filename);
         } else
         if (command == "WRITE_PROPAGATOR") {
            int polymerId, blockId, directionId, segmentId;
            readEcho(in, filename);
            in >> polymerId;
            in >> blockId;
            in >> directionId;
            in >> segmentId;
            Log::file() << Str("polymer ID   ", 21) << polymerId << "\n"
                        << Str("block ID   ", 21) << blockId << "\n"
                        << Str("direction ID ", 21) << directionId << "\n"
                        << Str("segment ID ", 21) << segmentId << std::endl;
            writePropagatorRGrid(filename, polymerId, blockId, directionId, segmentId);
         } else
         if (command == "WRITE_DATA") {
            readEcho(in, filename);
            writeData(filename);
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
         if (command == "CHECK_RGRID_SYMMETRY") {
            readEcho(in, inFileName);
            bool hasSymmetry;
            hasSymmetry = checkRGridFieldSymmetry(inFileName);
            if (hasSymmetry) {
               Log::file() << "Symmetry of r-grid file matches this space group." 
                           << std::endl;
            } else {
               Log::file() << "Symmetry of r-grid file does not match this space group." 
                           << std::endl;
            }
         } else
         if (command == "RHO_TO_OMEGA") {
            readEcho(in, inFileName);
            readEcho(in, outFileName);
            rhoToOmega(inFileName, outFileName);
         } else
         if (command == "OUTPUT_STARS") {
            readEcho(in, outFileName);
            outputStars(outFileName);
         } else
         if (command == "OUTPUT_WAVES") {
            readEcho(in, outFileName);
            outputWaves(outFileName);
         } else {
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

   /*
   * Compute Helmoltz free energy and pressure
   */
   template <int D>
   void System<D>::computeFreeEnergy()
   {
      UTIL_CHECK(hasWFields_);
      UTIL_CHECK(hasCFields_);
      UTIL_CHECK(hasSymmetricFields_);

      // Initialize to zero
      fHelmholtz_ = 0.0;
 
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
            if (phi > 1E-08) {
               fHelmholtz_ += phi*( mu - 1.0 )/length;
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
            if (phi > 1E-08) {
               fHelmholtz_ += phi*( mu - 1.0 )/size;
            }
         }
      }

      int nm  = mixture_.nMonomer();
      int nBasis = basis().nBasis();

      // Compute Legendre transform subtraction
      // Use expansion in symmetry-adapted orthonormal basis
      double temp = 0.0;
      for (int i = 0; i < nm; ++i) {
         for (int k = 0; k < nBasis; ++k) {
            temp += wFieldsBasis_[i][k] * cFieldsBasis_[i][k];
         }
      }

      // Compute contribution from external fields, if fields exist
      if (iterator().hasExternalField()) {
         for (int i = 0; i < nm; ++i) {
            for (int k = 0; k < nBasis; ++k) {
               temp -= iterator().externalField(i)[k] * cFieldsBasis_[i][k];
            }
         }
      }

      fHelmholtz_ -= temp;

      // Compute excess interaction free energy [ phi^{T}*chi*phi ]
      double chi;
      for (int i = 0; i < nm; ++i) {
         for (int j = i + 1; j < nm; ++j) {
            chi = interaction().chi(i,j);
            for (int k = 0; k < nBasis; ++k) {
               fHelmholtz_+= chi * cFieldsBasis_[i][k] * cFieldsBasis_[j][k];
            }
         }
      }

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

      // If the iterator has a mask, then the volume that should be
      // used to calculate free energy/pressure is the volume available
      // to the polymers, not the total unit cell volume. We thus divide
      // the free energy and pressure by (1 - maskField()[0]), the 
      // volume fraction of the unit cell that is occupied by the 
      // polymers. This properly scales them to the correct value.
      if (iterator().hasMask()) {
         fHelmholtz_ /= 1 - iterator().maskField()[0];
         pressure_ /= 1 - iterator().maskField()[0];
      }

   }

   template <int D>
   void System<D>::outputThermo(std::ostream& out) const
   {
      out << std::endl;
      out << "fHelmholtz    " << Dbl(fHelmholtz(), 18, 11) << std::endl;
      out << "pressure      " << Dbl(pressure(), 18, 11) << std::endl;
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

      // Allocate c_ work array, if necessary
      if (c_.isAllocated()) {
         UTIL_CHECK(c_.capacity() == nm);
      } else {
         c_.allocate(nm);
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
               c_[j] = 0.0;
            }
   
            // Compute clump sizes for all monomer types.
            nb = mixture_.polymer(i).nBlock(); 
            for (k = 0; k < nb; ++k) {
               Block<D>& block = mixture_.polymer(i).block(k);
               j = block.monomerId();
               c_[j] += block.length();
            }
    
            // Count the number of clumps of nonzero size
            nc = 0;
            for (j = 0; j < nm; ++j) {
               if (c_[j] > 1.0E-8) {
                  ++nc;
               }
            }
            homogeneous_.molecule(i).setNClump(nc);
    
            // Set clump properties for this Homogeneous::Molecule
            k = 0; // Clump index
            for (j = 0; j < nm; ++j) {
               if (c_[j] > 1.0E-8) {
                  homogeneous_.molecule(i).clump(k).setMonomerId(j);
                  homogeneous_.molecule(i).clump(k).setSize(c_[j]);
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
   * Read w-field in symmetry adapted basis format.
   */
   template <int D>
   void System<D>::readWBasis(const std::string & filename)
   {
      fieldIo().readFieldsBasis(filename, wFieldsBasis_, domain_.unitCell());
      fieldIo().convertBasisToRGrid(wFieldsBasis_, wFieldsRGrid_);
      hasWFields_ = true;
      hasSymmetricFields_ = true;
      hasCFields_ = false;
   }

   /*
   * Read w-fields in real-space grid (r-grid) format.
   */
   template <int D>
   void System<D>::readWRGrid(const std::string & filename)
   {
      fieldIo().readFieldsRGrid(filename, wFieldsRGrid_, domain_.unitCell());
      //fieldIo().convertRGridToBasis(wFieldsRGrid_, wFieldsBasis_);
      hasWFields_ = true;
      hasSymmetricFields_ = false;
      hasCFields_ = false;
   }

   /*
   * Set new w-field values.
   */
   template <int D>
   void System<D>::setWBasis(DArray< DArray<double> > const & fields)
   {
      // Update system wFields
      int nMonomer = mixture_.nMonomer();
      int nBasis = domain_.basis().nBasis();
      for (int i = 0; i < nMonomer; ++i) {
         DArray<double> const & f = fields[i];
         DArray<double> &       w = wFieldsBasis_[i];
         for (int j = 0; j < nBasis; ++j) {
            w[j] = f[j];
         }
      }

      // Update system wFieldsRGrid
      domain_.fieldIo().convertBasisToRGrid(wFieldsBasis_, wFieldsRGrid_);

      hasWFields_ = true;
      hasSymmetricFields_ = true;
      hasCFields_ = false;
   }

   /*
   * Set new w-field values, using r-grid fields as inputs.
   */
   template <int D>
   void System<D>::setWRGrid(DArray<Field> const & fields)
   {
      // Update system wFieldsRGrid
      int nMonomer = mixture_.nMonomer();
      int meshSize = domain_.mesh().size();
      for (int i = 0; i < nMonomer; ++i) {
         Field const & f = fields[i];
         Field& w = wFieldsRGrid_[i];
         for (int j = 0; j < meshSize; ++j) {
            w[j] = f[j];
         }
      }

      // Update system wFieldsRgrid
      // domain_.fieldIo().convertRGridToBasis(wFieldsRGrid_, wFieldsBasis_);

      hasWFields_ = true;
      hasSymmetricFields_ = false;
      hasCFields_ = false;
   }

   /*
   * Solve MDE for current w-fields, without iteration.
   */
   template <int D>
   void System<D>::compute(bool needStress)
   {
      UTIL_CHECK(hasWFields_);

      // Solve the modified diffusion equation (without iteration)
      mixture_.compute(wFieldsRGrid(), cFieldsRGrid_);
      hasCFields_ = true;

      if (hasSymmetricFields_) {

         // Convert c fields from r-grid to basis
         fieldIo().convertRGridToBasis(cFieldsRGrid_, cFieldsBasis_);

         if (needStress) {
            mixture_.computeStress();
         }
      }

   }

   /*
   * Iteratively solve a SCFT problem for specified parameters.
   */
   template <int D>
   int System<D>::iterate()
   {
      UTIL_CHECK(hasWFields_);
      UTIL_CHECK(hasSymmetricFields_);
      hasCFields_ = false;

      Log::file() << std::endl;
      Log::file() << std::endl;

      // Call iterator
      int error = iterator().solve();
      
      hasCFields_ = true;

      if (!error) {   
         if (!iterator().isFlexible()) {
            mixture().computeStress();
         }
         computeFreeEnergy();
         outputThermo(Log::file());
      }
      return error;
   }

   /*
   * Perform sweep along a line in parameter space.
   */
   template <int D>
   void System<D>::sweep()
   {
      UTIL_CHECK(hasWFields_);
      UTIL_CHECK(hasSymmetricFields_);
      Log::file() << std::endl;
      Log::file() << std::endl;

      // Perform sweep
      sweepPtr_->sweep();
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
   * Write w-fields in symmetry-adapted basis format. 
   */
   template <int D>
   void System<D>::writeWBasis(const std::string & filename) const
   {
      UTIL_CHECK(hasWFields_);
      fieldIo().writeFieldsBasis(filename, wFieldsBasis(), unitCell());
   }

   /*
   * Write w-fields in real space grid file format.
   */
   template <int D>
   void System<D>::writeWRGrid(const std::string & filename) const
   {
      UTIL_CHECK(hasWFields_);
      fieldIo().writeFieldsRGrid(filename, wFieldsRGrid(), unitCell());
   }

   /*
   * Write all concentration fields in symmetry-adapted basis format.
   */
   template <int D>
   void System<D>::writeCBasis(const std::string & filename) const
   {
      UTIL_CHECK(hasCFields_);
      fieldIo().writeFieldsBasis(filename, cFieldsBasis_, unitCell());
   }

   /*
   * Write all concentration fields in real space (r-grid) format.
   */
   template <int D>
   void System<D>::writeCRGrid(const std::string & filename) const
   {
      UTIL_CHECK(hasCFields_);
      fieldIo().writeFieldsRGrid(filename, cFieldsRGrid_, unitCell());
   }

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
   * Write the last time slice of the propagator.
   */
   template <int D>
   void System<D>::writePropagatorRGrid(const std::string & filename, 
                                        int polymerId, int blockId, int directionId, int segmentId) 
   const
   {
      RField<D> propField 
              = mixture_.polymer(polymerId).propagator(blockId, directionId).q(segmentId);
      fieldIo().writeFieldRGrid(filename, propField, unitCell());
   }

   /*
   * Write all data associated with the converged solution.
   */
   template <int D>
   void System<D>::writeData(const std::string & filename)
   {
      std::ofstream file;
      fileMaster().openOutputFile(filename, file);
      writeParam(file);
      outputThermo(file);
      file.close();
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
      fieldIo().readFieldsBasis(inFileName, tmpFields_, tmpUnitCell);
      fieldIo().convertBasisToRGrid(tmpFields_, tmpFieldsRGrid_);
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
      fieldIo().convertRGridToBasis(tmpFieldsRGrid_, tmpFields_);
      fieldIo().writeFieldsBasis(outFileName, tmpFields_, tmpUnitCell);
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
   * Convert fields from real-space grid to symmetry-adapted basis format.
   */
   template <int D>
   bool System<D>::checkRGridFieldSymmetry(const std::string & inFileName) 
   const
   {
      UnitCell<D> tmpUnitCell;
      fieldIo().readFieldsRGrid(inFileName, tmpFieldsRGrid_, tmpUnitCell);
      for (int i = 0; i < mixture_.nMonomer(); ++i) {
         bool symmetric = fieldIo().hasSymmetry(tmpFieldsRGrid_[i]);
         if (!symmetric) {
            return false;
         }
      }
      return true;
   }

   /*
   * Construct guess for omega (w-field) from rho (c-field).
   *
   * Modifies wFields and wFieldsRGrid and outputs wFields.
   */
   template <int D>
   void System<D>::rhoToOmega(const std::string & inFileName, 
                              const std::string & outFileName)
   {
      hasCFields_ = false;
      hasWFields_ = false;

      fieldIo().readFieldsBasis(inFileName, tmpFields_, domain_.unitCell());

      // Compute w fields from c fields
      for (int i = 0; i < basis().nBasis(); ++i) {
         for (int j = 0; j < mixture_.nMonomer(); ++j) {
            wFieldsBasis_[j][i] = 0.0;
            for (int k = 0; k < mixture_.nMonomer(); ++k) {
               wFieldsBasis_[j][i] += interaction().chi(j,k) * tmpFields_[k][i];
            }
         }
      }

      // Convert to r-grid format
      fieldIo().convertBasisToRGrid(wFieldsBasis_, wFieldsRGrid_);
      hasWFields_ = true;
      hasSymmetricFields_ = true;
      hasCFields_ = false;

      // Write w field in basis format
      fieldIo().writeFieldsBasis(outFileName, wFieldsBasis(), unitCell());
   }

   /*
   * Write description of symmetry-adapted stars and basis to file.
   */
   template <int D>
   void System<D>::outputStars(const std::string & outFileName) const
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
   void System<D>::outputWaves(const std::string & outFileName) const
   {
      std::ofstream outFile;
      fileMaster_.openOutputFile(outFileName, outFile);
      fieldIo().writeFieldHeader(outFile, mixture_.nMonomer(), 
                                 unitCell());
      basis().outputWaves(outFile);
   }

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

} // namespace Pspc
} // namespace Pscf
#endif
