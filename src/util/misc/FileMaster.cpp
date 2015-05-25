/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/global.h>

#include "FileMaster.h"
#ifdef UTIL_MPI
#include <util/mpi/MpiLoader.h>
#endif

#include <sstream>

namespace Util
{

   /*
   * Default constructor.
   */
   FileMaster::FileMaster()
    : paramFileName_(),
      commandFileName_(),
      inputPrefix_(),
      outputPrefix_(),
      directoryIdPrefix_(),
      rootPrefix_(),
      paramFilePtr_(0),
      commandFilePtr_(0),
      hasDirectoryId_(false),
      isCommonControl_(false)
   {  setClassName("FileMaster"); }

   /*
   * Copy constructor.
   */
   FileMaster::FileMaster(const FileMaster& other)
    : paramFileName_(),
      commandFileName_(),
      inputPrefix_(other.inputPrefix_),
      outputPrefix_(other.outputPrefix_),
      directoryIdPrefix_(other.directoryIdPrefix_),
      rootPrefix_(other.rootPrefix_),
      paramFilePtr_(0),
      commandFilePtr_(0),
      hasDirectoryId_(other.hasDirectoryId_),
      isCommonControl_(other.isCommonControl_)
   {}

   /*
   * Destructor.
   */
   FileMaster::~FileMaster()
   {
      if (paramFilePtr_) {
         paramFilePtr_->close();
         delete paramFilePtr_;
      }
      if (commandFilePtr_) {
         commandFilePtr_->close();
         delete commandFilePtr_;
      }
   }

   /*
   * Set root prefix for all path names.
   */
   void FileMaster::setRootPrefix(const std::string& rootPrefix)
   {  rootPrefix_ = rootPrefix; }

   /*
   * Set a directory Id prefix.
   */
   void FileMaster::setDirectoryId(int rank)
   {
      std::stringstream sMyId;
      sMyId << rank;
      directoryIdPrefix_ = sMyId.str();
      directoryIdPrefix_ += "/";
      hasDirectoryId_ = true;
   }

   /*
   * Set to use a common parameter file even if a directory id is set.
   */
   void FileMaster::setCommonControl()
   {  isCommonControl_ = true; }

   /*
   * Set the input prefix.
   */
   void FileMaster::setInputPrefix(const std::string& inputPrefix)
   {  inputPrefix_ = inputPrefix; }

   /*
   * Set the output prefix.
   */
   void FileMaster::setOutputPrefix(const std::string& outputPrefix)
   {  outputPrefix_ = outputPrefix; }

   /*
   * Set the parameter file name.
   */
   void FileMaster::setParamFileName(const std::string& paramFileName)
   {  paramFileName_ = paramFileName; }

   /*
   * Set the command file name.
   */
   void FileMaster::setCommandFileName(const std::string& commandFileName)
   {  commandFileName_ = commandFileName; }

   /*
   * Read parameters from file.
   */
   void FileMaster::readParameters(std::istream &in)
   {
      if (commandFileName_.empty()) {
         read<std::string>(in, "commandFileName",  commandFileName_);
      }
      readOptional<std::string>(in, "inputPrefix",  inputPrefix_);
      readOptional<std::string>(in, "outputPrefix", outputPrefix_);
   }

   /*
   * Get the default parameter stream.
   */
   void FileMaster::loadParameters(Serializable::IArchive &ar)
   {
      loadParameter<std::string>(ar, "inputPrefix",  inputPrefix_);
      loadParameter<std::string>(ar, "outputPrefix", outputPrefix_);
      #ifdef UTIL_MPI
      MpiLoader<Serializable::IArchive> loader(*this, ar);
      loader.load(isCommonControl_);
      #else
      ar >> isCommonControl_;
      #endif
   }

   /*
   * Save internal state to file.
   */
   void FileMaster::save(Serializable::OArchive &ar)
   {
      ar << inputPrefix_;
      ar << outputPrefix_;
      ar << isCommonControl_;
   }

   /*
   * Get the default parameter stream.
   */
   std::istream& FileMaster::paramFile()
   {
      if (paramFilePtr_) {
         return *paramFilePtr_;
      } else {
         if (paramFileName_.empty()) {
            if (!hasDirectoryId_ || isCommonControl_) {
               return std::cin;
            } else {
               paramFileName_ = "param";
            }
         } 
         paramFilePtr_ = new std::ifstream();
         // Log::file() << "Opening parameter file " 
         //            << paramFileName_  << std::endl;
         openControlFile(paramFileName_, *paramFilePtr_);
         return *paramFilePtr_;
      }
   }

   /*
   * Get the command input stream by reference.
   */
   std::istream& FileMaster::commandFile()
   {
      if (commandFilePtr_) {
         return *commandFilePtr_;
      } else {
         if (commandFileName_.empty()) {
            commandFileName_ = "commands";
         } 
         commandFilePtr_ = new std::ifstream();
         //Log::file() << "Opening command file " 
         //            << paramFileName_  << std::endl;
         openControlFile(commandFileName_, *commandFilePtr_);
         return *commandFilePtr_;
      }
   }

   /*
   * Open an ifstream with fully specified path and mode.
   */
   void
   FileMaster::open(const std::string& name, std::ifstream& in,
                    std::ios_base::openmode mode) const
   {
      in.open(name.c_str(), mode);
      if (in.fail()) {
         std::string message = "Error opening input file. Filename: ";
         message += name;
         UTIL_THROW(message.c_str());
      }
   }

   /*
   * Open an ofstream with fully specified path and mode.
   */
   void
   FileMaster::open(const std::string& name, std::ofstream& out,
                    std::ios_base::openmode mode) const
   {
      out.open(name.c_str(), mode);
      if (out.fail()) {
         std::string message = "Error opening output file. Filename: ";
         message += name;
         UTIL_THROW(message.c_str());
      }
   }

   /*
   * Open an input parameter or command control file.
   */
   void FileMaster::openControlFile(const std::string& name, 
                                   std::ifstream& in) const
   {
      std::string filename(rootPrefix_);
      if (hasDirectoryId_ && !isCommonControl_) {
         filename += directoryIdPrefix_;
      }
      filename += name;
      open(filename, in);
   }

   /*
   * Open an input restart file for reading.
   */
   void 
   FileMaster::openRestartIFile(const std::string& name, std::ifstream& in,
                                std::ios_base::openmode mode) const
   {
      std::string filename(rootPrefix_);
      if (hasDirectoryId_) {
         filename += directoryIdPrefix_;
      }
      filename += name;
      open(filename, in, mode);
   }

   /*
   * Open an output restart file for writing.
   */
   void 
   FileMaster::openRestartOFile(const std::string& name, std::ofstream& out,
                                std::ios_base::openmode mode) const
   {
      // Construct filename = outputPrefix_ + name
      std::string filename(rootPrefix_);
      if (hasDirectoryId_) {
         filename = directoryIdPrefix_;
      }
      filename += name;
      open(filename, out, mode);
   }

   /*
   * Open an input dat file with specified mode.
   */
   void
   FileMaster::openInputFile(const std::string& name, std::ifstream& in,
                             std::ios_base::openmode mode) const
   {
      // Construct filename = inputPrefix_ + name
      std::string filename(rootPrefix_);
      if (hasDirectoryId_) {
         filename += directoryIdPrefix_;
      }
      filename += inputPrefix_;
      filename += name;
      open(filename, in, mode);
   }

   /*
   * Open an output data file in specified mode.
   */
   void
   FileMaster::openOutputFile(const std::string& name, 
                              std::ofstream& out, 
                              std::ios_base::openmode mode) const
   {
      // Construct path
      std::string filename(rootPrefix_);
      if (hasDirectoryId_) {
         filename += directoryIdPrefix_;
      }
      filename += outputPrefix_;
      filename += name;

      open(filename, out, mode);
   }

   /*
   * Will paramFile() return std::cin ?
   */
   bool FileMaster::isCommonControl() const
   {  return ((!hasDirectoryId_) || isCommonControl_); }

}
