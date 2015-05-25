#ifndef UTIL_GLOBAL_H
#define UTIL_GLOBAL_H

/**
* \file  util/global.h
* \brief File containing preprocessor macros for error handling.
*/

//-----------------------------------------------------------------------------
// MPI Parallelization

#ifdef UTIL_MPI
#include <mpi.h>
#endif

//-----------------------------------------------------------------------------
// Log file

/**
* Include access to a Log file.
*/
#include <util/misc/Log.h>

//-----------------------------------------------------------------------------
// Errors: Assertions and and Exception macros

#ifndef UTIL_DEBUG

/**
* If defined, disable all C assert(...) statements. 
*/
#define NDEBUG        

#endif

#include "assert.h"
#include "misc/Exception.h"

/**
* Macro for the name of the current function (compiler dependent).
*/
#define UTIL_FUNC __PRETTY_FUNCTION__

/**
* Macro for throwing an Exception, reporting function, file and line number.
*/
#ifdef  UTIL_FUNC
  #ifndef UTIL_MPI
    #define UTIL_THROW(msg) throw Exception(UTIL_FUNC, msg, __FILE__, __LINE__)
  #else
    #define UTIL_THROW(msg) { \
      Exception e(UTIL_FUNC, msg, __FILE__, __LINE__); \
      MpiThrow(e); }
  #endif
#else
  #ifndef UTIL_MPI
    #define UTIL_THROW(msg) throw Exception(msg, __FILE__, __LINE__)
  #else
    #define UTIL_THROW(msg) { \
      Exception e(msg, __FILE__, __LINE__); \
      MpiThrow(e); }
  #endif
#endif

/**
* Assertion macro suitable for serial or parallel production code.
*/
#define UTIL_CHECK(condition) \
  if (!(condition)) { UTIL_THROW("Failed assertion: " #condition); }

/**
* Assertion macro suitable for debugging serial or parallel code.
*/
#ifdef NDEBUG
#define UTIL_ASSERT(condition) {}
#else
#define UTIL_ASSERT(condition) \
  if (!(condition)) { UTIL_THROW("Failed assertion: " #condition); }
#endif

#endif
