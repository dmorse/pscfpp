namespace Util
{

   /**
   * \defgroup Serialize_Module Serialization
   * \ingroup Util_NS_Module
   *
   * Serialization of C++ objects to/from file or memory.
   *
   * The code in this module provides a system for serializing sequences 
   * of C++ objects to a file or to random access memory. A serialization
   * of an object stores the full internal state and allows the object to 
   * be reconstructed.  The design is loosely based on that of the Boost 
   * serialization library,
   * http://www.boost.org/doc/libs/1_48_0/libs/serialization/doc/index.html
   * but is much simpler (and less powerful) than the Boost library. 
   *
   * \section Archives
   *
   * An archive stores serialized data, either in a file or in RAM. The 
   * definition of an archive used here is very similar to that used in the 
   * Boost serialization library. An archive class may model either a 
   * saving / output archive, to which data is saved, or a loading / input 
   * archive, from which data is loaded.  By convention, the names of 
   * saving/output archive classes end with the string OArchive and the 
   * names of loading/input archive classes end with the string IArchive. 
   *
   * Different archive classes store serialized objects in different forms. 
   * For example, TextFileOArchive and TextFileIArchive are saving and loading 
   * archive classes, respectively, that are wrappers for ofstream or ifstream 
   * file stream objects in which data is stored in a character representation.
   * BinaryFileOArchive and BinaryFileIArchive are saving/output and 
   * loading / input archives that store data in a binary format. 
   * MemoryOArchive and MemoryIArchive are saving and loading archives that 
   * stored data in binary form in a block of random-access memory. 
   *
   * \section Operators Overloaded IO operators
   *
   * Objects may be saved to a saving archive or loaded from a loading 
   * archive using overloaded operators, using the same syntax as that of the
   * Boost library.  Each saving archive class must define method templates
   * that overload the << (insertion) and & operators. These overloaded
   * operators must be equivalent, and must save an object to the archive.
   * If ar is an instance of a saving archive, such as BinaryFileOArchive, 
   * the expressions
   * \code
   *    ar << data; 
   *    ar &  data;
   * \endcode
   * are thus equivalent, and both save the state of variable data into
   * archive ar.  Each loading archive class must instead define template 
   * methods to overload the >> (extractor) and & operator, which must be 
   * equivalent, and which must load an object from the archive. If ar is 
   * an instance of a loading archive, such as BinaryFileIArchive, then 
   * the expressions
   * \code
   *    ar >> data;
   *    ar &  data;
   * \endcode
   * are equivalent, and both load the state of variable data from archive 
   * ar.
   *
   * \section Serialize Serialize Functions
   *
   * Objects of type T can be saved to or loaded from an instance of a
   * class Archive if and only if the compiler can find a function named 
   * serialize with the signature
   * \code
   *     void serialize(Archive& ar, T& data, unsigned int version)
   * \endcode
   * Here, "version" is an integer index that indicates the version of the
   * archive. This version id is normally given by an integer member of the 
   * archive class. The operator & for a class Archive is normally implemented 
   * by a method template 
   * \code
   *
   *  template <typename T>
   *  void Archive::operator & (T& data);
   *  {  serialize(*this, data, version_); }
   *
   * \endcode
   * that simply calls the appropiate serialize method. Here, version_ is
   * an integer member of the Archive class that stores the archive version 
   * id. Similar templates must be provided for the << or >> operator.
   *
   * Each archive class provides serialize functions for all of the built-in
   * C/C++ types, as well as few other common data types such as std::string. 
   * Definitions of the serialize function for saving archive types must save 
   * (write) data, and those for loading archive types must load (read) data. 
   * 
   * Instances of user-defined classes may also be serialized if an appropriate 
   * serialize function can be found by the compiler. Serialization of instances 
   * of a class T may be enabled by defining either:
   * 
   * - A global serialize function template, with a signature
   * \code
   *
   * template <class Archive>
   * inline void serialize(Archive& ar, T& data, const unsigned int version);
   * 
   * \endcode
   *
   * - A serialize method template in class T, with a signature
   * \code
   *
   * template <class Archive>
   * void T::serialize(Archive& ar, const unsigned int version);
   * 
   * \endcode
   * Note that, in either case, the archive type is normally a template
   * parameter, so that the same serialize function can work with multiple
   * types of archives. 
   *
   * In order to use this system, it is worth understanding how the compiler
   * finds an appropriate serialize method.  When the C++ compiler needs a 
   * serialize method for a particular archive type Archive and data type T, 
   * it will look first for a function serialize(Archive&, T&, unsigned int) 
   * with exactly the required signature, and then for an appropriate template. 
   * Such functions are provided for each archive classes for all of the
   * built-in C/C++ types, and are always used to serialize such types.
   * For class types, their is normally no such non-template function, and 
   * so the compiler will look for an appropriate template, giving priority
   * to templates in which fewer of the function parameters have types given
   * by template arguments, rather than explicit types. If the compiler has
   * access to a global serialize function template for class T with the 
   * signature described above, in which the archive type is a template 
   * parameter but the data type T is explicit, it will use this. If no such 
   * global serialize function template is found, the compiler will try to 
   * compile the following generic template, 
   * \code
   *
   * template <class Archive, typename T>
   * inline void serialize(Archive& ar, T& data, const unsigned int version)
   * {  data.serialize(ar, version); }
   *
   * \endcode
   * which is defined in the file src/util/serialize.h. This template
   * simply calls the serialize method of class T, and so will not 
   * compile if no such method exists. The compiler can thus use, in
   * decreasing order of priority: 1) An explicit serialize function
   * for type T and a specific archive type, 2) A serialize function
   * template for a specific type T in which the archive type is a
   * template parameter, or 3) A serialize method of class T in which
   * the archive type is a template parameter. If none of these are
   * accessible for class T, compilation will fail for any code that
   * attempts to serialize an instance of class T.
   *
   * The use of a single operator & to represent both output (when applied 
   * to a saving archive) and input (when applied to a loading archive), 
   * makes it possible to write a single serialize function template for 
   * each class that specifies how to order save or load instances of 
   * that class, by specifying the order in which members of the class
   * are serialized. For example, consider the following definition of 
   * a simple complex number class:
   * \code 
   *
   *   class Complex  {
   *   public:
   *
   *      A(double real, double imag) : real_(real), imag_(imag) {}
   * 
   *      template <class Archive>
   *      void serialize(Archive& ar, unsigned int version)
   *      { 
   *         ar & real_;
   *         ar & imag_:
   *      }
   *
   *   private:
   *
   *      double real_;
   *      double imag_;   
   * 
   *   } 
   *
   * \endcode
   * The serialize method template provides instructions for the order in 
   * which to either save the two floating point members of the class to
   * a saving archive, or to load them from a loading archive. The use of 
   * a template in which the archive type is a parameter allows a single 
   * serialize method to be used with any type of saving or loading archive.
   *
   * The most serious disadvantage of this system is that, if the serialize
   * method is defined by a template, it cannot also be a virtual method.
   * As a result, the serialize method template for a class cannot be 
   * accessed polymorphically, via a pointer or reference to a base class. 
   * This limitation becomes a problem in designs in which some objects are 
   * accessed only via base class pointers.  The Serializable abstract base
   * class, discussed below, partially solves this problem, by replacing 
   * the serialize method template by a pair of virtual save() and load() 
   * methods.
   *
   * \section Serializable Serializable Classes
   *
   * Serializable is an abstract base class that provides an alternate 
   * interface for serializing objects, using virtual functions rather than 
   * method templates.  Each subclass of Serializable must define virtual 
   * save() and load() methods with the following signatures:
   * \code
   * virtual void save(Serializable::OArchive& ar);
   * virtual void load(Serializable::IArchive& ar);
   * \endcode
   * The typenames Serializable::OArchive and Serializable::IArchive are
   * typedefs that define a pair of archive classes to be used for 
   * serialization. 
   *
   * The advantage of using virtual functions is that it allows these methods 
   * to be accessed polymorphically, via base class pointers or references.
   * The disadvantage is that it requires the hard-coding of a single type
   * type of saving and loading archive. To retain some flexibility, these
   * saving and loading types are defined in the Serializable class by a 
   * pair of typedefs. This allows the type of archives used with Serializable
   * objects to be changed throughout the code by changing these two typedefs
   * and recompiling. 
   *
   * In practice, a serialize method or function template should be defined 
   * for relatively simple, non-polymorphic classes, but polymorhpic classes
   * that are normally accessed via base class pointers need to be derived 
   * from Serializable, and must implement save and load methods.
   */

}
