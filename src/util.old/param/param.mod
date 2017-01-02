namespace Util{

   /**
   * \defgroup Param_Module Parameter File IO
   * \ingroup Util_NS_Module
   *
   * Classes used to read parameters from a parameter file. Any class 
   * that must read values of member variables from a file should be derived 
   * from ParamComposite, which provides methods for reading and writing a 
   * parameter file, using a programmatically defined file format.
   *
   * ParamComponent is an abstract base class. The classes ParamComposite, 
   * Parameter, Begin, End, and Blank are derived directly from ParamComponent. 
   * Parameter, Begin, End, and Blank are "leaf" notes it a tree structure.
   *
   * Each subclasses of Parameter represents a parameter associated with a
   * different type of C++ object. Such subclasses include class templates
   * ScalarParam, CArrayParam, DArrayParam, FArrayParam, CArray2DParam and
   * MatrixParam.  The template ScalarParam represents any parameter that 
   * is associated with either a primitive C type or a user type for which 
   * their exist overloaded "<<" and ">>" file IO operators. The templates
   * CArrayParam, DArrayParam, FArrayParam, CArray2DParam, and MatrixParam
   * difine parameter file formats for different types of 1D and 2D arrays.
   */

}
