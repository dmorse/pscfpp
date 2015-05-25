namespace Util{

   /**
   * \defgroup Boundary_Module Boundary
   * \ingroup Util_NS_Module
   *
   * \brief   Classes that model periodic boundary conditions.
   *
   * The classes OrthorhombicBoundary, MonoclinicBoundary etc. each define 
   * periodic boundary conditions for a periodic unit cell with a specific 
   * crystal system. These all share a generic interface.
   *
   * The choice of which type of boundary is actually used is defined by the 
   * typedef Boundary, which is defined in the file Boundary.h. This 
   * typedef defines the name Boundary to be a synonym for one of the 
   * concrete boundary classes, e.g., for OrthorhombicBoundary.
   *
   * The rest of the Simpatico source code refers to the boundary class only 
   * via the typedef name Boundary.
   */
 
}
