#ifndef PARAM_TEST_CLASSES_H
#define PARAM_TEST_CLASSES_H

#include <util/space/Vector.h>
#include <util/space/IntVector.h>
#include <util/containers/Matrix.h>
#include <util/param/ParamComposite.h>
#include <util/archives/serialize.h>

#include <string>

   // Base class for Factory and Manager. 
   class A : public ParamComposite 
   {
   public:

      A()
      { setClassName("A"); }

   };


   // First subclass of A.
   class B : public A 
   {

   public:

      B()
      {  setClassName("B"); }

      virtual ~B()
      { } 

      virtual void readParameters(std::istream& in) 
      {
         readOptional<double>(in, "x", x_); // optional
         readOptional<int>(in, "opt", opt_); // optional
         read<int>(in, "m", m_);
      }

      virtual void loadParameters(Serializable::IArchive& ar) 
      {
         loadParameter<double>(ar, "x", x_, false); // optional
         loadParameter<int>(ar, "opt", opt_, false); // optional
         loadParameter<int>(ar, "m", m_);
      }

      virtual void save(Serializable::OArchive& ar) 
      {
         // ParamComposite::save(ar);
         Parameter::saveOptional(ar, x_, false);
         Parameter::saveOptional(ar, opt_, false);
         ar << m_;
      }

   private:

      double x_;
      int    opt_;
      int    m_;

   };


   // Second subclass of A.
   class C : public A
   {

   public:

      C()
      {  setClassName("C"); }

      virtual ~C()
      { } 

      virtual void readParameters(std::istream& in) 
      {  read<int>(in, "m", m_); }

      virtual void loadParameters(Serializable::IArchive& in) 
      {  loadParameter<int>(in, "m", m_); }

      virtual void save(Serializable::OArchive& ar) 
      {  ar << m_; }

   private:

      int    m_;

   };

   // Third subclass of A.
   class D : public A
   {

   public:

      D()
      { setClassName("D"); }

      virtual ~D()
      { } 

      virtual void readParameters(std::istream& in) 
      {  read<double>(in, "d", d_); }

      virtual void loadParameters(Serializable::IArchive& in) 
      {  loadParameter<double>(in, "d", d_); }

      virtual void save(Serializable::OArchive& ar) 
      {  ar << d_; }

   private:

      double  d_;

   };


   // Another class
   class E  : public ParamComposite
   {

   public:

      E()
      { setClassName("E"); }

      virtual ~E()
      { } // std::cout << "D destructor" << std::endl; 

      virtual void readParameters(std::istream& in) 
      {  read<double>(in, "e", e_); }

      virtual void loadParameters(Serializable::IArchive& in) 
      {  loadParameter<double>(in, "e", e_); }

      virtual void save(Serializable::OArchive& ar) 
      {  ar << e_; }

      double e()
      {  return e_; }

   private:

      double  e_;

   };

   class AFactory : public Factory<A>
   {

   public:

      ~AFactory()
      { } // std::cout << "AFactory destructor" << std::endl; 

      virtual A* factory(const std::string& classname) const 
      {
         A* ptr = 0;

         ptr = trySubfactories(classname);
         if (ptr) return ptr;

         if (classname == "B") {
            ptr = new B();
         } else
         if (classname == "C") {
            ptr = new C();
         } 
         return ptr;
      }

   };

   class CustomAFactory : public AFactory
   {

   public:

      ~CustomAFactory()
      { } 

      virtual A* factory(const std::string& classname) const
      {
         A* ptr = 0;

         // Try subfactories
         ptr = trySubfactories(classname);
         if (ptr) return ptr;

         if (classname == "D") {
            ptr = new D();
         }
 
         return ptr;
      }

   };

   class AManager : public Manager<A>
   {

   public:

      AManager()
      { setClassName("AManager"); }

      ~AManager()
      {} 

      Factory<A>* newDefaultFactory() const
      { return new AFactory(); }

   };


   class AComposite : public ParamComposite
   {
   
   public:
   
      AComposite()
      {  
         setClassName("AComposite");
         value6_.allocate(4); 
         value9_.allocate(2, 2); 
      }
   
      virtual void readParameters(std::istream& in)
      {
         read<int>(in, "value0", value0_);
         readOptional<int>(in, "optInt", optInt_); // optional
         read<long>(in, "value1", value1_);
         read<double>(in, "value2", value2_);
         read<std::string>(in, "str", str_);
         readCArray<int>(in, "value3", value3_, 3);
         readCArray<double>(in, "value4", value4_, 3);
         readCArray2D<double>(in, "value5", value5_[0], 2, 2, 2);
         //readDArray<double>(in, "value6", value6_, 4);
         readOptionalDArray<double>(in, "value6", value6_, 4); // optional
         readBlank(in);
         read<Vector>(in, "value7", value7_);
         read<IntVector>(in, "value8", value8_);
         readDMatrix<double>(in, "value9", value9_, 2, 2);
         readParamComposite(in, e_);
         readParamComposite(in, manager_);
      }
   
      virtual void loadParameters(Serializable::IArchive& ar)
      {
         loadParameter<int>(ar, "value0", value0_);
         loadParameter<int>(ar, "optInt", optInt_, false); // optional
         loadParameter<long>(ar, "value1", value1_);
         loadParameter<double>(ar, "value2", value2_);
         loadParameter<std::string>(ar, "str", str_);
         loadCArray<int>(ar, "value3", value3_, 3);
         loadCArray<double>(ar, "value4", value4_, 3);
         loadCArray2D<double>(ar, "value5", value5_[0], 2, 2, 2);
         loadDArray<double>(ar, "value6", value6_, 4, false); // optional
         addBlank();
         loadParameter<Vector>(ar, "value7", value7_);
         loadParameter<IntVector>(ar, "value8", value8_);
         loadDMatrix<double>(ar, "value9", value9_, 2, 2);
         loadParamComposite(ar, e_);
         //loadParamComposite(ar, manager_);
      }

      virtual void save(Serializable::OArchive& ar) 
      {
         ar & value0_;
         Parameter::saveOptional(ar, optInt_, false);
         ar & value1_;
         ar & value2_;
         ar & str_;
         ar.pack(value3_, 3);
         ar.pack(value4_, 3);
         ar.pack(value5_[0], 2, 2, 2);
         // ar & value6_;
         Parameter::saveOptional(ar, value6_, true);
         ar & value7_;
         ar & value8_;
         ar & value9_;
         e_.save(ar);
         //manager_.save(ar);
      }

   private:
   
      int     value0_;
      int     optInt_;
      long    value1_;
      double  value2_;
      std::string str_;
      int     value3_[3];
      double  value4_[3];
      double  value5_[2][2];
      DArray<double> value6_;
      Vector    value7_;
      IntVector value8_;
      DMatrix<double> value9_;
  
      E         e_; 
      AManager  manager_;
   
   };

   class BComposite : public ParamComposite
   {
   
   public:
   
      BComposite()
      {  
         setClassName("BComposite");
      }
   
      virtual void readParameters(std::istream& in)
      {
         read<int>(in, "value0", value0_);
         readOptional<int>(in, "optInt", optInt_); // optional
         read<long>(in, "value1", value1_);
         read<double>(in, "value2", value2_);
         read<std::string>(in, "str", str_);
      }
   
      virtual void loadParameters(Serializable::IArchive& ar)
      {
         loadParameter<int>(ar, "value0", value0_);
         loadParameter<int>(ar, "optInt", optInt_, false); // optional
         loadParameter<long>(ar, "value1", value1_);
         loadParameter<double>(ar, "value2", value2_);
         loadParameter<std::string>(ar, "str", str_);
      }

      virtual void save(Serializable::OArchive& ar) 
      {
         ar & value0_;
         Parameter::saveOptional(ar, optInt_, false);
         ar & value1_;
         ar & value2_;
         ar & str_;
      }

   private:
   
      int     value0_;
      int     optInt_;
      long    value1_;
      double  value2_;
      std::string str_;
   
   };

#endif
