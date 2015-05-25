#ifndef SERIALIZE_TEST_CLASS_H
#define SERIALIZE_TEST_CLASS_H

struct SerializeTestClass 
{

   int    i;
   double d;

   SerializeTestClass()
     : i(0),
       d(0.0)
   {}

   template <class Archive>
   void serialize(Archive& ar, const unsigned int version)
   {
      ar & i;
      ar & d;
   }

};
#endif
