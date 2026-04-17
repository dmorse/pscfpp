#ifndef PSCF_SORT_TEST_H
#define PSCF_SORT_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <pscf/math/Sort.h>

#include <vector>
#include <iostream>

using namespace Pscf;

class SortTest : public UnitTest 
{

private:

   std::vector< Sort::Item<double> > items;
   std::vector< Sort::Bunch > bunches;
   int n;
   
public:

   void setUp()
   {
      n = 10;
      items.reserve(n);
      Sort::Item<double> item;
      for (int i = 0; i < n; ++i) {
         item.value = 0.0;
         item.id = i;
         items.push_back(item); 
      }
      items[0].value = 0.8;
      items[1].value = 0.5;
      items[2].value = 0.1;
      items[3].value = 0.2;
      items[4].value = 0.3;
      items[5].value = 0.0;
      items[6].value = 0.5;
      items[7].value = 0.1;
      items[8].value = 0.3;
      items[9].value = 0.8;
   }

   void tearDown()
   {
      items.clear();
      bunches.clear();
      n = 0;
   }
  
   void testItem() 
   {
      printMethod(TEST_FUNC);
      Sort::Item<double> item;
      item.value = 0.3;
      item.id = 4;
      std::cout << std::endl << item.value << "  " << item.id;
   }

   void testSort()
   {
      printMethod(TEST_FUNC);
      Sort::sort(items);
      
      //std::cout << std::endl <<  0 << "  "
      //          << items[0].value << "  " << items[0].id;
      for (int i=1; i < n; ++i) {
         TEST_ASSERT(!(items[i].value < items[i-1].value));
         //std::cout << std::endl << i << "  "
         //          << items[i].value << "  " << items[i].id;
      }
   }

   void testBunch()
   {
      printMethod(TEST_FUNC);
      Sort::sort(items);
      double epsilon = 1.0E-8;
      Sort::findBunches(items, bunches, epsilon);
      int nb = bunches.size();
      int begin, end, size;
      double oldVal, newVal;
      for (int i=0; i < nb; ++i) {
         begin = bunches[i][0];
         end = bunches[i][1];
         size = end - begin;
         TEST_ASSERT(size > 0);
         //std::cout << std::endl << begin << "  " << end;
         if (size > 1) {
            oldVal = items[begin].value;
            for (int j=begin; j < end; ++j) {
               newVal = items[j].value;
               TEST_ASSERT(std::abs(newVal - oldVal) < epsilon);
            }
         }
      }
   }

};

TEST_BEGIN(SortTest)
//TEST_ADD(SortTest, testItem)
TEST_ADD(SortTest, testSort)
TEST_ADD(SortTest, testBunch)
TEST_END(SortTest)

#endif
