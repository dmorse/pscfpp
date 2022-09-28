#ifndef PSPG_R_FIELD_COMPARISON_TPP
#define PSPG_R_FIELD_COMPARISON_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "RFieldComparison.h"

namespace Pscf {
namespace Pspg {

   // Default Constructor
   template <int D>
   RFieldComparison<D>::RFieldComparison()
   {};

   // Comparator for individual fields
   template <int D>
   double RFieldComparison<D>::compare(RDField<D> const& a, RDField<D> const& b)
   {
      // Copy data to host.
      int nPoints = a.capacity();
      cudaReal* temp_a = new cudaReal[nPoints];
      cudaReal* temp_b = new cudaReal[nPoints];
      cudaMemcpy(temp_a, a.cDField(), nPoints*sizeof(cudaReal), 
                                      cudaMemcpyDeviceToHost);
      cudaMemcpy(temp_b, b.cDField(), nPoints*sizeof(cudaReal), 
                                      cudaMemcpyDeviceToHost);

      DArray< cudaReal > h_a, h_b;
      h_a.allocate(nPoints);
      h_b.allocate(nPoints);
      for (int j = 0; j < nPoints; j++) {
         h_a[j] = temp_a[j];
         h_b[j] = temp_b[j];
      }
      fieldComparison_.compare(h_a,h_b);
      compared_ = true;

      return fieldComparison_.maxDiff();
   }

   // Comparator for arrays of fields
   template <int D>
   double RFieldComparison<D>::compare(DArray<RDField<D>> const& a, 
                                       DArray<RDField<D>> const& b)
   {
      int nFields = a.capacity();
      int nPoints = a[0].capacity();

      // Array on host (used for cudaMemcpy)
      DArray< cudaReal* > temp_a, temp_b;
      temp_a.allocate(nFields);
      temp_b.allocate(nFields);

      // DArrays on host (for compatibility with FieldComparison)
      DArray< DArray< cudaReal > > h_a, h_b;
      h_a.allocate(nFields);
      h_b.allocate(nFields);
      
      for (int i = 0; i < nFields; i++) {
         temp_a[i] = new cudaReal[nPoints];
         temp_b[i] = new cudaReal[nPoints];
         cudaMemcpy(temp_a[i], a[i].cDField(), nPoints*sizeof(cudaReal), 
                                               cudaMemcpyDeviceToHost);
         cudaMemcpy(temp_b[i], b[i].cDField(), nPoints*sizeof(cudaReal), 
                                               cudaMemcpyDeviceToHost);

         h_a[i].allocate(nPoints);
         h_b[i].allocate(nPoints);

         for (int j = 0; j < nPoints; j++) {
            h_a[i][j] = temp_a[i][j];
            h_b[i][j] = temp_b[i][j];
         }
      }

      // Perform comparison
      fieldComparison_.compare(h_a,h_b);
      compared_ = true;

      return fieldComparison_.maxDiff();
   }

}
}
#endif
