#ifndef PSSP_BLOCK_TPP
#define PSSP_BLOCK_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Block.h"
#include <pscf/mesh/Mesh.h>
#include <pscf/mesh/MeshIterator.h>
#include <pscf/crystal/UnitCell.h>
#include <pscf/crystal/shiftToMinimum.h>
#include <util/containers/FMatrix.h>      // member template
#include <util/containers/DArray.h>      // member template
#include <util/containers/FArray.h>      // member template

namespace Pscf { 
namespace Pssp {

   using namespace Util;

   /*
   * Constructor.
   */
   template <int D>
   Block<D>::Block()
    : meshPtr_(0),
      kMeshDimensions_(0),
      ds_(0.0),
      ns_(0)
   {
      propagator(0).setBlock(*this);
      propagator(1).setBlock(*this);
   }

   /*
   * Destructor.
   */
   template <int D>
   Block<D>::~Block()
   {}

   template <int D>
   void Block<D>::setDiscretization(double ds, const Mesh<D>& mesh)
   {  
      UTIL_CHECK(mesh.size() > 1);
      UTIL_CHECK(ds > 0.0);

      // Set association to mesh
      meshPtr_ = &mesh;

      // Set contour length discretization
      ns_ = floor(length()/ds + 0.5) + 1;
      if (ns_%2 == 0) {
         ns_ += 1;
      }
      ds_ = length()/double(ns_ - 1);

      // Compute Fourier space kMeshDimensions_ 
      for (int i = 0; i < D; ++i) {
         if (i < D - 1) {
            kMeshDimensions_[i] = mesh.dimensions()[i];
         } else {
            kMeshDimensions_[i] = mesh.dimensions()[i]/2 + 1;
       }
      }

      // Allocate work arrays
      expKsq_.allocate(kMeshDimensions_);
      expW_.allocate(mesh.dimensions());
      expKsq2_.allocate(kMeshDimensions_);
      expW2_.allocate(mesh.dimensions());
      qr_.allocate(mesh.dimensions());
      qk_.allocate(mesh.dimensions());
      qr2_.allocate(mesh.dimensions());
      qk2_.allocate(mesh.dimensions());
      qf_.allocate(mesh.dimensions());

      q1.allocate(mesh.dimensions());
      q2.allocate(mesh.dimensions());
      q1p.allocate(mesh.dimensions());
      q2p.allocate(mesh.dimensions());

      propagator(0).allocate(ns_, mesh);
      propagator(1).allocate(ns_, mesh);
      cField().allocate(mesh.dimensions());

   }

   /*
   * Setup data that depend on the unit cell parameters.
   */
   template <int D>
   void 
   Block<D>::setupUnitCell(const UnitCell<D>& unitCell)
   {
      MeshIterator<D> iter;
      // std::cout << "kDimensions = " << kMeshDimensions_ << std::endl;
      iter.setDimensions(kMeshDimensions_);
      IntVec<D> G, Gmin;
      double Gsq;
      double factor = -1.0*kuhn()*kuhn()*ds_/6.0;
      // std::cout << "factor      = " << factor << std::endl;
      int i;
      for (iter.begin(); !iter.atEnd(); ++iter) {
         i = iter.rank(); 
         G = iter.position();
         Gmin = shiftToMinimum(G, mesh().dimensions(), unitCell);
         Gsq = unitCell.ksq(Gmin);
         expKsq_[i] = exp(Gsq*factor);
         expKsq2_[i] = exp(Gsq*factor*0.5);
         //std::cout << i    << "  " 
         //         << Gmin << "  " 
         //          << Gsq  << "  "
         //          << expKsq_[i] << std::endl;
      }
   }
      
   /*
   * Setup the contour length step algorithm.
   */
   template <int D>
   void 
   Block<D>::setupSolver(Block<D>::WField const& w)
   {
      // Preconditions
      int nx = mesh().size();
      UTIL_CHECK(nx > 0);
      
      // Populate expW_
      int i;
      // std::cout << std::endl;
      for (i = 0; i < nx; ++i) {
         expW_[i] = exp(-0.5*w[i]*ds_);
         expW2_[i] = exp(-0.5*0.5*w[i]*ds_);
         // std::cout << "i = " << i 
         //           << " expW_[i] = " << expW_[i]
         //          << std::endl;
      }

      #if 0
      MeshIterator<D> iter;
      IntVec<D> G;
      IntVec<D> Gmin;
      double Gsq;
      double factor = -1.0*kuhn()*kuhn()*ds_/6.0;
      // std::cout << "kDimensions = " << kMeshDimensions_ << std::endl;
      // std::cout << "factor      = " << factor << std::endl;
      iter.setDimensions(kMeshDimensions_);
      for (iter.begin(); !iter.atEnd(); ++iter) {
         i = iter.rank(); 
         G = iter.position();
         Gmin = shiftToMinimum(G, mesh().dimensions(), unitCell);
         Gsq = unitCell.ksq(Gmin);
         expKsq_[i] = exp(Gsq*factor);
         expKsq2_[i] = exp(Gsq*factor*0.5);
         //std::cout << i    << "  " 
         //          << Gmin << "  " 
         //          << Gsq  << "  "
         //          << expKsq_[i] << std::endl;
      }
      #endif
      
   }

   /*
   * Integrate to calculate monomer concentration for this block
   */
   template <int D>
   void Block<D>::computeConcentration(double prefactor)
   {
      // Preconditions
      int nx = mesh().size();
      UTIL_CHECK(nx > 0);
      UTIL_CHECK(ns_ > 0);
      UTIL_CHECK(ds_ > 0);
      UTIL_CHECK(propagator(0).isAllocated());
      UTIL_CHECK(propagator(1).isAllocated());
      UTIL_CHECK(cField().capacity() == nx) 

      // Initialize cField to zero at all points
      int i;
      for (i = 0; i < nx; ++i) {
         cField()[i] = 0.0;
      }

      Propagator<D> const & p0 = propagator(0);
      Propagator<D> const & p1 = propagator(1);

      // Evaluate unnormalized integral
      for (i = 0; i < nx; ++i) {
         cField()[i] += 0.5*p0.q(0)[i]*p1.q(ns_ - 1)[i];
      }
      for (int j = 1; j < ns_ - 1; ++j) {
         for (i = 0; i < nx; ++i) {
            cField()[i] += p0.q(j)[i]*p1.q(ns_ - 1 - j)[i];
         }
      }
      for (i = 0; i < nx; ++i) {
         cField()[i] += 0.5*p0.q(ns_ - 1)[i]*p1.q(0)[i];
      }

      // Normalize
      prefactor *= ds_;
      for (i = 0; i < nx; ++i) {
         cField()[i] *= prefactor;
      }
      #if 0
      #endif

   }

   /*
   * Integrate to Stress exerted by the chain for this block
   */
   template <int D>
   void Block<D>::computeStress(Basis<D>& basis, double prefactor)
   {   
      // Preconditions
      int nx = mesh().size();
      UTIL_CHECK(nx > 0); 
      UTIL_CHECK(ns_ > 0); 
      UTIL_CHECK(ds_ > 0); 
      UTIL_CHECK(propagator(0).isAllocated());
      UTIL_CHECK(propagator(1).isAllocated());

      double dels, normal, increment;
      int r,c;
 
      normal = 3.0*6.0;

      r = (basis.nStar())-1;
      c = 6;

      DArray<double> q1s;
      DArray<double> q2s;
      DArray<double> q3s;
      q1s.allocate(basis.nStar());
      q2s.allocate(basis.nStar());
      q3s.allocate(basis.nStar());

      FArray<double, 6> dQ;

      //FArray<double, basis.nStar()> q1s;
      //FArray<double, basis.nStar()> q2s;

      // Initialize work array and pStress to zero at all points
      int i;
      for (i = 0; i < 6; ++i) {
         dQ [i] = 0.0;
         pStress [i] = 0.0;
      }   

      Propagator<D> const & p0 = propagator(0);
      Propagator<D> const & p1 = propagator(1);

      //Evaluate unnormalized integral   
      for (int j = 0; j < ns_ ; ++j) {
           
          // fft_.forwardTransform(p0.q(j), q1);
         
           q1p = p0.q(j);
           fft_.forwardTransform(q1p, q1);
           basis.convertFieldDftToComponents(q1 , q1s);
           
           //fft_.forwardTransform(p1.q(ns_ - 1 - j), q2);  
           
           q2p = p1.q(ns_ - 1 - j);
           //q2p = p1.q(j);
           fft_.forwardTransform(q2p, q2); 
           basis.convertFieldDftToComponents(q2 , q2s); 

          // for (int d = 0; d < c; ++d){
            //  std::cout<<"q1s("<<d<<")="<<"\t"<<q1s[d]<<"\n";
             // std::cout<<"q2s("<<d<<")="<<"\t"<<q2s[d]<<"\n";
          // }
           


           dels = ds_;

           if (j != 0 && j != ns_ - 1) {
              if (j % 2 == 0) {
                 dels = dels*4.0;
              } else {
                 dels = dels*2.0;
              }           
           }

           for (int n = 0; n < r ; ++n) {
              increment = 0;

              for (int m = 0; m < c ; ++m) {
                 q3s [m] = q1s [m] * basis.star(n).dEigen[m];  
                 increment += q3s [m]*q2s [m]; 
              }
              increment = (increment * kuhn() * kuhn() * dels)/normal;
              dQ [n] = dQ[n]-increment; 
           }    
      }   
      
      // Normalize
      for (i = 0; i < r; ++i) {
         pStress[i] = pStress[i] - (dQ[i] * prefactor);
      }   
      //#if 0
      //#endif

   }

   /*
   * Propagate solution by one step.
   */
   template <int D>
   void Block<D>::step(const QField& q, QField& qNew)
   {
      // Check real-space mesh sizes`
      int nx = mesh().size();
      UTIL_CHECK(nx > 0);
      UTIL_CHECK(q.isAllocated());
      UTIL_CHECK(qNew.isAllocated());
      UTIL_CHECK(q.capacity() == nx);
      UTIL_CHECK(qNew.capacity() == nx);
      UTIL_CHECK(qr_.capacity() == nx);
      UTIL_CHECK(expW_.capacity() == nx);

      // Fourier-space mesh sizes
      int nk = qk_.capacity();
      UTIL_CHECK(expKsq_.capacity() == nk);

      // Apply pseudo-spectral algorithm
      int i;
      for (i = 0; i < nx; ++i) {
         qr_[i] = q[i]*expW_[i];
         qr2_[i] = q[i]*expW2_[i];
      }
      fft_.forwardTransform(qr_, qk_);
      fft_.forwardTransform(qr2_, qk2_);
      for (i = 0; i < nk; ++i) {
         qk_[i][0] *= expKsq_[i];
         qk_[i][1] *= expKsq_[i];
         qk2_[i][0] *= expKsq2_[i];
         qk2_[i][1] *= expKsq2_[i];
      }
      fft_.inverseTransform(qk_, qr_);
      fft_.inverseTransform(qk2_, qr2_);
      for (i = 0; i < nx; ++i) {
         qf_[i] = qr_[i]*expW_[i];
         qr2_[i] = qr2_[i]*expW_[i];
      }

      fft_.forwardTransform(qr2_, qk2_);
      for (i = 0; i < nk; ++i) {
         qk2_[i][0] *= expKsq2_[i];
         qk2_[i][1] *= expKsq2_[i];
      }
      fft_.inverseTransform(qk2_, qr2_);
      for (i = 0; i < nx; ++i) {
         qr2_[i] = qr2_[i]*expW2_[i];
      }
      for (i = 0; i < nx; ++i) {
         qNew[i] = (4.0*qr2_[i] - qf_[i])/3.0;
      }
   }

}
}
#endif
