#ifndef PSSP_BASIS_TPP
#define PSSP_BASIS_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Basis.h"
#include <pscf/crystal/shiftToMinimum.h>
#include <vector>


namespace Pscf {
namespace Pssp
{

   /*
   * Constructor.
   */
   template <int D>
   Basis<D>::Basis()
    : nWave_(0), 
      nStar_(0), 
      unitCellPtr_(0), 
      meshPtr_(0)
   {}

   /*
   * Construct basis for pseudo-spectral scft.
   */
   template <int D>
   void Basis<D>::makeBasis(const Mesh<D>& mesh, 
                            const UnitCell<D>& unitCell,
                            std::string groupName)
   {
      // Save pointers to mesh and unit cell
      meshPtr_ = &mesh;
      unitCellPtr_ = &unitCell;

      // Allocate arrays
      nWave_ = mesh.size();
      waves_.allocate(nWave_);
      waveId_.allocate(nWave_); 

      // Make sorted array of waves
      makeWaves();

      // Create identity group
      // Temporary stub: Replace with code to read a group from file,
      // or pass the group to makeBasis and generate it elsewhere.
      SpaceGroup<D> group;
      group.makeCompleteGroup();

      // Identify stars of waves that are related by symmetry
      makeStars(group);

   }

   /*
   * Construct ordered list of waves.
   * 
   * On exit:
   *  - Array waves_ contains list of waves ordered by sqNorm.
   *  - Each wave has indicesDft, indicesBz and sqNorm set.
   *  - stars_ array is still empty.
   */
   template <int D>
   void Basis<D>::makeWaves()
   {
      IntVec<D> meshDimensions = mesh().dimensions();

      // Generate dft mesh of waves, store in local std::vector twaves_
      std::vector<NWave> twaves;
      twaves.reserve(nWave_);
      {
         NWave w;
         IntVec<D> v;

         // Loop over dft mesh to generate all waves
         MeshIterator<D> itr(mesh().dimensions());
         for (itr.begin(); !itr.atEnd(); ++itr) {
            w.indicesDft = itr.position();
            v = shiftToMinimum(w.indicesDft, meshDimensions, *unitCellPtr_);
            w.indicesBz = v;
            w.sqNorm = unitCellPtr_->ksq(v);
            twaves.push_back(w);
         }
      }

      // Sort temporary container twaves by wavevector norm
      {
         // Define function object that sorts based on sqNorm
         struct NormComparator{
            bool operator () (NWave a, NWave b)
            {  return (a.sqNorm < b.sqNorm); }
         };
         // Sort twaves
         NormComparator comp;
         std::sort(twaves.begin(), twaves.end(), comp);
      }

      // Copy temporary array twaves into member variable waves_
      for (int i = 0; i < nWave_; ++i) {
         waves_[i].sqNorm = twaves[i].sqNorm;
         waves_[i].indicesDft = twaves[i].indicesDft;
         waves_[i].indicesBz = twaves[i].indicesBz;
      }

   }

   template <int D>
   bool Basis<D>::NWaveComp::operator() (Basis::NWave a, Basis::NWave b) const
   {
       if (a.indicesDft == b.indicesDft) {
          return false;
       } else 
       if (a.indicesBz > b.indicesBz) {
          return true;
       } else
       if (a.indicesBz < b.indicesBz) {
          return false;
       } else {
          return (a.indicesDft < b.indicesDft);
       }
   }

   template <int D>
   void Basis<D>::makeStars(const SpaceGroup<D>& group)
   {
      std::set<NWave, NWaveComp> list;  // set of vectors of equal norm
      std::set<NWave, NWaveComp> star;  // set of symmetry-related vectors 
      GArray<NWave> work;               // temporary list, ordered by star
      typename std::set<NWave, NWaveComp>::iterator rootItr;
      typename std::set<NWave, NWaveComp>::iterator setItr;
      NWave wave;
      Star newStar;

      double Gsq;
      double Gsq_max = 1.0;
      double epsilon = 1.0E-8;
      IntVec<D> rootVec;
      IntVec<D> vec;
      int listId = 0;      // id for this list
      int listBegin = 0;   // id of first wave in this list
      int listEnd = 0;     // (id of last wave in this list) + 1
      int starId = 0;      // id for this star
      int starBegin = 0;   // id of first wave in this star
      int starSize = 0;    // size (number of waves) in this star
      int i, j;
      for (i = 0; i <= nWave_; ++i) {

         // Determine if this wave begins a new list
         bool newList = false;
         if (i == nWave_) {
            listEnd = i;
            newList = true;
         } else {
            Gsq = waves_[i].sqNorm;
            if (Gsq > Gsq_max + epsilon) {
               Gsq_max = Gsq;
               listEnd = i;
               newList = true;
            }
         }

         // Process new list
         if (newList && listEnd > 0) {

            // Copy waves of equal norm into container "list"
            list.clear();
            for (j = listBegin; j < listEnd; ++j) {
               wave.indicesDft = waves_[j].indicesDft;
               wave.indicesBz = waves_[j].indicesBz;
               wave.sqNorm = waves_[j].sqNorm;
               list.insert(wave);
            }

            // Loop over stars within this list
            work.clear();
            IntVec<D> nVec;
            rootItr = list.begin();
            int nextInvert = 1;
            while (list.size() > 0) {

               // Construct a star from root wave *rootItr.
               rootVec = rootItr->indicesBz;
               Gsq = unitCell().ksq(rootVec);
               star.clear();
               for (j = 0; j < group.size(); ++j) {
                  vec = rootVec*group[j];
                  wave.sqNorm = Gsq;
                  wave.indicesBz = vec;
                  wave.indicesDft = vec;
                  mesh().shift(wave.indicesDft);
                  star.insert(wave);
                  work.append(wave);
                  list.erase(wave);
               }
               starSize = star.size();

               // TODO: Replace above by code that generates
               // a star for an arbitrary space group.

               // Process the star
               newStar.beginId = starBegin;
               newStar.eigen = Gsq;
               if (nextInvert == -1) {

                  // If this star is second of pair related by symmetry,
                  // then set root of next to beginning of remaining list.

                  newStar.invertFlag = -1;
                  rootItr = list.begin();
                  nextInvert = 1;

               } else {

                  // If this star is not the second of a pair, 
                  // then determine if it is closed under inversion.

                  // Compute negation of root vector, shift to DFT mesh
                  nVec.negate(rootItr->indicesBz);
                  (*meshPtr_).shift(nVec);
       
                  // Search for negation of root vector within star.
                  bool negationFound = false;
                  setItr = star.begin(); 
                  for ( ; setItr != star.end(); ++setItr) {
                     if (nVec == setItr->indicesDft) {
                        negationFound = true;
                        break;
                     }
                  }

                  if (negationFound) {

                     // If star is closed under inversion, then root
                     // of next star is first vector of remaining list.

                     newStar.invertFlag = 0;
                     rootItr = list.begin();
                     nextInvert = 1;

                  } else {

                     newStar.invertFlag = 1;
                     nextInvert = -1;

                     // If start is not closed, find negation of root in
                     // the remaining list, use negation as root of the
                     // next star.

                     setItr = list.begin(); 
                     for ( ; setItr != list.end(); ++setItr) {
                        if (nVec == setItr->indicesDft) {
                           negationFound = true;
                           rootItr = setItr;
                           break;
                        }
                     }
                     // On exit after break, rootItr = &nVec
                     if (!negationFound) {
                        UTIL_THROW("Negative of root vector not found");
                     }
                  }

               }
               newStar.size = star.size();
               newStar.endId = newStar.beginId + star.size();
               stars_.append(newStar);

               ++starId;
               starBegin = newStar.endId;
            }
            UTIL_CHECK(list.size() == 0);
            UTIL_CHECK(work.size() == listEnd - listBegin);

            // Copy work container into corresponding section of waves_
            int k;
            for (j = 0; j < work.size(); ++j) {
               k = j + listBegin;
               waves_[k].indicesDft = work[j].indicesDft;
               waves_[k].indicesBz = work[j].indicesBz;
               waves_[k].sqNorm = work[j].sqNorm;
            }

            ++listId;
            listBegin = listEnd;
         }
      }
      nStar_ = stars_.size();

      // Final processing of stars
      int waveId;
      for (i = 0; i < nStar_; ++i) {

         // Set characteristic wavevector waveBz for each star
         if (stars_[i].invertFlag == -1) {
           waveId = stars_[i].endId - 1;
         } else {
           waveId = stars_[i].beginId;
         }
         stars_[i].waveBz = waves_[waveId].indicesBz;

         // Set starId for all associated waves
         for (j = stars_[i].beginId; j < stars_[i].endId; ++j) {
            waves_[j].starId = i;
            waves_[j].coeff = std::complex<double>(1.0, 0.0);
         }

         // Set coeff for all associated waves
         for (j = stars_[i].beginId; j < stars_[i].endId; ++j) {
            waves_[j].coeff = std::complex<double>(1.0, 0.0);
         }
      }

      #if 1
      // Output all waves
      std::cout << std::endl;
      std::cout << "Waves:" << std::endl;
      for (i = 0; i < nWave_; ++i) {
         std::cout << Int(i,4);
         std::cout << Int(waves_[i].starId, 4) << " |";
         for (j = 0; j < D; ++j) {
            std::cout << Int(waves_[i].indicesBz[j], 4);
         }
         std::cout << " | " << Dbl(waves_[i].sqNorm, 12);
         std::cout << std::endl;
      }
      #endif
 
      #if 1 
      // Output all stars
      std::cout << std::endl;
      std::cout << "Stars:" << std::endl;
      for (i = 0; i < nStar_; ++i) {
         std::cout << i 
                   << Int(stars_[i].beginId, 3)
                   << Int(stars_[i].endId, 3)
                   << Int(stars_[i].invertFlag, 3);
         std::cout << " |";
         for (j = 0; j < D; ++j) {
            std::cout << Int(stars_[i].waveBz[j], 4);
         }
         std::cout << " | " << Dbl(stars_[i].eigen, 12);
         std::cout << std::endl;
      }
      #endif
  
   }

   template <int D>
   void Basis<D>::convertFieldComponentsToDft(DArray<double>& components, 
                                              RFieldDft<D>& dft)
   {
      std::complex<double> coeff;
      double z;
      double z1;
      double z2;
      IntVec<D> G2;
      int partnerId;
      //loop through a factor of 2 more indices. 
      //Difficult to optimize without specialization or making new objects
      for (int i = 0; i < meshPtr_->size(); ++i) {
         
         //if last indice of waves is not > n3/2+1
         IntVec<D> waveId = meshPtr_->position(i);
         int starId = wave(waveId).starId;
         if (!wave(waveId).implicit) {
            coeff = wave(waveId).coeff;
            //determining the rank of RFieldDft
            IntVec<D> offsets;
            offsets[D-1] = 1;
            for (int j = D-1; j > 0; --j) {
               if (j == D-1) {
                  offsets[j-1] = offsets[j]*(meshPtr_->dimension(j)/2 + 1);
               } else {
                  offsets[j-1] = offsets[j]*meshPtr_->dimension(j);
               }
            }
            int rank = 0;
            for (int j = 0; j < D; ++j) {
               rank += wave(waveId).indicesDft[j] * offsets[j];
            }
            //std::cout<<rank<<std::endl;

            switch (stars_[starId].invertFlag) {
               case 0 :
                  z = components[starId];
                  coeff = z*coeff;
                  dft[rank][0] = coeff.real();
                  dft[rank][1] = coeff.imag();
                  break;
               case 1 :
                  z1 = components[starId];
                  for(int j = 0; j < D; ++j){
                     G2[j] = -waveId[j];
                  }
                  meshPtr_->shift(G2);
                  partnerId = meshPtr_->rank(G2);
                  z2 = components[partnerId];
                  /*//debug code
                  if(stars_[partnerId].invertFlag != -1) {
                     std::cout<<"This star is "<<starId<<" and my partner is "
                     <<partnerId<<std::endl;
                     std::cout << "I have invertFlag of "
                               << stars_[starId].invertFlag
                               << " and he has "
                               << stars_[partnerId].invertFlag<<std::endl;
                     std::cout << " WaveBz "<<stars_[starId].waveBz 
                               << " and he "
                               << stars_[partnerId].waveBz<<std::endl;
                     std::cout << "Wave dft indices"
                               << waves_[stars_[starId].beginId].indicesDft
                               << "and he is" 
                               << waves_[stars_[partnerId].beginId].indicesDft
                               << std::endl;
                     std::cout << "My previous wave has implicit of "
                               << waves_[stars_[starId].beginId - 1].implicit
                               << std::endl;
                     std::cout<<"And position of "
                     <<waves_[stars_[starId].beginId-1].indicesDft<<std::endl;
                  }
                  //end debug*/
                  UTIL_CHECK(stars_[partnerId].invertFlag == -1);
                  coeff = std::complex<double>(z1,-z2)*coeff/sqrt(2);
                  dft[rank][0] = coeff.real();
                  dft[rank][1] = coeff.imag();
                  break;
               case -1 :
                  z1 = components[starId];
                  for(int j = 0; j < D; ++j){
                     G2[j] = -waveId[j];
                  }
                  meshPtr_->shift(G2);
                  partnerId = meshPtr_->rank(G2);
                  z2 = components[partnerId];
                  UTIL_CHECK(stars_[partnerId].invertFlag == 1);
                  coeff = std::complex<double>(z2, z1) * coeff / sqrt(2);
                  dft[rank][0] = coeff.real();
                  dft[rank][1] = coeff.imag();
                  break;
            }

         } else {
            //do nothing. Dft does not need allocation for this
         }
      }
   }

   template <int D>
   void Basis<D>::convertFieldDftToComponents(RFieldDft<D>& dft, 
                                              DArray<double>& components)
   {
      //if its not too distrupting maybe use consistent names for logical size?
      int nStars = nStar_;
      UTIL_CHECK(nStars == components.capacity());
      IntVec<D> indiceBz;
      std::complex<double> coeff;
      std::complex<double> temp;
      fftw_complex z;
      for (int i = 0; i < nStars; ++i) {
         indiceBz = stars_[i].waveBz;
         coeff = wave(indiceBz).coeff;

         //check if wave exists in dft
         if (!wave(indiceBz).implicit) {
            meshPtr_->shift(indiceBz);
            //std::cout<<"Yoohoo"<<std::endl;
            //might be good to rethink how elements of RFieldDft is accessed
            IntVec<D> offsets;
            offsets[D-1] = 1;
            for (int j = D-1; j > 0; --j) {
               if(j == D-1) {
                  offsets[j-1] = offsets[j]*(meshPtr_->dimension(j)/2 + 1);
               } else {
                  offsets[j-1] = offsets[j]*meshPtr_->dimension(j);
               }
            }
            int rank = 0;
            for (int j = 0; j < D; j++) {
               rank += indiceBz[j] * offsets[j];
            }
            //std::cout<<offsets<<std::endl;
            //std::cout<<indiceBz<<std::endl;
            //std::cout<<rank<<std::endl;
            z[0] = dft[rank][0];
            z[1] = dft[rank][1];
            //Note: the code implies that we do not need indicesDft from Wave? 
         } else {
         //wave does not exists in dft. have to find explicit pair
            meshPtr_->shift(indiceBz);
            for (int j = 0; j < D; ++j) {
               indiceBz[j] = -indiceBz[j];
            }
            meshPtr_->shift(indiceBz);


            IntVec<D> offsets;
            offsets[D-1] = 1;
            for (int j = D-1; j > 0; --j) {
               if(j == D-1) {
                  offsets[j-1] = offsets[j]*(meshPtr_->dimension(j)/2 + 1);
               } else {
                  offsets[j-1] = offsets[j]*meshPtr_->dimension(j);
               }
            }

            int rank = 0;
            for (int j = 0; j < D; j++) {
               rank += indiceBz[j] * offsets[j];
            }
            //std::cout<<"Yoohoo"<<std::endl;
            //std::cout<<rank<<std::endl;
            z[0] = dft[rank][0];
            z[1] = dft[rank][1];
            z[1] = -z[1];
         }
         
         //std::cout<<i<<std::endl;
         //reintepret cast is not used since coding standards is old
         temp = std::complex<double>(z[0],z[1]);
         //std::cout<<temp<<std::endl;
         //assign value to components
         switch(stars_[i].invertFlag) {
            case 0 :
               components[i] = (temp/coeff).real();
               break;
            case 1 :
               components[i] = (temp/coeff).real() * sqrt(2);
               break;
            case -1 :
               components[i] = (temp/coeff).imag() * sqrt(2);
               break;
         }
      }
   }

   template <int D>
   int Basis<D>::nBasis() const
   {
      int count = 0;
      for (int i = 0; i < stars_.capacity(); ++i) {
         if (!stars_[i].cancel) {
            count++;
         }
      }
      return count;
   }

   template <int D>
   typename Basis<D>::Wave& Basis<D>::wave(IntVec<D> vector)
   {
      if (!meshPtr_->isInMesh(vector)) {
         meshPtr_->shift(vector);
      }
      return waves_[waveId_[meshPtr_->rank(vector)]];
   }

   template <int D>
   void Basis<D>::makedksq(const UnitCell<D>& unitCell)
   {

      unitCellPtr_ = &unitCell;

      // Initialize all elements to zero
      int i, j, p, q;
      for (i = 0; i < unitCellPtr_->nParams(); ++i) {
         for (j=0; j < nStar_; ++j){
            dksq(i,j)=0.0;
         }
      }
	
      for (i = 0; i < unitCellPtr_->nParams(); ++i) {
         for (j=0; j < nStar_; ++j){
            for (p=0; p < D; ++p){
               for (q=0; q < D; ++q){

                  dksq(i,j) = dksq(i,j) + (stars_[j].waveBz[p]*stars_[j].waveBz[q]*(unitCellPtr_->dkkBasis(i, p,q)));
  
               }
            }
         }
      }

   }

}
}

#endif
