#include <pssp/basis/Basis.h>

using namespace Pscf;
using namespace Util;
using namespace Pscf::Pssp;

int dimer[3] = {10, 10, 10};

template <int D>
int cheat_helper(int N, IntVec<D>& vector)
{
  if(N - 1 == 0)
    return 0;
  else {
    IntVec<D-1> temp;
    for(int i = 0; i < D - 1; i++) {
      temp[i] = vector[i];
    }
    return vector[N-1] + cheat_helper(N-1, temp) * dimer[N-1];
  }
}
template <int D>
int island(IntVec<D> monkey){
  if(D-1 == 0)
    return 0;
  return cheat_helper(D, monkey);
}

int main(){
  IntVec<3> monkey;
  monkey[0] = 3;
  monkey[1] = 4;
  monkey[2] = 5;
  std::cout<<island(monkey)<<std::endl;
  return 0;
}
