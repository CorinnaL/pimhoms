#include <functional>
#include "CLinComb.hpp"
extern bool latex;
using Hom = std::function<CLinComb (CLinComb)>;
std::function<CLinComb (CLinComb)> operator*(std::function<CLinComb (CLinComb)> hom1,
    std::function<CLinComb (CLinComb)> hom2);
void checkThree(std::function<CLinComb (CLinComb)> hom1,
        std::function<CLinComb (CLinComb)> hom2,
        std::function<CLinComb (CLinComb)> hom3,
        CLinComb x,std::string firstcomp, std::string secondcomp,
        std::string thirdcomp);
void checkTwo(Hom hom1,
      Hom  hom2,
      CLinComb x,std::string firstcomp, std::string secondcomp);
