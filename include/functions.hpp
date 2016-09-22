#include <functional>
#include "CLinComb.hpp"


enum class comptype {iij,ijj};

extern bool latex;
using Hom = std::function<CLinComb (CLinComb)>;
Hom operator*(Hom hom1,
    Hom hom2);
void checkThree(Hom hom1, Hom hom2, Hom hom3, CLinComb x,
        std::string i, std::string j, std::string k,
        std::string firstlam, std::string secondlam1, std::string secondlam2,
        std::string secondlam3, std::string thirdlam, comptype twoitwoj);
void checkTwo(Hom hom1, Hom  hom2, CLinComb x,
        std::string i, std::string j, std::string k, std::string firstlam,
        std::string secondlam1,std::string secondlam2, std::string thirdlam);

std::string makecalctitleij(std::string i,std::string j,std::string k,
        std::string firstlam,std::string secondlam, std::string thirdlam);
std::string makecalctitleji(std::string i,std::string j,std::string k,
        std::string firstlam,std::string secondlam, std::string thirdlam);
