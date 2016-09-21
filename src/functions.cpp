#include <string>
#include "CTensor.hpp"
#include "CSymm.hpp"
#include "functions.hpp"
#include <iostream>

std::function<CLinComb (CLinComb)> operator*(std::function<CLinComb (CLinComb)> hom1,
    std::function<CLinComb (CLinComb)> hom2)
{
  return [hom1,hom2] (CLinComb a) {return hom1(hom2(a));};
}
void checkThree(std::function<CLinComb (CLinComb)> hom1,
        std::function<CLinComb (CLinComb)> hom2,
        std::function<CLinComb (CLinComb)> hom3,
        CLinComb x,std::string firstcomp, std::string secondcomp,
        std::string thirdcomp)
{
    if(latex)
    {
        std::cout << ""<<firstcomp <<"\\\\"<<std::endl;
        std::cout << x.print() << std::endl;
        CLinComb c1 = hom1(x);
        std::cout << ""<<secondcomp <<"\\\\"<<std::endl;
        std::cout << x.print() << std::endl;
        CLinComb c2 = hom2(x);
        std::cout << ""<<thirdcomp <<"\\\\"<<std::endl;
        std::cout << x.print() << std::endl;
        CLinComb c3 = hom3(x);
        std::cout << "\\end{align*}\n\\begin{align*}"<<std::endl;
    }
    else
    {
        std::cout << x.print() << std::endl;
        std::cout << "***************"<<firstcomp<<"***************" <<std::endl;
        CLinComb c1 = hom1(x);
        std::cout << "***************"<<secondcomp<<"***************" <<std::endl;
        CLinComb c2 = hom2(x);
        std::cout << "***************"<<thirdcomp<<"***************" <<std::endl;
        CLinComb c3 = hom3(x);
        if(c1.print() != c2.print() || c1.print() != c3.print())
        {
            std::cout << "----------------------------------------------------------"<<std::endl;
            std::cout << "!!!!!!!Check!!!!!!!!!"<<std::endl;
        }
        std::cout << "------------------------------------------------------------"<<std::endl;
    }
}

void checkTwo(Hom hom1,
      Hom  hom2,
      CLinComb x,std::string firstcomp, std::string secondcomp)
{
    if(latex)
    {
        std::cout << ""<<firstcomp <<"\\\\"<<std::endl;
        std::cout << x.print()<<std::endl;
        CLinComb c1 = hom1(x);
        std::cout << ""<<secondcomp <<"\\\\"<<std::endl;
        std::cout << x.print()<<std::endl;
        CLinComb c2 = hom2(x);
        std::cout << "\\end{align*}\n\\begin{align*}"<<std::endl;
    }
    else
    {
        std::cout << x.print() << std::endl;
        std::cout << "***************"<<firstcomp<<"***************" <<std::endl;
        CLinComb c1 = hom1(x);
        std::cout << "***************"<<secondcomp<<"***************" <<std::endl;
        CLinComb c2 = hom2(x);
        if(c1.print() != c2.print())
        {
          std::cout << "----------------------------------------------------------"<<std::endl;
          std::cout << "!!!!!!!Check!!!!!!!!!"<<std::endl;
        }
        std::cout << "------------------------------------------------------------"<<std::endl;
    }
}

std::string findreplace(std::string input, std::string find, std::string replace)
{
    std::size_t findsize = find.size();
    while(input.find(find)!=std::string::npos)
    {
        input.replace(input.find(find),findsize,replace);
    }
    return input;
}


CTensor operator *(CTensor l,CSymm r)
{
    std::string symmstring = r.get();
    CTensor result = l;
    while(symmstring!="")
    {
        if(symmstring.front()=='s')
        {
            result.setComponents(result.second(), result.third(), result.first());
        }
        else if(symmstring.front()=='t')
        {
            result.setComponents(result.second(), result.first(), result.third());
        }
        symmstring = symmstring.substr(1,symmstring.size());
    }
    return result;
}

