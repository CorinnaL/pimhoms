#include "CTensor.hpp"
extern bool latex;

CTensor::CTensor(std::string first, std::string second, std::string third):
    firstcomp(first), secondcomp(second),thirdcomp(third)
{

}

CTensor::CTensor():
    firstcomp(""), secondcomp(""),thirdcomp("")
{

}
void CTensor::setComponents(std::string first, std::string second, std::string third)
{
    firstcomp = first;
    secondcomp = second;
    thirdcomp = third;
}

std::string CTensor::print() const
{
  if(isZero()) return "0";
  if(latex)
  {
      return "("+firstcomp+"\\otimes "+secondcomp+"\\otimes "+thirdcomp+")";
  }
  else
  {
      return "("+firstcomp+","+secondcomp+","+thirdcomp+")";
  }
}

bool operator <(const CTensor l,const CTensor r)
{
    if(l.first()<r.first()) return true;
    else if(l.first()>r.first()) return false;
    else if(l.second()<r.second()) return true;
    else if(l.second()>r.second()) return false;
    else if(l.third()<r.third()) return true;
    else if(l.third()>r.third()) return false;
    else return false;
}

bool operator ==(const CTensor l, const CTensor r)
{
    return !(l<r) && !(r<l);
}
