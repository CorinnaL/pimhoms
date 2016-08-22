#include <assert.h>

#include "CWreathModuleElement.hpp"
#include <iostream>

CTensor operator *(CTensor l,CSymm r);

CWreathModuleElement::CWreathModuleElement(subgroup type, CSymm coset, int multiplicity,
    CTensor tensor, specht lambda, spechtbasis x=spechtbasis::x1):
    mtype(type), mcoset(coset), mmultiplicity(multiplicity), mtensor(tensor),
    mlambda(lambda), mx(x)
{
    this->setCoset(coset);
}

//CWreathModuleElement::CWreathModuleElement(subgroup type):mtype(type), msignum(false)
//{
//    mtensor = CTensor("","","");
//    mcoset = "";
//    mmultiplicity = 1;
//}

std::string CWreathModuleElement::print() const
{
    if(mtensor.isZero()|| mmultiplicity == 0)
    {
        return "";
    }
    std::string output;
    if(this->multiplicity() == -1) output +="-";
    else if(this->multiplicity()!=1) output+=std::to_string(mmultiplicity);
    output+=mtensor.print();
    //output+="\otimes"+mcoset.print();
    output+=mcoset.print();
    if(mlambda == specht::r)
    {
      output+= "_";
      output+=mx==spechtbasis::x1?"x1":"x2";
    }
    return output;
}

CSymm subgroupToSymm(subgroup x)
{
  switch(x)
  {
    case subgroup::i:
      return CSymm("");
    case subgroup::t:
      return CSymm("t");
    case subgroup::ts:
      return CSymm("ts");
    case subgroup::tss:
      return CSymm("tss");
    default:
      {
      std::cout << "weird stuff" << std::endl;
      return "";
      }
  }
}

void CWreathModuleElement::setCoset(CSymm coset)
{
  mcoset = coset;
}

int subgroupSize(subgroup sub)
{
  if(sub == subgroup::i)
    return 1;
  if(sub == subgroup::S3)
    return 6;
  if(sub == subgroup::s || sub == subgroup::ss)
    return 3;
  else
    return 2;
}

bool operator <(CWreathModuleElement lhs, CWreathModuleElement rhs)
{
  if(lhs.coset() < rhs.coset()) return true;
  if(rhs.coset() < lhs.coset()) return false;
  if(lhs.x() == spechtbasis::x1 && rhs.x() == spechtbasis::x2) return true;
  if(lhs.x() == spechtbasis::x2 && rhs.x() == spechtbasis::x1) return false;
  if(lhs.tensor() < rhs.tensor()) return true;
  else return false;
}

subgroup symmToSubgroup(CSymm sy)
{
    if(sy=="")
      return subgroup::i;
    if(sy=="t")
      return subgroup::t;
    if(sy=="ts")
      return subgroup::ts;
    if(sy=="tss")
      return subgroup::tss;
    else
      {
      std::cout << "weird stuff" << std::endl;
      return subgroup::s;
      }
}
