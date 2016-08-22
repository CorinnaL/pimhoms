#include "CSymm.hpp"
#include <iostream>
#include <sstream>

std::string findreplace(std::string indput, std::string find, std::string replace);

CSymm::CSymm()
{
    //ctor
}
void CSymm::Simplify()
{
    element = findreplace(element, "st","tss");
    element = findreplace(element, "sss","");
    element = findreplace(element, "tt","");
}

const CSymm operator *(const CSymm& l,const CSymm& r)
{
    CSymm result(l.get()+r.get());
    return result;
}

void CSymm::setElement(std::string val)
{
    element = val;
    Simplify();
}

CSymm::CSymm(std::string s):element(s)
{
    Simplify();
}

CSymm::CSymm(char const* s):element(s)
{
    Simplify();
}

bool operator <(const CSymm& l,const CSymm& r)
{
    return l.get()<r.get();
}

bool operator ==(const CSymm& l,const CSymm& r)
{
    return l.get()==r.get();
}

std::string CSymm::print() const
{
    if(element == "") return "id";
    else if(element == "t") return "(12)";
    else if(element == "s") return "(123)";
    else if(element == "ts") return "(23)";
    else if(element == "ss") return "(132)";
    else if(element == "tss") return "(13)";
    else return "invalid";
}

std::list<CSymm> factorize(const CSymm& r)
{
    std::list<CSymm> out;
    for(char c : r.get())
    {
      std::stringstream ss;
      std::string s;
      ss << c;
      ss >> s;
        out.push_back(CSymm(s));
    }
    return out;
}

int signum(const CSymm& r)
{
  int sg = 1;
  for(const auto& fac : factorize(r))
  {
    if(fac == "t")
    {
      sg*= -1;
    }
  }
  return sg;
}

CSymm inverse(const CSymm& l)
{
  // this is a little hack using that all elements in S3\A3 are transpositions
  if(signum(l) == -1 || l=="") return l;
  if(l=="s") return "ss";
  else return "s";
}
