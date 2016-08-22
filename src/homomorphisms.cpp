#include "CLinComb.hpp"
#include "CWreathModuleElement.hpp"
#include "homomorphism.hpp"
std::string homname = "g";
CTensor operator *(CTensor l,CSymm r);


CLinComb changeFirstGen(CLinComb in,subgroup type, specht lambda)
{
  CLinComb output(type,lambda);
  for(const auto& summand : in.summands())
  {
    std::string a = summand.tensor().first();
    std::string b = summand.tensor().second();
    std::string c = summand.tensor().third();
    CSymm pi = summand.coset();
    int m = summand.multiplicity();
    spechtbasis x = summand.x();
    output.CreateAndAdd(homname+a,b,c,pi,m,x);
  }
  std::cout << output.print()<<std::endl;
  return output;
}


CLinComb fi1ij(CLinComb in,specht lam)
  // Here we assume that the input is from iij(2) / iij(1,1).
  // Then there are different cases: either i1 != j,
  // then the result will be of type 1, or i1 == j,
  // then the result will be of type 2. In this function
  // we differentiate between the two cases by passing
  // the type of the output as an argument.
{
  CLinComb output(subgroup::tss,lam);
  if(lam == specht::t1)
  {
    output = CLinComb(subgroup::i,lam);
  }
  for(const auto& summand : in.summands())
  {
    std::string a = summand.tensor().first();
    std::string b = summand.tensor().second();
    std::string c = summand.tensor().third();
    CSymm pi = summand.coset();
    int m = summand.multiplicity();
    switch(in.type())
    {
    case subgroup::t:
      {
        if(in.lambda() == specht::t2)
        {
          output.CreateAndAdd(homname+a,b,c,pi,1*m);
          output.CreateAndAdd(homname+b,a,c,"t"*pi,1*m);
          break;
        }
        if(in.lambda() == specht::s2)
        {
          output.CreateAndAdd(homname+a,b,c,pi,1*m);
          output.CreateAndAdd(homname+b,a,c,"t"*pi,-1*m);
          break;
        }
      }
    default:
      {
        std::cout<<"The element" << in.print() << "is incompatible" <<
          "with the homomorphism fi1ij. Returning input."<<std::endl;
        return in;
      }
    }
  }
  std::cout << output.print()<<std::endl;
  return output;
}

//The function only assumes that the output is of type 1. There are
// two different cases where this can occur. Either the input is of
// type 2 of the form iij with i1!=j or the input is already
// of type 1.
CLinComb fijk(CLinComb in)
{
  switch(in.type())
  {
  case subgroup::i:
    {
      return changeFirstGen(in,in.type(),in.lambda());
    }
  case subgroup::t:
    {
      return fi1ij(in,specht::t1);
    }
  default:
    {
      std::cout<<"The element" << in.print() << "is incompatible" <<
        "with the homomorphism fijk. Returning input."<<std::endl;
      return in;
    }
  }
}

// The output being of type 2 of the form iki happens again in
// two cases. Either the input is of type 1 of the form iji1
// or it is of type 2 of the form iii1.
// The two following functions are first for the output
// in V(i,k,i;(2)) and then in V(i,k,i;(1,1)).
CLinComb fikit(CLinComb in)
{
  switch(in.type())
  {
  case subgroup::i:
    {
      return changeFirstGen(in,subgroup::tss,specht::t2);
    }
  case subgroup::t:
    {
      return fi1ij(in,specht::t2);
    }
  default:
    {
      std::cout<<"The element" << in.print() << "is incompatible" <<
        "with the homomorphism fikit. Returning input."<<std::endl;
      return in;
    }
  }
}

CLinComb fikis(CLinComb in)
{
  switch(in.type())
  {
  case subgroup::i:
    {
      return changeFirstGen(in,subgroup::tss,specht::s2);
    }
  case subgroup::t:
    {
      return fi1ij(in,specht::s2);
    }
  default:
    {
      std::cout<<"The element" << in.print() << "is incompatible" <<
        "with the homomorphism fikis. Returning input."<<std::endl;
      return in;
    }
  }
}

// We get an element of type 2 of the form kii whenever
// the input is of type 3. We also get an element of
// this form if the input is of type kii and k1!=i.
CLinComb fkiis(CLinComb in)
{
  CLinComb output(subgroup::ts,specht::s2);
  for(const auto& summand : in.summands())
  {
    std::string a = summand.tensor().first();
    std::string b = summand.tensor().second();
    std::string c = summand.tensor().third();
    CSymm pi = summand.coset();
    int m = summand.multiplicity();
    switch(in.type())
    {
    case subgroup::ts:
      {
        if(in.lambda() == specht::s2)
        {
          return changeFirstGen(in,in.type(),in.lambda());
        }
      }

    case subgroup::S3:
      {
        if(in.lambda() == specht::r)
        {
          if(summand.x() == spechtbasis::x1)
          {
            output.CreateAndAdd(homname+c,a,b,"s"*pi,2*m);
            output.CreateAndAdd(homname+a,b,c,pi,-1*m);
            output.CreateAndAdd(homname+b,c,a,"ss"*pi,-1*m);
            break;
          }
          if(summand.x() == spechtbasis::x2)
          {
            output.CreateAndAdd(homname+a,b,c,pi,2*m);
            output.CreateAndAdd(homname+b,c,a,"ss"*pi,-1*m);
            output.CreateAndAdd(homname+c,a,b,"s"*pi,-1*m);
            break;
          }
        }

        if(in.lambda() == specht::s3)
        {
          output.CreateAndAdd(homname+a,b,c,""*pi,1*m);
          output.CreateAndAdd(homname+b,c,a,"ss"*pi,1*m);
          output.CreateAndAdd(homname+c,a,b,"s"*pi,1*m);
          break;
        }
        if(in.lambda() == specht::t3)
        {
          std::cout<<"The element" << in.print() << "is incompatible" <<
            "with the homomorphism fkiis. Returning input."<<std::endl;
          return in;
        }
      }
    default:
      {
        std::cout<<"The element" << in.print() << "is incompatible" <<
          "with the homomorphism fkiis. Returning input."<<std::endl;
        return in;
      }
    }
  }
  std::cout << output.print()<<std::endl;
  return output;
}


CLinComb fkiit(CLinComb in)
{
  CLinComb output(subgroup::ts,specht::t2);
  for(const auto& summand : in.summands())
  {
    std::string a = summand.tensor().first();
    std::string b = summand.tensor().second();
    std::string c = summand.tensor().third();
    CSymm pi = summand.coset();
    int m = summand.multiplicity();
    switch(in.type())
    {
    case subgroup::ts:
      {
        if(in.lambda() == specht::t2)
        {
          return changeFirstGen(in,in.type(),in.lambda());
        }
      }
    case subgroup::S3:
      {
        if(in.lambda() == specht::r)
        {
          if(summand.x() == spechtbasis::x1)
          {
            output.CreateAndAdd(homname+a,b,c,pi,1*m);
            output.CreateAndAdd(homname+b,c,a,"ss"*pi,-1*m);
            break;
          }
          if(summand.x() == spechtbasis::x2)
          {
            output.CreateAndAdd(homname+b,c,a,"ss"*pi,1*m);
            output.CreateAndAdd(homname+c,a,b,"s"*pi,-1*m);
            break;
          }
        }

        if(in.lambda() == specht::t3)
        {
          output.CreateAndAdd(homname+a,b,c,""*pi,1*m);
          output.CreateAndAdd(homname+b,c,a,"ss"*pi,1*m);
          output.CreateAndAdd(homname+c,a,b,"s"*pi,1*m);
          break;
        }
        if(in.lambda() == specht::s3)
        {
          std::cout<<"The element" << in.print() << "is incompatible" <<
            "with the homomorphism fkiit. Returning input."<<std::endl;
          return in;
        }
      }
    default:
      {
        std::cout<<"The element" << in.print() << "is incompatible" <<
          "with the homomorphism fkiit. Returning input."<<std::endl;
        return in;
      }
    }
  }
  std::cout << output.print()<<std::endl;
  return output;
}

CLinComb fiiir(CLinComb in)
{
  CLinComb output(subgroup::S3,specht::r);
  for(const auto& summand : in.summands())
  {
    std::string a = summand.tensor().first();
    std::string b = summand.tensor().second();
    std::string c = summand.tensor().third();
    CSymm pi = summand.coset();
    int m = summand.multiplicity();
    switch(in.type())
    {
      case subgroup::ts:
        {
          if(in.lambda() == specht::t2)
          {
            output.CreateAndAdd(homname+a,b,c,pi,2*m,spechtbasis::x1);
            output.CreateAndAdd(homname+a,b,c,pi,1*m,spechtbasis::x2);
            break;
          }
          if(in.lambda() == specht::s2)
          {
            output.CreateAndAdd(homname+a,b,c,pi,1*m,spechtbasis::x2);
            break;
          }
        }
      case subgroup::t:
        {
          auto in2 = isokii_iik(in);
          return fiiir(in2);
        }
      default:
        {
          std::cout<<"The element" << in.print() << "is incompatible" <<
            "with the homomorphism fiiir. Returning input."<<std::endl;
          return in;
        }
    }
  }
  std::cout << output.print()<<std::endl;
  return output;
}

CLinComb fiiit(CLinComb in)
{
  return changeFirstGen(in,subgroup::S3,specht::t3);
}
CLinComb fiiis(CLinComb in)
{
  return changeFirstGen(in,subgroup::S3,specht::s3);
}

// We only get an element of the form iik if the input
// is of type 1
CLinComb fiiks(CLinComb in)
{
  if(in.type() == subgroup::i)
  {
    return changeFirstGen(in,subgroup::t,specht::s2);
  }
  else
  {
    std::cout<<"The element" << in.print() << "is incompatible" <<
      "with the homomorphism fiikt. Returning input."<<std::endl;
    return in;
  }
}

CLinComb fiikt(CLinComb in)
{
  if(in.type() == subgroup::i)
  {
    return changeFirstGen(in,subgroup::t,specht::t2);
  }
  else
  {
    std::cout<<"The element" << in.print() << "is incompatible" <<
      "with the homomorphism fiikt. Returning input."<<std::endl;
    return in;
  }
}

CLinComb iso(CLinComb in, CSymm r)
{
  if(in.lambda() == specht::r)
  {
    // Handling the case with lambda=(2,1) would increase the
    // complexity of the function and we don't need this case.
    std::cout <<
      "The iso function is not designed to work with V(i,i,i;(2,1))."<<
      "Returning input unchanged." <<std::endl;
  }
  subgroup newtype = in.type();
  
  // type: if in.type() == subgroup::S3 or in.type() == subgroup::i the
  // type won't change. Otherwise we get the new type as
  // i.type() conjugate by r.
  if(in.type() != subgroup::S3 && in.type() != subgroup::i)
  {
    newtype = symmToSubgroup(inverse(r)*subgroupToSymm(in.type())*r);
  }
  CLinComb output(newtype,in.lambda());
  for(auto summand : in.summands())
  {
    CSymm pi = summand.coset();
    int m = summand.multiplicity();
    output.CreateAndAdd(summand.tensor()*r,inverse(r)*pi,m);
  }
  std::cout << output.print()<<std::endl;
  return output;
}

CLinComb isoiik_kii(CLinComb in)
{
  return iso(in,"tss");
}

CLinComb isoiki_kii(CLinComb in)
{
  return iso(in,"t");
}

CLinComb isokii_iki(CLinComb in)
{
  return iso(in,"t");
}

CLinComb isoiik_iki(CLinComb in)
{
  return iso(in,"ts");
}

CLinComb isokii_iik(CLinComb in)
{
  return iso(in,"tss");
}

CLinComb isoiki_iik(CLinComb in)
{
  return iso(in,"ts");
}


CLinComb isojik_ijk(CLinComb in)
{
  return iso(in,"t");
}

CLinComb isoikj_ijk(CLinComb in)
{
  return iso(in,"ts");
}

CLinComb isokji_ijk(CLinComb in)
{
  return iso(in,"tss");
}

CLinComb isojki_ijk(CLinComb in)
{
  return iso(in,"s");
}

CLinComb isokij_ijk(CLinComb in)
{
  return iso(in,"ss");
}
