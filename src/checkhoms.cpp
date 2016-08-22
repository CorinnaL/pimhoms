#include "homomorphism.hpp"
#include "CLinComb.hpp"

void checkHomomorphism(CLinComb (*hom)(CLinComb),CLinComb x)
{
  auto hx = hom(x);
  auto xs = x*"s";
  auto xt = x*"t";
  std::cout << "x=" << x.print() << std::endl;
  std::cout << "x*s=" << xs.print() << std::endl;
  std::cout << "x*t=" << xt.print() << std::endl;
  std::cout << "f(x)=" << hx.print() << std::endl;
  auto hxs = hx*"s";
  std::cout << "f(x*s)=" << hom(xs).print() << std::endl;
  std::cout << "f(x)*s=" << hxs.print() << std::endl;
  auto hxt = hx*"t";
  std::cout << "f(x*t)=" << hom(xt).print() << std::endl;
  std::cout << "f(x)*t=" << hxt.print() << std::endl;
  if(hxs.print() == hom(xs).print() &&
     hxt.print() == hom(xt).print())
  {
    std::cout<<"The tested map is a homomorphism on the tested element."<<std::endl;
  }
  else
  {
    std::cout<<"The tested map is NOT a homomorphism."<<std::endl;
  }
}

void checkAll()
{
  auto s = CSymm("ss");
  auto t = factorize(s);
  for(auto x : t)
  {
    std::cout << x.get() <<std::endl;
  }
  CLinComb x_iiir1(subgroup::S3,specht::r);
  x_iiir1.CreateAndAdd("a","b","c","",1,spechtbasis::x1);
  // a\otimes b\otimes c\otimes x_1\in P(i,i,i;(2,1))
  CLinComb x_iiir2(subgroup::S3,specht::r);
  x_iiir2.CreateAndAdd("a","b","c","",1,spechtbasis::x2);
  // a\otimes b\otimes c\otimes x_2\in P(i,i,i;(2,1))
  CLinComb x_iiit(subgroup::S3,specht::t3);
  x_iiit.CreateAndAdd("a","b","c","",1);
  // a\otimes b\otimes c\in P(i,i,i;(3))
  CLinComb x_iiis(subgroup::S3,specht::s3);
  x_iiis.CreateAndAdd("a","b","c","",1);
  // a\otimes b\otimes c\in P(i,i,i;(1,1,1))
  CLinComb x_iikt(subgroup::t,specht::t2);
  x_iikt.CreateAndAdd("a","b","c","",1);
  // a\otimes b\otimes c\in P(i,i,k;(2))
  CLinComb x_iiks(subgroup::t,specht::s2);
  x_iiks.CreateAndAdd("a","b","c","",1);
  // a\otimes b\otimes c\in P(i,i,k;(1,1))
  CLinComb x_ijk(subgroup::i,specht::t1);
  x_ijk.CreateAndAdd("a","b","c","",1);
  // a\otimes b\otimes c\in P(i,j,k;(1))
  std::cout<< "---------------------------------------------------" <<std::endl;
  std::cout<< "P(i,i,i;(2,1)" <<std::endl;
  std::cout<< "---------------------------------------------------" <<std::endl;
  checkHomomorphism(fkiit,x_iiir1);
  checkHomomorphism(fkiit,x_iiir2);
  checkHomomorphism(fkiis,x_iiir1);
  checkHomomorphism(fkiis,x_iiir2);
  std::cout<< "---------------------------------------------------" <<std::endl;
  std::cout<< "P(i,i,k;(2)" <<std::endl;
  std::cout<< "---------------------------------------------------" <<std::endl;
  checkHomomorphism(fijk,x_iikt);
  checkHomomorphism(fikit,x_iikt);
  checkHomomorphism(fikis,x_iikt);
  checkHomomorphism(fkiit,isokii_iik(x_iikt));
  checkHomomorphism(fiiir,isokii_iik(x_iikt));
  std::cout<< "---------------------------------------------------" <<std::endl;
  std::cout<< "P(i,i,k;(1,1)" <<std::endl;
  std::cout<< "---------------------------------------------------" <<std::endl;
  checkHomomorphism(fijk,x_iiks);
  checkHomomorphism(fikit,x_iiks);
  checkHomomorphism(fikis,x_iiks);
  checkHomomorphism(fkiis,isokii_iik(x_iiks));
  checkHomomorphism(fiiir,isokii_iik(x_iiks));
  std::cout<< "---------------------------------------------------" <<std::endl;
  std::cout<< "P(i,j,k;(1)" <<std::endl;
  std::cout<< "---------------------------------------------------" <<std::endl;
  checkHomomorphism(fijk,x_ijk);
  checkHomomorphism(fiikt,x_ijk);
  checkHomomorphism(fiiks,x_ijk);
}
