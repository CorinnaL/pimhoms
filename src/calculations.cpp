#include "homomorphism.hpp"
#include "CLinComb.hpp"
#include <functional>
void checkHomomorphism(CLinComb (*hom)(CLinComb),CLinComb x);
using Hom = std::function<CLinComb (CLinComb)>;
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

void checkTwo(Hom hom1,
      Hom  hom2,
      CLinComb x,std::string firstcomp, std::string secondcomp)
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
void checkTwo(std::function<CLinComb (CLinComb)> hom11,std::function<CLinComb (CLinComb)> hom12,
              CLinComb (*hom13)(CLinComb),
              CLinComb (*hom21)(CLinComb),CLinComb (*hom22)(CLinComb),
              CLinComb (*hom23)(CLinComb),
              CLinComb x,std::string firstcomp, std::string secondcomp)
{
  std::cout << x.print() << std::endl;
  std::cout << "***************"<<firstcomp<<"***************" <<std::endl;
  CLinComb c1 = hom13(hom12(hom11(x)));
  std::cout << "***************"<<secondcomp<<"***************" <<std::endl;
  CLinComb c2 = hom23(hom22(hom21(x)));
  if(c1.print() != c2.print())
  std::cout << "!!!!!!!Check!!!!!!!!!"<<std::endl;
  std::cout << "------------------------------------------------------------"<<std::endl;
}
//{
//  auto hx = hom(x);
//  auto xs = x*"s";
//  auto xt = x*"t";
//  std::cout << "x=" << x.print() << std::endl;
//  std::cout << "x*s=" << xs.print() << std::endl;
//  std::cout << "x*t=" << xt.print() << std::endl;
//  std::cout << "f(x)=" << hx.print() << std::endl;
//  auto hxs = hx*"s";
//  std::cout << "f(x*s)=" << hom(xs).print() << std::endl;
//  std::cout << "f(x)*s=" << hxs.print() << std::endl;
//  auto hxt = hx*"t";
//  std::cout << "f(x*t)=" << hom(xt).print() << std::endl;
//  std::cout << "f(x)*t=" << hxt.print() << std::endl;
//  if(hxs.print() == hom(xs).print() &&
//     hxt.print() == hom(xt).print())
//  {
//    std::cout<<"The tested map is a homomorphism on the tested element."<<std::endl;
//  }
//  else
//  {
//    std::cout<<"The tested map is NOT a homomorphism."<<std::endl;
//  }
//}

void calculateStuff()
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
  CLinComb x_kiit(subgroup::ts,specht::t2);
  x_kiit.CreateAndAdd("a","b","c","",1);
  // a\otimes b\otimes c\in P(k,i,i;(2))
  CLinComb x_iikt(subgroup::t,specht::t2);
  x_iikt.CreateAndAdd("a","b","c","",1);
  // a\otimes b\otimes c\in P(i,i,k;(2))
  CLinComb x_iiks(subgroup::t,specht::s2);
  x_iiks.CreateAndAdd("a","b","c","",1);
  // a\otimes b\otimes c\in P(i,i,k;(1,1))
  CLinComb x_ijk(subgroup::i,specht::t1);
  x_ijk.CreateAndAdd("a","b","c","",1);
  // a\otimes b\otimes c\in P(i,j,k;(1))
  checkTwo(fkiit,isoiik_kii,fikit,fkiis,isoiik_kii,fikit,x_iiir1,
     "iii(2,1) -> i-ii(2) = iii-(2) -> i-ii-(2)",
     "iii(2,1) -> i-ii(1,1) = iii-(1,1) -> i-ii-(2)");
  checkTwo(fkiit,isoiik_kii,fikit,fkiis,isoiik_kii,fikit,x_iiir2,
     "iii(2,1) -> i-ii(2) = iii-(2) -> i-ii-(2)",
     "iii(2,1) -> i-ii(1,1) = iii-(1,1) -> i-ii-(2)");
  checkTwo(fkiit,isoiik_kii,fikis,fkiis,isoiik_kii,fikis,x_iiir1,
     "iii(2,1) -> i-ii(2) = iii-(2) -> i-ii-(1,1)",
     "iii(2,1) -> i-ii(1,1) = iii-(1,1) -> i-ii-(1,1)");
  checkTwo(fkiit,isoiik_kii,fikis,fkiis,isoiik_kii,fikis,x_iiir2,
     "iii(2,1) -> i-ii(2) = iii-(2) -> i-ii-(1,1)",
     "iii(2,1) -> i-ii(1,1) = iii-(1,1) -> i-ii-(1,1)");
  //those are all possibilities for P(i,i,i;lambda) as there
  //is only one non-zero descending path from P(i,i,i;lambda)
  //if lambda is (3) or (1,1,1)
  checkThree(Hom(fiikt)*Hom(isokji_ijk)*Hom(fijk),
             Hom(fkiit)*Hom(fiiit)*Hom(isokii_iik),
             Hom(fkiit)*Hom(fiiir)*Hom(isokii_iik),
             x_iikt,
     "iii+(2) -> i-ii+ = i+ii- -> iii-(2)",
     "iii+(2) = i+ii(2) -> iii(3) -> i-ii(2)",
     "iii+(2) = i+ii(2) -> iii(2,1) -> i-ii(2)");
  checkTwo(Hom(fiiks)*Hom(isokji_ijk)*Hom(fijk),
             Hom(fkiis)*Hom(fiiir)*Hom(isokii_iik),
             x_iikt,
     "iii+(2) -> i-ii+ = i+ii- -> iii-(1,1)",
     "iii+(2) = i+ii(2) -> iii(2,1) -> i-ii(1,1)");
  checkThree(Hom(isokii_iik)*Hom(fiiks)*Hom(isokji_ijk)*Hom(fijk),
             Hom(fkiis)*Hom(fiiis)*Hom(isokii_iik),
             Hom(fkiis)*Hom(fiiir)*Hom(isokii_iik),
             x_iiks,
     "iii+(1,1) -> i-ii+ = i+ii- -> iii-(1,1) = i-ii(1,1)",
     "iii+(1,1) = i+ii(1,1) -> iii(1,1,1) -> i-ii(1,1)",
     "iii+(1,1) = i+ii(1,1) -> iii(2,1) -> i-ii(1,1)");
  checkTwo(Hom(isokii_iik)*Hom(fiikt)*Hom(isokji_ijk)*Hom(fijk),
             Hom(fkiit)*Hom(fiiir)*Hom(isokii_iik),
             x_iiks,
     "iii+(1,1) -> i-ii+ = i+ii- -> iii-(2) = i-ii(2)",
     "iii+(1,1) = i+ii(1,1) -> iii(2,1) -> i-ii(2)");
  // note that there is only one non-zero path iii+mu -> i-i-i+mu
  // as the intermediate module is always i-ii+
  checkThree(Hom(isokij_ijk)*Hom(fijk)*Hom(isoiik_kii)*Hom(fkiit),
             Hom(fijk)*Hom(isoiik_iki)*Hom(fikit)*Hom(isoiik_kii),
             Hom(fijk)*Hom(isoiik_iki)*Hom(fikis)*Hom(isoiik_kii),
             x_kiit,
     "ii+i+(2) -> i-i+i+(2) = i+i+i-(2) -> ii+i- = i-ii+",
     "ii+i+(2) = i+i+i(2) -> ii+i(2) = iii+(2) ->i-ii+",
     "ii+i+(2) = i+i+i(2) -> ii+i(1,1) = iii+(2) ->i-ii+");
  checkTwo( Hom(fiiir)*Hom(isokii_iki)*Hom(fikit)*Hom(isoiik_kii),
            Hom(fiiir)*Hom(isokii_iki)*Hom(fikis)*Hom(isoiik_kii),
             x_kiit,
     "ii+i+(2) = i+i+i(2) -> ii+i(2) = i+ii(2) -> iii(2,1)",
     "ii+i+(2) = i+i+i(2) -> ii+i(1,1) = i+ii(1,1) -> iii(2,1)");
  // There is only one descending non-zero path ii+i+mu -> iiilambda
  checkTwo( Hom(fijk)*Hom(isoiik_kii)*Hom(fkiit),
            Hom(isokji_ijk)*Hom(fijk)*Hom(isokji_ijk)*Hom(fijk)*Hom(isoiik_kii),
             x_kiit,
     "jii(2) -> j-ii(2) = iij-(2) -> i-ij",
     "jii(2) = iij(2) -> i-ij = jii- -> j-ii- -> i-ij");
  // again there is only one non-zero path jiimu -> ji-i-mu
  // as the intermediate module is always ji-i


  checkThree(Hom(fkiit)*Hom(isokii_iik)*Hom(fiikt)*Hom(isojik_ijk),
             Hom(isokii_iki)*Hom(fikit)*Hom(fiikt)*Hom(isokji_ijk),
             Hom(isokii_iki)*Hom(fikit)*Hom(fiiks)*Hom(isokji_ijk),
             x_ijk,
     "ii+i+2 = i+ii+2 -> iii+2(2) = i+2ii(2) -> i+ii(2)",
     "ii+i+2 = i+2i+i -> i+i+i(2) -> ii+i(2) =  i+ii(2)",
     "ii+i+2 = i+2i+i -> i+i+i(1,1) -> ii+i(2) =  i+ii(2)");

  checkThree(Hom(fkiis)*Hom(isokii_iik)*Hom(fiiks)*Hom(isojik_ijk),
             Hom(isokii_iki)*Hom(fikis)*Hom(fiikt)*Hom(isokji_ijk),
             Hom(isokii_iki)*Hom(fikis)*Hom(fiiks)*Hom(isokji_ijk),
             x_ijk,
     "ii+i+2 = i+ii+2 -> iii+2(1,1) = i+2ii(1,1) -> i+ii(1,1)",
     "ii+i+2 = i+2i+i -> i+i+i(2) -> ii+i(1,1) =  i+ii(1,1)",
     "ii+i+2 = i+2i+i -> i+i+i(1,1) -> ii+i(1,1) =  i+ii(1,1)");

  //TODO: this looks generic.
  checkTwo(Hom(isojik_ijk)*Hom(fijk)*Hom(isojik_ijk)*Hom(fijk),
             Hom(fijk)*Hom(isojik_ijk)*Hom(fijk)*Hom(isojik_ijk),
             x_ijk,
     "jii+ -> j-ii+ = ij-i+ -> i-j-i+ -> j-i-i+",
     "jii+ = iji+ -> i-ji+ = ji-i+ ->j-i-i+");
  checkTwo(Hom(isokii_iki)*Hom(fikit)*Hom(isokij_ijk)*Hom(fijk),
             Hom(fkiit)*Hom(isokii_iki)*Hom(fikit)*Hom(isokij_ijk),
             x_ijk,
     "jii+ -> j-ii+ = i+j-i -> ij-i(2) -> j-ii(2)",
     "jii+ = i+ji -> iji(2) = jii(2) ->j-ii(2)");
  checkTwo(Hom(isokii_iki)*Hom(fikis)*Hom(isokij_ijk)*Hom(fijk),
             Hom(fkiis)*Hom(isokii_iki)*Hom(fikis)*Hom(isokij_ijk),
             x_ijk,
     "jii+ -> j-ii+ = i+j-i -> ij-i(1,1) -> j-ii(1,1)",
     "jii+ = i+ji -> iji(1,1) = jii(1,1) ->j-ii(1,1)");
}
