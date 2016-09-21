#include "homomorphism.hpp"
#include "functions.hpp"
void checkHomomorphism(CLinComb (*hom)(CLinComb),CLinComb x);

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
  checkTwo(Hom(fikit)*Hom(isoiik_kii)*Hom(fkiit),
      Hom(fikit)*Hom(isoiik_kii)*Hom(fkiis),
      x_iiir1,
     "iii(2,1) -> i-ii(2) = iii-(2) -> i-ii-(2)",
     "iii(2,1) -> i-ii(1,1) = iii-(1,1) -> i-ii-(2)");
  checkTwo(Hom(fikit)*Hom(isoiik_kii)*Hom(fkiit),
      Hom(fikit)*Hom(isoiik_kii)*Hom(fkiis),
      x_iiir2,
     "iii(2,1) -> i-ii(2) = iii-(2) -> i-ii-(2)",
     "iii(2,1) -> i-ii(1,1) = iii-(1,1) -> i-ii-(2)");
  checkTwo(Hom(fikis)*Hom(isoiik_kii)*Hom(fkiit),
      Hom(fikis)*Hom(isoiik_kii)*Hom(fkiis),
      x_iiir1,
     "iii(2,1) -> i-ii(2) = iii-(2) -> i-ii-(1,1)",
     "iii(2,1) -> i-ii(1,1) = iii-(1,1) -> i-ii-(1,1)");
  checkTwo(Hom(fikis)*Hom(isoiik_kii)*Hom(fkiit),
      Hom(fikis)*Hom(isoiik_kii)*Hom(fkiit),
      x_iiir2,
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
