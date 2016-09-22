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
  CLinComb x_kiis(subgroup::ts,specht::s2);
  x_kiis.CreateAndAdd("a","b","c","",1);
  // a\otimes b\otimes c\in P(k,i,i;(2))
  CLinComb x_iikt(subgroup::t,specht::t2);
  x_iikt.CreateAndAdd("a","b","c","",1);
  // a\otimes b\otimes c\in P(i,i,k;(2))
  CLinComb x_iiks(subgroup::t,specht::s2);
  x_iiks.CreateAndAdd("a","b","c","",1);
  // a\otimes b\otimes c\in P(i,i,k;(1,1))
  CLinComb x_ikit(subgroup::tss,specht::t2);
  x_ikit.CreateAndAdd("a","b","c","",1);
  // a\otimes b\otimes c\in P(i,i,k;(2))
  CLinComb x_ikis(subgroup::tss,specht::s2);
  x_ikis.CreateAndAdd("a","b","c","",1);
  // a\otimes b\otimes c\in P(i,i,k;(1,1))
  CLinComb x_ijk(subgroup::i,specht::t1);
  x_ijk.CreateAndAdd("a","b","c","",1);
  // a\otimes b\otimes c\in P(i,j,k;(1))
  checkTwo(Hom(fikit)*Hom(isoiik_kii)*Hom(fkiit),
      Hom(fikit)*Hom(isoiik_kii)*Hom(fkiis),
      x_iiir1,
      "i","i","i","(2,1)","(2)","(1,1)","(2)");
  checkTwo(Hom(fikit)*Hom(isoiik_kii)*Hom(fkiit),
      Hom(fikit)*Hom(isoiik_kii)*Hom(fkiis),
      x_iiir2,
      "i","i","i","(2,1)","(2)","(1,1)","(2)");
  checkTwo(Hom(fikis)*Hom(isoiik_kii)*Hom(fkiit),
      Hom(fikis)*Hom(isoiik_kii)*Hom(fkiis),
      x_iiir1,
      "i","i","i","(2,1)","(2)","(1,1)","(1,1)");
  checkTwo(Hom(fikis)*Hom(isoiik_kii)*Hom(fkiit),
      Hom(fikis)*Hom(isoiik_kii)*Hom(fkiit),
      x_iiir2,
      "i","i","i","(2,1)","(2)","(1,1)","(1,1)");
  //those are all possibilities for P(i,i,i;lambda) as there
  //is only one non-zero descending path from P(i,i,i;lambda)
  //if lambda is (3) or (1,1,1)
  checkThree(Hom(isojik_ijk)*Hom(fiikt)*Hom(isojik_ijk)*Hom(fijk),
             Hom(fkiit)*Hom(fiiit)*Hom(isokii_iki),
             Hom(fkiit)*Hom(fiiir)*Hom(isokii_iki),
             x_ikit,
             "i","i*","i","(2)","","(3)","(2,1)","(2)",comptype::ijj);
  checkTwo(Hom(isojik_ijk)*Hom(fiiks)*Hom(isojik_ijk)*Hom(fijk),
             Hom(fkiis)*Hom(fiiir)*Hom(isokii_iki),
             x_ikit,
             "i","i*","i","(2)","","(2,1)","(1,1)");

  checkThree(Hom(isojik_ijk)*Hom(fiiks)*Hom(isojik_ijk)*Hom(fijk),
             Hom(fkiis)*Hom(fiiis)*Hom(isokii_iki),
             Hom(fkiis)*Hom(fiiir)*Hom(isokii_iki),
             x_ikis,
             "i","i*","i","(1,1)","","(1,1,1)","(2,1)","(1,1)",comptype::ijj);
  checkTwo(Hom(isojik_ijk)*Hom(fiiks)*Hom(isojik_ijk)*Hom(fijk),
             Hom(fkiis)*Hom(fiiir)*Hom(isokii_iki),
             x_ikis,
             "i","i*","i","(1,1)","","(2,1)","(2)");
  // note that there is only one non-zero path iii+mu -> i-i-i+mu
  // as the intermediate module is always i-ii+
  checkThree(Hom(isojik_ijk)*Hom(fijk)*Hom(isoiki_kii)*Hom(fkiit),
             Hom(fijk)*Hom(fiikt)*Hom(isoiki_kii),
             Hom(fijk)*Hom(fiiks)*Hom(isoiki_kii),
             x_kiit,
             "i","i*","i*","(2)","(2)","(2)","(1,1)","",comptype::ijj);
  checkTwo( Hom(fiiir)*Hom(isokii_iki)*Hom(fikit),
            Hom(fiiir)*Hom(isokii_iki)*Hom(fikis),
             x_iikt,
             "i*","i*","i","(2)","(2)","(1,1)","(2,1)");
  // There is only one descending non-zero path ii+i+mu -> iiilambda
  checkTwo( Hom(isojik_ijk)*Hom(fijk)*Hom(isoiik_kii)*Hom(fkiit),
            Hom(isokji_ijk)*Hom(fijk)*Hom(isokji_ijk)*Hom(fijk)*Hom(isoiki_kii),
             x_kiit,
             "j","i","i","(2)","(2)","","");
  // again there is only one non-zero path jiimu -> ji-i-mu
  // as the intermediate module is always ji-i
  checkThree(Hom(isojik_ijk)*Hom(fijk)*Hom(isoiki_kii)*Hom(fkiis),
             Hom(fijk)*Hom(fiiks)*Hom(isoiki_kii),
             Hom(fijk)*Hom(fiikt)*Hom(isoiki_kii),
             x_kiis,
             "i","i*","i*","(1,1)","(1,1)","(1,1)","(2)","",comptype::ijj);
  checkTwo( Hom(fiiir)*Hom(isokii_iki)*Hom(fikis),
            Hom(fiiir)*Hom(isokii_iki)*Hom(fikit),
             x_iiks,
             "i*","i*","i","(1,1)","(1,1)","(2)","(2,1)");
  // There is only one descending non-zero path ii+i+mu -> iiilambda
  checkTwo( Hom(isojik_ijk)*Hom(fijk)*Hom(isoiik_kii)*Hom(fkiit),
            Hom(fijk)*Hom(isokji_ijk)*Hom(fijk)*Hom(isoiki_kii),
             x_kiit,
             "j","i","i","(2)","(2)","","");
  // again there is only one non-zero path jiimu -> ji-i-mu
  // as the intermediate module is always ji-i



  checkThree(Hom(isoiki_kii)*Hom(fkiit)*Hom(isokii_iki)*Hom(fikit),
             Hom(fikit)*Hom(fiikt)*Hom(isojik_ijk),
             Hom(fikit)*Hom(fiiks)*Hom(isojik_ijk),
             x_ijk,
             "i","i*","i'","","(2)","(2)","(1,1)","(2)",comptype::ijj);

  checkThree(Hom(isoiki_kii)*Hom(fkiis)*Hom(isokii_iki)*Hom(fikis),
             Hom(fikis)*Hom(fiiks)*Hom(isojik_ijk),
             Hom(fikis)*Hom(fiikt)*Hom(isojik_ijk),
             x_ijk,
             "i","i*","i'","","(1,1)","(1,1)","(2)","(1,1)",comptype::ijj);

  checkTwo(Hom(isokii_iki)*Hom(fikit)*Hom(isojik_ijk)*Hom(fijk),
             Hom(fkiit)*Hom(isokii_iki)*Hom(fikit)*Hom(isojik_ijk),
             x_ijk,
             "j","i*","i","","","(2)","(2)");
  checkTwo(Hom(isokii_iki)*Hom(fikis)*Hom(isojik_ijk)*Hom(fijk),
             Hom(fkiis)*Hom(isokii_iki)*Hom(fikis)*Hom(isojik_ijk),
             x_ijk,
             "j","i*","i","","","(1,1)","(1,1)");

  checkThree(Hom(isojik_ijk)*Hom(fijk)*Hom(isojik_ijk)*Hom(fijk),
             Hom(fijk)*Hom(fiikt)*Hom(isojik_ijk),
             Hom(fijk)*Hom(fiiks)*Hom(isojik_ijk),
             x_ijk,
             "i","i*","j","","","(2)","(1,1)","",comptype::ijj);

  checkTwo(Hom(isojik_ijk)*Hom(fijk)*Hom(isojik_ijk)*Hom(fijk),
             Hom(fijk)*Hom(isojik_ijk)*Hom(fijk)*Hom(isojik_ijk),
             x_ijk,
             "i","j","k","","","","");
}
