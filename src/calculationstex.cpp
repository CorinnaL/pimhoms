#include "homomorphism.hpp"
#include "functions.hpp"
void checkHomomorphism(CLinComb (*hom)(CLinComb),CLinComb x);
void calculateStufftex()
{
    std::cout<<"\\begin{align*}\\\\"<<std::endl;
  auto s = CSymm("ss");
  auto t = factorize(s);
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
  // a\otimes b\otimes c\in P(k,i,i;(1,1))
  CLinComb x_iikt(subgroup::t,specht::t2);
  x_iikt.CreateAndAdd("a","b","c","",1);
  // a\otimes b\otimes c\in P(i,i,k;(2))
  CLinComb x_iiks(subgroup::t,specht::s2);
  x_iiks.CreateAndAdd("a","b","c","",1);
  // a\otimes b\otimes c\in P(i,i,k;(1,1))
  CLinComb x_ijk(subgroup::i,specht::t1);
  x_ijk.CreateAndAdd("a","b","c","",1);
  // a\otimes b\otimes c\in P(i,j,k;(1))



  //Idea: i' is either i+1 or i-1 which means
  //we cover both the ascending and the
  //descending paths.
  checkTwo(Hom(fikit)*Hom(isoiik_kii)*Hom(fkiit),
      Hom(fikit)*Hom(isoiik_kii)*Hom(fkiis),
      x_iiir1,
     "\\ol{P'}(i,i,i;(2,1)) &\\rightarrow \\ol{P'}(i',i,i;(2)) \\cong \\ol{P'}(i,i,i';(2)) \\rightarrow \\ol{P'}(i',i,i';(2))",
     "\\ol{P'}(i,i,i;(2,1)) &\\rightarrow \\ol{P'}(i',i,i;(1,1)) \\cong \\ol{P'}(i,i,i';(1,1)) \\rightarrow \\ol{P'}(i',i,i';(2))");
  checkTwo(Hom(fikit)*Hom(isoiik_kii)*Hom(fkiit),
      Hom(fikit)*Hom(isoiik_kii)*Hom(fkiis),
      x_iiir2,
     "\\ol{P'}(i,i,i;(2,1)) &\\rightarrow \\ol{P'}(i',i,i;(2)) \\cong \\ol{P'}(i,i,i';(2)) \\rightarrow \\ol{P'}(i',i,i';(2))",
     "\\ol{P'}(i,i,i;(2,1)) &\\rightarrow \\ol{P'}(i',i,i;(1,1)) \\cong \\ol{P'}(i,i,i';(1,1)) \\rightarrow \\ol{P'}(i',i,i';(2))");
  checkTwo(Hom(fikis)*Hom(isoiik_kii)*Hom(fkiit),
      Hom(fikis)*Hom(isoiik_kii)*Hom(fkiis),
      x_iiir1,
     "\\ol{P'}(i,i,i;(2,1)) &\\rightarrow \\ol{P'}(i',i,i;(2)) \\cong \\ol{P'}(i,i,i';(2)) \\rightarrow \\ol{P'}(i',i,i';(1,1))",
     "\\ol{P'}(i,i,i;(2,1)) &\\rightarrow \\ol{P'}(i',i,i;(1,1)) \\cong \\ol{P'}(i,i,i';(1,1)) \\rightarrow \\ol{P'}(i',i,i';(1,1))");
  checkTwo(Hom(fikis)*Hom(isoiik_kii)*Hom(fkiit),
      Hom(fikis)*Hom(isoiik_kii)*Hom(fkiit),
      x_iiir2,
     "\\ol{P'}(i,i,i;(2,1)) &\\rightarrow \\ol{P'}(i',i,i;(2)) \\cong \\ol{P'}(i,i,i';(2)) \\rightarrow \\ol{P'}(i',i,i';(1,1))",
     "\\ol{P'}(i,i,i;(2,1)) &\\rightarrow \\ol{P'}(i',i,i;(1,1)) \\cong \\ol{P'}(i,i,i';(1,1)) \\rightarrow \\ol{P'}(i',i,i';(1,1))");
  //those are all possibilities for P(i,i,i;lambda) as there
  //is only one non-zero descending path from P(i,i,i;lambda)
  //if lambda is (3) or (1,1,1)

  //i* = i+1 if i' = i-1 and i*=i-1 if i'=i+1, so i*'=i
  checkThree(Hom(isokii_iik)*Hom(fiikt)*Hom(isokji_ijk)*Hom(fijk),
             Hom(fkiit)*Hom(fiiit)*Hom(isokii_iik),
             Hom(fkiit)*Hom(fiiir)*Hom(isokii_iik),
             x_iikt,
     "\\ol{P'}(i,i,i^*;(2)) &\\rightarrow \\ol{P'}(i',i,i^*) \\cong \\ol{P'}(i^*,i,i') \\rightarrow \\ol{P'}(i,i,i';(2))\\cong\\ol{P'}(i',i,i;(2))",
     "\\ol{P'}(i,i,i^*;(2)) &\\cong \\ol{P'}(i^*,i,i;(2)) \\rightarrow \\ol{P'}(i,i,i;(3)) \\rightarrow \\ol{P'}(i',i,i;(2))",
     "\\ol{P'}(i,i,i^*;(2)) &\\cong \\ol{P'}(i^*,i,i;(2)) \\rightarrow \\ol{P'}(i,i,i;(2,1)) \\rightarrow \\ol{P'}(i',i,i;(2))");
  checkTwo(Hom(isokii_iik)*Hom(fiiks)*Hom(isokji_ijk)*Hom(fijk),
             Hom(fkiis)*Hom(fiiir)*Hom(isokii_iik),
             x_iikt,
     "\\ol{P'}(i,i,i^*;(2)) &\\rightarrow \\ol{P'}(i',i,i^*) \\cong \\ol{P'}(i^*,i,i') \\rightarrow \\ol{P'}(i,i,i';(1,1))\\cong\\ol{P'}(i',i,i;(1,1))",
     "\\ol{P'}(i,i,i^*;(2)) &\\cong \\ol{P'}(i^*,i,i;(2)) \\rightarrow \\ol{P'}(i,i,i;(2,1)) \\rightarrow \\ol{P'}(i',i,i;(1,1))");
  checkThree(Hom(isokii_iik)*Hom(fiiks)*Hom(isokji_ijk)*Hom(fijk),
             Hom(fkiis)*Hom(fiiis)*Hom(isokii_iik),
             Hom(fkiis)*Hom(fiiir)*Hom(isokii_iik),
             x_iiks,
     "\\ol{P'}(i,i,i^*;(1,1)) &\\rightarrow \\ol{P'}(i',i,i^*) \\cong \\ol{P'}(i^*,i,i') \\rightarrow \\ol{P'}(i,i,i';(1,1)) \\cong \\ol{P'}(i',i,i;(1,1))",
     "\\ol{P'}(i,i,i^*;(1,1)) &\\cong \\ol{P'}(i^*,i,i;(1,1)) \\rightarrow \\ol{P'}(i,i,i;(1,1,1)) \\rightarrow \\ol{P'}(i',i,i;(1,1))",
     "\\ol{P'}(i,i,i^*;(1,1)) &\\cong \\ol{P'}(i^*,i,i;(1,1)) \\rightarrow \\ol{P'}(i,i,i;(2,1)) \\rightarrow \\ol{P'}(i',i,i;(1,1))");
  checkTwo(Hom(isokii_iik)*Hom(fiikt)*Hom(isokji_ijk)*Hom(fijk),
             Hom(fkiit)*Hom(fiiir)*Hom(isokii_iik),
             x_iiks,
     "\\ol{P'}(i,i,i^*;(1,1)) &\\rightarrow \\ol{P'}(i',i,i^*) \\cong \\ol{P'}(i^*,i,i') \\rightarrow \\ol{P'}(i,i,i';(2)) \\cong \\ol{P'}(i',i,i;(2))",
     "\\ol{P'}(i,i,i^*;(1,1)) &\\cong \\ol{P'}(i^*,i,i;(1,1)) \\rightarrow \\ol{P'}(i,i,i;(2,1)) -> \\ol{P'}(i',i,i;(2))");
  // note that there is only one non-zero path iii*mu -> i'i'i*mu
  // as the intermediate module is always \\ol{P'}(i',i,i^*)
  checkThree(Hom(isokij_ijk)*Hom(fijk)*Hom(isoiik_kii)*Hom(fkiit),
             Hom(fijk)*Hom(isoiik_iki)*Hom(fikit)*Hom(isoiik_kii),
             Hom(fijk)*Hom(isoiik_iki)*Hom(fikis)*Hom(isoiik_kii),
             x_kiit,
     "\\ol{P'}(i,i^*,i^*;(2)) &\\rightarrow \\ol{P'}(i',i^*,i^*;(2)) \\cong \\ol{P'}(i^*,i^*,i';(2)) \\rightarrow \\ol{P'}(i,i^*,i') \\cong \\ol{P'}(i',i,i^*)",
     "\\ol{P'}(i,i^*,i^*;(2)) &\\cong \\ol{P'}(i^*,i^*,i;(2)) \\rightarrow \\ol{P'}(i,i^*,i;(2)) \\cong \\ol{P'}(i,i,i^*;(2)) \\rightarrow\\ol{P'}(i',i,i^*)",
     "\\ol{P'}(i,i^*,i^*;(2)) &\\cong \\ol{P'}(i^*,i^*,i;(2)) \\rightarrow \\ol{P'}(i,i^*,i;(1,1)) \\cong \\ol{P'}(i,i,i^*;(1,1)) \\rightarrow\\ol{P'}(i',i,i^*)");
  checkTwo( Hom(fiiir)*Hom(isokii_iki)*Hom(fikit)*Hom(isoiik_kii),
            Hom(fiiir)*Hom(isokii_iki)*Hom(fikis)*Hom(isoiik_kii),
             x_kiit,
     "\\ol{P'}(i,i^*,i^*;(2)) &\\cong \\ol{P'}(i^*,i^*,i;(2)) \\rightarrow \\ol{P'}(i,i^*,i;(2)) \\cong \\ol{P'}(i^*,i,i;(2)) \\rightarrow \\ol{P'}(i,i,i;(2,1))",
     "\\ol{P'}(i,i^*,i^*;(2)) &\\cong \\ol{P'}(i^*,i^*,i;(2)) \\rightarrow \\ol{P'}(i,i^*,i;(1,1)) \\cong \\ol{P'}(i^*,i,i;(1,1)) \\rightarrow \\ol{P'}(i,i,i;(2,1))");
  // There is only one descending non-zero path ii*i*mu -> iiilambda
  checkTwo( Hom(fijk)*Hom(isoiik_kii)*Hom(fkiit),
            Hom(isokji_ijk)*Hom(fijk)*Hom(isokji_ijk)*Hom(fijk)*Hom(isoiik_kii),
             x_kiit,
     "\\ol{P}(j,i,i;(2)) &\\rightarrow \\ol{P'}(j',i,i;(2)) \\cong \\ol{P'}(i,i,j';(2)) \\rightarrow \\ol{P'}(i',i,j')",
     "\\ol{P'}(j,i,i;(2)) &\\cong \\ol{P'}(i,i,j;(2)) \\rightarrow \\ol{P'}(i',i,j) \\cong \\ol{P'}(j,i,i') \\rightarrow \\ol{P'}(j',i,i') \\rightarrow \\ol{P'}(i',i,j')");
  // again there is only one non-zero path jiimu -> ji'i'mu
  // as the intermediate module is always ji'i

  checkThree(Hom(isokij_ijk)*Hom(fijk)*Hom(isoiik_kii)*Hom(fkiis),
             Hom(fijk)*Hom(isoiik_iki)*Hom(fikit)*Hom(isoiik_kii),
             Hom(fijk)*Hom(isoiik_iki)*Hom(fikis)*Hom(isoiik_kii),
             x_kiis,
     "\\ol{P'}(i,i^*,i^*;(1,1)) &\\rightarrow \\ol{P'}(i',i^*,i^*;(1,1)) \\cong \\ol{P'}(i^*,i^*,i';(1,1)) \\rightarrow \\ol{P'}(i,i^*,i') \\cong \\ol{P'}(i',i,i^*)",
     "\\ol{P'}(i,i^*,i^*;(1,1)) &\\cong \\ol{P'}(i^*,i^*,i;(1,1)) \\rightarrow \\ol{P'}(i,i^*,i;(2)) \\cong \\ol{P'}(i,i,i^*;(2)) \\rightarrow\\ol{P'}(i',i,i^*)",
     "\\ol{P'}(i,i^*,i^*;(1,1)) &\\cong \\ol{P'}(i^*,i^*,i;(1,1)) \\rightarrow \\ol{P'}(i,i^*,i;(1,1)) \\cong \\ol{P'}(i,i,i^*;(1,1)) \\rightarrow\\ol{P'}(i',i,i^*)");
  checkTwo( Hom(fiiir)*Hom(isokii_iki)*Hom(fikit)*Hom(isoiik_kii),
            Hom(fiiir)*Hom(isokii_iki)*Hom(fikis)*Hom(isoiik_kii),
             x_kiis,
     "\\ol{P'}(i,i^*,i^*;(1,1)) &\\cong \\ol{P'}(i^*,i^*,i;(1,1)) \\rightarrow \\ol{P'}(i,i^*,i;(2)) \\cong \\ol{P'}(i^*,i,i;(2)) \\rightarrow \\ol{P'}(i,i,i;(2,1))",
     "\\ol{P'}(i,i^*,i^*;(1,1)) &\\cong \\ol{P'}(i^*,i^*,i;(1,1)) \\rightarrow \\ol{P'}(i,i^*,i;(1,1)) \\cong \\ol{P'}(i^*,i,i;(1,1)) \\rightarrow \\ol{P'}(i,i,i;(2,1))");
  // There is only one descending non-zero path ii*i*mu -> iiilambda
  checkTwo( Hom(fijk)*Hom(isoiik_kii)*Hom(fkiis),
            Hom(isokji_ijk)*Hom(fijk)*Hom(isokji_ijk)*Hom(fijk)*Hom(isoiik_kii),
             x_kiis,
     "\\ol{P}(j,i,i;(1,1)) &\\rightarrow \\ol{P'}(j',i,i;(1,1)) \\cong \\ol{P'}(i,i,j';(1,1)) \\rightarrow \\ol{P'}(i',i,j')",
     "\\ol{P'}(j,i,i;(1,1)) &\\cong \\ol{P'}(i,i,j;(1,1)) \\rightarrow \\ol{P'}(i',i,j) \\cong \\ol{P'}(j,i,i') \\rightarrow \\ol{P'}(j',i,i') \\rightarrow \\ol{P'}(i',i,j')");
  // again there is only one non-zero path jiimu -> ji'i'mu
  // as the intermediate module is always ji'i



  checkThree(Hom(fkiit)*Hom(isokii_iik)*Hom(fiikt)*Hom(isojik_ijk),
             Hom(isokii_iki)*Hom(fikit)*Hom(fiikt)*Hom(isokji_ijk),
             Hom(isokii_iki)*Hom(fikit)*Hom(fiiks)*Hom(isokji_ijk),
             x_ijk,
     "\\ol{P'}(i',i,i^*) &\\cong \\ol{P'}(i,i',i^*) \\rightarrow \\ol{P'}(i',i',i^*;(2)) \\cong \\ol{P'}(i^*,i',i';(2)) \\rightarrow \\ol{P'}(i,i',i';(2))",
     "\\ol{P'}(i',i,i^*) &\\cong \\ol{P'}(i^*,i,i') \\rightarrow \\ol{P'}(i,i,i';(2)) \\rightarrow \\ol{P'}(i',i,i';(2)) \\cong  \\ol{P'}(i,i',i';(2))",
     "\\ol{P'}(i',i,i^*) &\\cong \\ol{P'}(i^*,i,i') \\rightarrow \\ol{P'}(i,i,i';(1,1)) \\rightarrow \\ol{P'}(i',i,i';(2)) \\cong  \\ol{P'}(i,i',i';(2))");

  checkThree(Hom(fkiis)*Hom(isokii_iik)*Hom(fiiks)*Hom(isojik_ijk),
             Hom(isokii_iki)*Hom(fikis)*Hom(fiikt)*Hom(isokji_ijk),
             Hom(isokii_iki)*Hom(fikis)*Hom(fiiks)*Hom(isokji_ijk),
             x_ijk,
     "\\ol{P'}(i',i,i^*) &\\cong \\ol{P'}(i,i',i^*) \\rightarrow \\ol{P'}(i',i',i^*;(1,1)) \\cong \\ol{P'}(i^*,i',i';(1,1)) \\rightarrow \\ol{P'}(i,i',i';(1,1))",
     "\\ol{P'}(i',i,i^*) &\\cong \\ol{P'}(i^*,i,i') \\rightarrow \\ol{P'}(i,i,i';(2)) \\rightarrow \\ol{P'}(i',i,i';(1,1)) \\cong  \\ol{P'}(i,i',i';(1,1))",
     "\\ol{P'}(i',i,i^*) &\\cong \\ol{P'}(i^*,i,i') \\rightarrow \\ol{P'}(i,i,i';(1,1)) \\rightarrow \\ol{P'}(i',i,i';(1,1)) \\cong  \\ol{P'}(i,i',i';(1,1))");

  checkTwo(Hom(isokii_iki)*Hom(fikit)*Hom(isokij_ijk)*Hom(fijk),
             Hom(fkiit)*Hom(isokii_iki)*Hom(fikit)*Hom(isokij_ijk),
             x_ijk,
     "\\ol{P'}(j,i,i^*) &\\rightarrow \\ol{P'}(j',i,i^*) \\cong \\ol{P'}(i^*,j',i) \\rightarrow \\ol{P'}(i,j',i;(2)) \\cong \\ol{P'}(j',i,i;(2))",
     "\\ol{P'}(j,i,i^*) &\\cong \\ol{P'}(i^*,j,i) \\rightarrow \\ol{P'}(i,j,i;(2)) \\cong \\ol{P'}(j,i,i;(2)) \\rightarrow\\ol{P'}(j',i,i;(2))");
  checkTwo(Hom(isokii_iki)*Hom(fikis)*Hom(isokij_ijk)*Hom(fijk),
             Hom(fkiis)*Hom(isokii_iki)*Hom(fikis)*Hom(isokij_ijk),
             x_ijk,
     "\\ol{P'}(j,i,i^*) &\\rightarrow \\ol{P'}(j',i,i^*) \\cong \\ol{P'}(i^*,j',i) \\rightarrow \\ol{P'}(i,j',i;(1,1)) \\cong \\ol{P'}(j',i,i;(1,1))",
     "\\ol{P'}(j,i,i^*) &\\cong \\ol{P'}(i^*,j,i) \\rightarrow \\ol{P'}(i,j,i;(1,1)) \\cong \\ol{P'}(j,i,i;(1,1)) \\rightarrow\\ol{P'}(j',i,i;(1,1))");
  checkThree(Hom(isojik_ijk)*Hom(fijk)*Hom(isojik_ijk)*Hom(fijk),
             Hom(fijk)*Hom(fiikt)*Hom(isojik_ijk),
             Hom(fijk)*Hom(fiiks)*Hom(isojik_ijk),
             x_ijk,
     "\\ol{P'}(i,i^*,j) &\\rightarrow \\ol{P'}(i',i^*,j) \\cong \\ol{P'}(i^*,i',j) \\rightarrow \\ol{P'}(i,i',j) \\cong \\ol{P'}(i',i,j)",
     "\\ol{P'}(i,i^*,j) &\\cong \\ol{P'}(i^*,i,j) \\rightarrow \\ol{P'}(i,i,j;(2)) \\rightarrow\\ol{P'}(i',i,j)",
     "\\ol{P'}(i,i^*,j) &\\cong \\ol{P'}(i^*,i,j) \\rightarrow \\ol{P'}(i,i,j;(1,1)) \\rightarrow\\ol{P'}(i',i,j)");
  checkTwo(Hom(isojik_ijk)*Hom(fijk)*Hom(isojik_ijk)*Hom(fijk),
             Hom(fijk)*Hom(isojik_ijk)*Hom(fijk)*Hom(isojik_ijk),
             x_ijk,
     "\\ol{P'}(i,j,k) &\\rightarrow \\ol{P'}(i',j,k) \\cong \\ol{P'}(j,i',k) \\rightarrow \\ol{P'}(j',i',k) \\cong \\ol{P'}(i',j',k)",
     "\\ol{P'}(i,j,k) &\\cong \\ol{P'}(j,i,k) \\rightarrow \\ol{P'}(j',i,k) \\cong \\ol{P'}(i,j',k) \\rightarrow\\ol{P'}(i',j',k)");
  std::cout << makecalctitleij("i","j","k","","","")<<std::endl;
  std::cout << makecalctitleji("i","j","k","","","")<<std::endl;
  std::cout << makecalctitleij("i","i","k","(2)","","")<<std::endl;
}
