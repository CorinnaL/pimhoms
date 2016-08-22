#ifndef HOMOMORPHISM_HPP_INCLUDED
#define HOMOMORPHISM_HPP_INCLUDED
#include "CLinComb.hpp"
#include "CWreathModuleElement.hpp"

CLinComb fiiir(CLinComb in);
CLinComb fiiis(CLinComb in);
CLinComb fiiit(CLinComb in);
CLinComb fkiit(CLinComb in);
CLinComb fkiis(CLinComb in);
CLinComb fikis(CLinComb in);
CLinComb fikit(CLinComb in);
CLinComb fiiks(CLinComb in);
CLinComb fiikt(CLinComb in);
CLinComb fijk(CLinComb in);

CLinComb isokii_iik(CLinComb in);
CLinComb isoiki_iik(CLinComb in);

CLinComb isokii_iki(CLinComb in);
CLinComb isoiik_iki(CLinComb in);

CLinComb isoiki_kii(CLinComb in);
CLinComb isoiik_kii(CLinComb in);

CLinComb isojik_ijk(CLinComb in);
CLinComb isoikj_ijk(CLinComb in);
CLinComb isokji_ijk(CLinComb in);
CLinComb isokij_ijk(CLinComb in);
CLinComb isojki_ijk(CLinComb in);

CLinComb changeFirstGen(CLinComb in);
#endif
