#ifndef CLINCOMB_H
#define CLINCOMB_H
#include <list>
#include <string>
#include <iostream>
#include "CWreathModuleElement.hpp"

class CLinComb
{
    public:
        CLinComb(subgroup type,specht lambda);
        subgroup type() const { return mtype; }
        specht lambda() const { return mlambda; }
        std::list<CWreathModuleElement> summands() const {return msummands;}
        void add(CWreathModuleElement);
        std::string print();
        void CreateAndAdd(CTensor ten, CSymm pi, int mult, 
            spechtbasis x = spechtbasis::x1);
        void CreateAndAdd(std::string first, std::string second,
            std::string third, CSymm pi, int mult,
            spechtbasis x = spechtbasis::x1);
        void collect();
        void mult(CSymm r);
    private:
        subgroup mtype;
        specht mlambda;
        std::list<CWreathModuleElement> msummands;
};

std::list<CWreathModuleElement> normalize(CWreathModuleElement element);
CLinComb operator*(const CLinComb& l,const CSymm& r);
#endif // CLINCOMB_H
