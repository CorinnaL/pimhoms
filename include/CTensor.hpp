#ifndef CTENSOR_H
#define CTENSOR_H
#include <string>
/*
 * An object of type CTensor is our representation of an element of
 * a tensor product of three modules. The intended interpretation
 * on the variables is such that the element is
 * firstcomp (X) secondcomp (X) thirdcomp.
 *
 * The only structural assumption we make is that the empty string
 * represents the zero-element of the underlying vector spaces.
 *
 * The print function returns the LaTeX code describing the
 * intended tensor product.
*/

class CTensor
{
    public:

        CTensor(std::string first, std::string second, std::string third);
        CTensor();
        void setComponents(std::string first, std::string second, std::string third);
        std::string first() const {return firstcomp;}
        std::string second() const {return secondcomp;}
        std::string third() const {return thirdcomp;}
        std::string print() const;
        bool isZero() const {return firstcomp=="" || secondcomp == "" || thirdcomp == "";}
    private:
        std::string firstcomp;
        std::string secondcomp;
        std::string thirdcomp;
};

bool operator <(const CTensor l,const CTensor r);
bool operator ==(const CTensor l, const CTensor r);
#endif // CTENSOR_H
