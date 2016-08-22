#ifndef CSYMM_H
#define CSYMM_H
#include <string>
#include <list>

/*
 * An object of type CSymm describes an element of the symmetric group
 * S3. The element is represented by a string where the empty string
 * is the identity element, "s" describes the permutation (1,2,3)
 * and "t" the transposition (1,2).
 *
 * As the actual string is private it can only be changed via the
 * set Element function. This function calls the Simplify() function
 * which uses the relations st = ts^2, s^3=1 and t^2=1 to transform
 * string to take the form t^as^b where a is 0 or 1 and b is between 0 and 2.
 *
 * The constructors also call the Simplify() function so it can always be
 * assumed that the string element is of the form above.
 *
 * The print function returns the permutation in cycle notation
 * except the identity function which is printed as "id".
*/


class CSymm
{
    public:
        CSymm();
        CSymm(std::string s);
        CSymm(char const* s);
        std::string get() const { return element; }
        void setElement(std::string val);
        void Simplify();
        std::string print() const;
    private:
        std::string element;
};

int signum(const CSymm&);
const CSymm operator *(const CSymm&,const CSymm&);
std::list<CSymm> factorize(const CSymm&);
CSymm inverse(const CSymm& l);
bool operator <(const CSymm& l,const CSymm& r);
bool operator ==(const CSymm& l,const CSymm& r);
#endif // CSYMM_H
