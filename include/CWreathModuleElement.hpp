
#ifndef CWREATHMODULEELEMENT_H
#define CWREATHMODULEELEMENT_H
#include <map>
#include "CSymm.hpp"
#include "CTensor.hpp"
// bool operator <(CTensor l,CTensor r);

/*
 * An object of type CWreathModuleElement represents a pure tensor of
 * a module V(i,j,k;lam) as described in the thesis. This module is
 * constructed by first taking the tensor product of a module V of a
 * subgroup of the wreath product G~H for a subgroup H of S3 and
 * an H-module S and then inducing the resulting module to G~S3.
 *
*/

enum class subgroup {i, t, ts, tss, s, ss, S3};
enum class specht {t1, t2, t3, s2, s3, r};
enum class spechtbasis {x1,x2};
CSymm subgroupToSymm(subgroup sub);
subgroup symmToSubgroup(CSymm sy);
int subgroupSize(subgroup sub);

class CWreathModuleElement
{
    public:
        CWreathModuleElement(subgroup type, CSymm coset, int multiplicity,
            CTensor tensor, specht lambda, spechtbasis x);
        subgroup type() const {return mtype;}
        specht lambda() const {return mlambda;}
        spechtbasis x() const {return mx;}
        CSymm coset() const {return mcoset;}
        int multiplicity() const {return mmultiplicity;}
        CTensor tensor() const {return mtensor;}

        void setCoset(CSymm coset);
        void setMultiplicity(int m) {mmultiplicity = m;}
        void setTensor(CTensor tensor) {mtensor = tensor;}
        void setX(spechtbasis x) {mx = x;}

        std::string print() const;
    private:
        subgroup mtype;
        CSymm mcoset;
        int mmultiplicity;
        CTensor mtensor;
        specht mlambda;
        spechtbasis mx;
};

bool operator <(CWreathModuleElement lhs, CWreathModuleElement rhs);

/*
 * The subgroup H is describe by the enumeration subgroup.
 * As all non-trivial subgroups of S3 are cyclic, we represent them
 * by their generator, e.g. subgroup::ts represents G\wr<ts>.
 *
 * The module S is described by an element mlambda of the enumeration
 * specht. Note that S will always be a tensor product of specht
 * modules of subgroups of S3 and at most one of those will be of
 * a subgroup other than S1. If S=S^(1) (x) S^(1) (x) S^(1) then
 * lambda will be t1 (meaning the partition describing the trivial
 * module of S1). Otherwise lambda will describe the partition of
 * the module of Sn with n>1. Here the partitions have the following
 * meanings:
 * t1 = (1), t2 = (2), t3 = (3) (trivial representation for S1,S2,S3)
 * s2 = (1,1), s3 = (1,1,1) (signum representation for S2,S3)
 * r = (2,1) (regular representation of S3)
 *
 * The spechtbasis object mx represents a basis vector of S. As
 * S will always be at most 2-dimensional, the enumeration only contains
 * two elements.
 *
 * The rest of the data consists of a CTensor object which describes the
 * element of the underlying G^3-module V, a permutation mcoset of type
 * CSymm representing the forth tensor product component and
 * a scalar factor mmultiplicity which we here assume to be an integer.
 *
 */

#endif
