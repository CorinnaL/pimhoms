#include "CLinComb.hpp"
#include <iostream>
#include <assert.h>

CTensor operator *(CTensor l,CSymm r);


CLinComb::CLinComb(subgroup type, specht lambda):mtype(type),mlambda(lambda)
{
}

void CLinComb::collect()
{
    std::list<CWreathModuleElement> newsummands;
    for(auto& oldsummand : msummands)
    {
      bool added = false;
      for(auto& summand : newsummands)
      {
          if(summand.coset() == oldsummand.coset() && summand.tensor() == oldsummand.tensor()
              && summand.x() == oldsummand.x())
          {
              summand.setMultiplicity(summand.multiplicity()+oldsummand.multiplicity());
              added = true;
          }
      }
      if(!added)
      {
        newsummands.push_back(oldsummand);
      }
    }
    msummands = newsummands;
    msummands.sort();
}

void CLinComb::add(CWreathModuleElement addsummand)
{
    std::list<CWreathModuleElement> summandlist = normalize(addsummand);
    for(auto& newsummand : summandlist)
    {
      if(newsummand.type()!=mtype || newsummand.lambda()!=mlambda)
      {
          std::cout<<"Incompatible element. Nothing was added."<<std::endl;
          return;
      }
      msummands.push_back(newsummand);
    }
    this->collect();
}


std::string CLinComb::print()
{
    std::string output;
    for(auto& summand : msummands)
    {
        if(summand.multiplicity()<0&&output!="")
        {
            output.pop_back();
        }
        if(summand.multiplicity()!=0)
        {
            output+=summand.print();
            output+='+';
        }
    }
    if(!output.empty()) output.pop_back();
    else return "0";
    return output;
}

void CLinComb::CreateAndAdd(CTensor ten, CSymm pi, int mult, spechtbasis x)
{
        CWreathModuleElement newsummand(mtype, pi, mult, ten, mlambda, x);
        this->add(newsummand);
}

void CLinComb::CreateAndAdd(std::string first, std::string second,
    std::string third, CSymm pi, int mult, spechtbasis x)
{
        CTensor newtensor(first,second,third);
        this->CreateAndAdd(newtensor,pi,mult,x);
}


std::list<CWreathModuleElement> normalize(CWreathModuleElement element)
{
    std::list<CWreathModuleElement> out;
    CWreathModuleElement temp = element;
    subgroup type = element.type();
    specht lambda = element.lambda();
    if(type == subgroup::S3)
    {
        if(lambda == specht::r)
        {
          if(element.coset() == "")
          {
            out.push_back(temp);
          }
          else if(element.coset() =="s")
          {
            if(element.x() == spechtbasis::x1)
            {
                temp.setTensor(element.tensor()*element.coset());
                temp.setMultiplicity(element.multiplicity()*-1);
                temp.setCoset("");
                out.push_back(temp);
                temp.setX(spechtbasis::x2);
                out.push_back(temp);
            }
            if(element.x() == spechtbasis::x2)
            {
                temp.setTensor(element.tensor()*element.coset());
                temp.setX(spechtbasis::x1);
                temp.setCoset("");
                out.push_back(temp);
            }
          }
          else if(element.coset() =="t")
          {
            if(element.x() == spechtbasis::x1)
            {
                temp.setTensor(element.tensor()*element.coset());
                temp.setMultiplicity(element.multiplicity()*-1);
                temp.setCoset("");
                out.push_back(temp);
            }
            if(element.x() == spechtbasis::x2)
            {
                temp.setTensor(element.tensor()*element.coset());
                temp.setCoset("");
                out.push_back(temp);
                temp.setX(spechtbasis::x1);
                out.push_back(temp);
            }
          }
          else
          {
            auto factors = factorize(element.coset());
            std::list<CWreathModuleElement> out2;
            out2.push_back(temp);
            for(auto& fac : factors)
            {
              out = out2;
              out2.clear();
              for(auto& elt : out)
              {
                elt.setCoset(fac);
                for(auto& elt2 : normalize(elt))
                {
                  out2.push_back(elt2);
                }
              }
            }
            out = out2;
          }
        }
        else
        {
          if(lambda==specht::s3)
          {
            temp.setMultiplicity(element.multiplicity()*signum(element.coset()));
          }
          temp.setTensor(element.tensor()*element.coset());
          temp.setCoset("");
          out.push_back(temp);
        }
    }
    else if(type == subgroup::i)
    {
      out.push_back(temp);
    }
    else  // type is a transposition
    {
        assert(type == subgroup::t || type == subgroup::ts || type == subgroup::tss);
        if( element.coset() == "" || element.coset() == "s" || element.coset() == "ss")
        {
            out.push_back(temp);
        }
        else
        {
            temp.setCoset(subgroupToSymm(type)*element.coset());
            temp.setTensor(element.tensor()*subgroupToSymm(type));
            if(lambda == specht::s2 || lambda == specht::s3)
            {
              temp.setMultiplicity(element.multiplicity()*-1);
            }
            out.push_back(temp);
          /*
           * If the module is induced from G~<tau>, then we want
           * to bring the elements to a normal form where the tensor
           * factor of the induction is always given by s or s^2.
           * Since type represents the generator (and only non-trivial
           * element) of the subgroup from which we induce and the 
           * transpositions are selft-inverse
           * a (x) b (x) c (x) x1 (x) coset =
           * a (x) b (x) c (x) x1 (x) type*type*coset
           * by interpreting type as an element of S3.
           * Since we assumed that coset \notin {i,s,ss}, we know that
           * type*coset \in {i,s,ss}.
           * The only thing left to note is that if lambda represents
           * a signum representation, the action on S will multiply the
           * element by -1 (as x1*type = -x1).
           * Therefore the above calculation
           * gives us the intended element in the chosen normal form.
           */
        }
    }
    return out;
}

void CLinComb::mult(CSymm r)
{
    std::list<CWreathModuleElement> newsummands;
    for(auto & summand : msummands)
    {
      summand.setCoset(summand.coset()*r);
      auto nomlist = normalize(summand);
      for( auto& sum : nomlist)
      {
        newsummands.push_back(sum);
      }
    }
    msummands = newsummands;
    this->collect();
}

CLinComb operator*(const CLinComb& l,const CSymm& r)
{
    CLinComb output(l.type(),l.lambda());
    auto oldsummands = l.summands();
    for(auto summand : oldsummands)
    {
      summand.setCoset(summand.coset()*r);
      output.add(summand);
    }
    return output;
}
