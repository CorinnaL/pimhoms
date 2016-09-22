#include <string>
#include "CTensor.hpp"
#include "CSymm.hpp"
#include "functions.hpp"
#include <iostream>

extern bool latex;
std::string texify(std::string i)
{
    if(latex && i.size()==2 && i.substr(1,2)=="*")
    {
        return i.substr(0,1)+"^*";
    }
    else
    {
        return i;
    }

}
std::string iprime(std::string i)
{
    if(i.size() ==1)
        return i+"'";
    else if(i.size()==2 &&i.substr(1,2)=="*")
        return i.substr(0,1);
    else if(i.size()==3 &&i.substr(1,3)=="^*")
        return i.substr(0,1);
    std::cout<<"No method to apply ' to this string, beware!"<<std::endl;
    return i;
}
std::string projectiveIndec(std::string i, std::string j, std::string k, std::string lam)
{
    if(latex)
    {
        if(lam != "")
        {
            return "\\ol{P'}("+i+","+j+","+k+";"+lam+")";
        }
        return "\\ol{P'}("+i+","+j+","+k+")";
    }
    return i+j+k+lam;
}
std::string makecalctitleij(std::string i,std::string j,std::string k,
        std::string firstlam,std::string secondlam, std::string thirdlam)
{
    std::string firstind = projectiveIndec(i,j,k,firstlam);
    std::string secondind = projectiveIndec(iprime(i),j,k,secondlam);
    std::string thirdind = projectiveIndec(j,iprime(i),k,secondlam);
    std::string fourthind = projectiveIndec(iprime(j),iprime(i),k,thirdlam);
    std::string fifthind = projectiveIndec(iprime(i),iprime(j),k,thirdlam);
    std::string arrow = "->";
    std::string cong = "=";
    if(latex)
    {
        arrow = "\\rightarrow";
        cong = "\\cong";
    }
    if(i==j)
    {
        return firstind+arrow+secondind+cong+thirdind+arrow+fourthind;
    }
    return firstind+arrow+secondind+cong+thirdind+arrow+fourthind+cong+fifthind;
}

std::string makecalctitleji(std::string i,std::string j,std::string k,
        std::string firstlam,std::string secondlam, std::string thirdlam)
{
    std::string firstind = projectiveIndec(i,j,k,firstlam);
    std::string secondind = projectiveIndec(j,i,k,firstlam);
    std::string thirdind = projectiveIndec(iprime(j),i,k,secondlam);
    std::string fourthind = projectiveIndec(i,iprime(j),k,secondlam);
    std::string fifthind = projectiveIndec(iprime(i),iprime(j),k,thirdlam);
    std::string arrow = "->";
    std::string cong = "=";
    if(latex)
    {
        arrow = "\\rightarrow";
        cong = "\\cong";
    }
    return firstind+cong+secondind+arrow+thirdind+cong+fourthind+arrow+fifthind;
}

Hom operator*(Hom hom1,
    Hom hom2)
{
  return [hom1,hom2] (CLinComb a) {return hom1(hom2(a));};
}

void checkThree(Hom hom1,
        Hom hom2,
        Hom hom3,
        CLinComb x,
        std::string i, std::string j, std::string k,
        std::string firstlam, std::string secondlam1, std::string secondlam2,
        std::string secondlam3, std::string thirdlam, comptype twoitwoj)
{
    std::string firstcomp;
    std::string secondcomp;
    std::string thirdcomp;
    if(twoitwoj==comptype::iij)
    {
        firstcomp = makecalctitleij(i,j,k,firstlam,secondlam1,thirdlam);
        secondcomp = makecalctitleij(i,j,k,firstlam,secondlam2,thirdlam);
        thirdcomp = makecalctitleji(i,j,k,firstlam,secondlam3,thirdlam);
    }
    else
    {
        firstcomp = makecalctitleij(i,j,k,firstlam,secondlam1,thirdlam);
        secondcomp = makecalctitleji(i,j,k,firstlam,secondlam2,thirdlam);
        thirdcomp = makecalctitleji(i,j,k,firstlam,secondlam3,thirdlam);
    }

    std::string linestart = "-------------------------";
    std::string lineend = linestart;
    std::string calcend = "***************************************************************";
    if(latex)
    {
        linestart = "";
        lineend = "\\\\";
        calcend = "\\end{align*}\n\\begin{align*}";
    }
        std::cout << linestart <<firstcomp<<lineend <<std::endl;
        std::cout << x.print() << std::endl;
        CLinComb c1 = hom1(x);
        std::cout << linestart <<secondcomp<<lineend <<std::endl;
        std::cout << x.print() << std::endl;
        CLinComb c2 = hom2(x);
        std::cout << linestart <<thirdcomp<<lineend <<std::endl;
        std::cout << x.print() << std::endl;
        CLinComb c3 = hom3(x);
        std::cout << calcend <<std::endl;
}



void checkTwo(Hom hom1, Hom  hom2, CLinComb x,
        std::string i, std::string j, std::string k, std::string firstlam,
        std::string secondlam1,std::string secondlam2, std::string thirdlam)
{
    std::string firstcomp;
    std::string secondcomp;
    if(i==j)
    {
        firstcomp = makecalctitleij(i,j,k,firstlam,secondlam1,thirdlam);
        secondcomp = makecalctitleij(i,j,k,firstlam,secondlam2,thirdlam);
    }
    else
    {
        firstcomp = makecalctitleij(i,j,k,firstlam,secondlam1,thirdlam);
        secondcomp = makecalctitleji(i,j,k,firstlam,secondlam2,thirdlam);
    }

    std::string linestart = "--------------------------";
    std::string lineend = linestart;
    std::string calcend = "***************************************************************";
    if(latex)
    {
        linestart = "";
        lineend = "\\\\";
        calcend = "\\end{align*}\n\\begin{align*}";
    }
        std::cout << linestart <<firstcomp<<lineend <<std::endl;
        std::cout << x.print() << std::endl;
        CLinComb c1 = hom1(x);
        std::cout << linestart <<secondcomp<<lineend <<std::endl;
        std::cout << x.print() << std::endl;
        CLinComb c2 = hom2(x);
        std::cout << calcend <<std::endl;
}


std::string findreplace(std::string input, std::string find, std::string replace)
{
    std::size_t findsize = find.size();
    while(input.find(find)!=std::string::npos)
    {
        input.replace(input.find(find),findsize,replace);
    }
    return input;
}


CTensor operator *(CTensor l,CSymm r)
{
    std::string symmstring = r.get();
    CTensor result = l;
    while(symmstring!="")
    {
        if(symmstring.front()=='s')
        {
            result.setComponents(result.second(), result.third(), result.first());
        }
        else if(symmstring.front()=='t')
        {
            result.setComponents(result.second(), result.first(), result.third());
        }
        symmstring = symmstring.substr(1,symmstring.size());
    }
    return result;
}

