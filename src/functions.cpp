#include <string>
#include "CTensor.hpp"
#include "CSymm.hpp"
#include <iostream>

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

