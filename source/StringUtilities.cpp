#include "StringUtilities.h"




bool StringUtilities::ends_with(string const & value, string const & ending)
{
    if (ending.size() > value.size())
        return false;
    return equal(ending.rbegin(), ending.rend(), value.rbegin());
}




vector<string> StringUtilities::split(string myString, char delimiter)
{
   vector<string> parts;
   string part;

   istringstream myStream(myString);
   while(getline(myStream, part, delimiter)) 
   {
      parts.push_back(part);
   }

  return parts;
}


string StringUtilities::dtos(double value, bool scientific)
{
    stringstream os;

    if (scientific)
    {
        os << std::scientific;
    }
    else
    {
        os << fixed;
    }
    
    os << showpoint;
    os << setprecision(6);
    os << value;

    return os.str();
}



void StringUtilities::print( std::vector <std::string> & vector )
{
    for (size_t n = 0; n < vector.size(); n++)
        std::cout << "\"" << vector[ n ] << "\"" << std::endl;
    std::cout << std::endl;
}
