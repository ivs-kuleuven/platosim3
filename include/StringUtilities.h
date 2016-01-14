//
//  StringUtilities.h
//  New PLATO Simulator
//
//  Created by Rik Huygen on 27/11/15.
//  Copyright © 2015 KU Leuven. All rights reserved.
//

#ifndef STRING_UTILITIES_H
#define STRING_UTILITIES_H

#include <iostream>
#include <string>
#include <vector>
#include <cstddef>


static bool ends_with(std::string const & value, std::string const & ending)
{
    if (ending.size() > value.size())
        return false;
    return std::equal(ending.rbegin(), ending.rend(), value.rbegin());
}


struct split
{
    enum empties_t { empties_ok, no_empties };
};

template <typename Container> Container& split(
    Container& result,
    const typename Container::value_type& s,
    const typename Container::value_type& delimiters,
    split::empties_t empties = split::empties_ok )
{
    result.clear();
    size_t current;
    size_t next = -1;
    do
    {
        if (empties == split::no_empties)
        {
            next = s.find_first_not_of( delimiters, next + 1 );
            if (next == Container::value_type::npos) 
                break;
            next -= 1;
        }
        current = next + 1;
        next = s.find_first_of( delimiters, current );
        result.push_back( s.substr( current, next - current ) );
    }
    while (next != Container::value_type::npos);
    return result;
}


static void print( std::vector <std::string> & vector )
{
    for (size_t n = 0; n < vector.size(); n++)
        std::cout << "\"" << vector[ n ] << "\"" << std::endl;
    std::cout << std::endl;
}


#endif /* STRING_UTILITIES_H */
