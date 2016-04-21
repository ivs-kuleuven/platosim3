#include "StringUtilities.h"



namespace StringUtilities
{

    bool ends_with(string const & value, string const & ending)
    {
        if (ending.size() > value.size())
            return false;
        return equal(ending.rbegin(), ending.rend(), value.rbegin());
    }
    
    vector<string> split(string myString, char delimiter)
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
    
    
    string dtos(double value, bool scientific, int precision)
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
        os << setprecision(precision);
        os << value;
    
        return os.str();
    }
    

    void print( std::vector <std::string> & vector )
    {
        for (size_t n = 0; n < vector.size(); n++)
            std::cout << "\"" << vector[ n ] << "\"" << std::endl;
        std::cout << std::endl;
    }

    
    /**
     * @brief      Replace the ENV['var'] pattern with the value of the environment variable 'var'.
     * 
     * @details    The given string might contain a pattern that represents an environment variable
     *             which needs to be replaced by its actual value. The pattern that is supported is 
     *             ENV['var'] where var is the name of the environment variable. 
     *             
     *             The name can only contain uppercase characters and the underscore.
     *             
     *             Whenever an error occurs, the original string will be returned and a warning message
     *             will be logged, no replacements will be done.
     *             
     *             The inputString can only contain one pattern to match.
     *             
     * @example    PLATOSIM_HOME=/Users/rik/Git/PlatoSim3
     *             
     *             the ENV['PLATOSIM_HOME']/inputfiles will become /Users/rik/Git/PlatoSim3/inputfiles
     *
     * @param[in]  inputString  the string containing the pattern
     * 
     * @return     a copy of the string where the pattern has been replaced
     */
    string replaceEnvironmentVariable(const string inputString)
    {
        string outputString;
        string preMatch;       // contains the part of inputString before ENV['var']
        string match;          // contains the var in ENV['var']
        string postMatch;      // contains the part of inputString after ENV['var']
    
        const string prefix("ENV['");
        const string suffix("']");
    
        string::size_type idx_prefix = inputString.find(prefix);
        string::size_type idx_suffix = inputString.find(suffix);
    
        if (idx_prefix == string::npos)
        {
            // inputString does not contain an environment variable
            
            outputString = inputString;
        }
        else
        {
            if (idx_suffix == string::npos)
            {
                // inputString does not contain the closing part of the ENV['var'] pattern
    
                Log.warning("replaceEnvironmentVariable: closing part of the ENV['var'] pattern is not present, returning original string.\n"
                    "inputString: " + inputString);
    
                outputString = inputString;
            }
            preMatch = inputString.substr(0, idx_prefix);
            match = inputString.substr(idx_prefix + prefix.size(), idx_suffix - (idx_prefix + prefix.size()));
            postMatch = inputString.substr(idx_suffix + suffix.size());
    
            // getenv will return a NULL when the environment variable doesn't exist
            // We need to catch this before trying to convert to string (as the
            // conversion will call strlen which will crash on NULL).
            
            const char * result = getenv(match.c_str());
            if (!result)
            {
                Log.warning("replaceEnvironmentVariable: Environment variable " + match + " doesn't exist.");
                return inputString;
            }
    
            string envReplacement = string(result);
    
            outputString = preMatch + envReplacement + postMatch;
        }
    
    
        return outputString;
    }


//  This is the implementation of replaceEnvironmentVariable using regular expressions,
//  but the code failed on GCC v4.8.5. This only works as off GCC v4.9.2.
// 
//  string replaceEnvironmentVariable(const string inputString)
//  {
//      const regex re(R"(ENV\['([A-Z_]+)'\])");
//      smatch match;
//      string outputString;
//      
//      // Check if we have a match, i.e. if inputString contains the required pattern
///      regex_search(inputString, match, re);
///      if (match.empty())
//          return inputString;
//  
//      if (match.size() != 2)
//      {
//          Log.warning("replaceEnvironmentVariable: expected one match for an environment variable, got " + to_string(match.size()-1) + " matches.");
//          return inputString;
//      }
//  
//      // match[0] will contain the complete match including ENV[' and ']
//      // match[1] will contain only the name of the environment variable
//  
//      string envName = match[1];
//  
//      // getenv will return a NULL when the environment variable doesn't exist
//      // We need to catch this before trying to convert to string (as the
//      // conversion will call strlen which will crash on NULL).
//      
//      const char * result = getenv(envName.c_str());
//      if (!result)
//      {
//          Log.warning("replaceEnvironmentVariable: Environment variable " + envName + " doesn't exist.");
//          return inputString;
//      }
//  
//      string envReplacement = string(result);
//  
//      outputString = regex_replace(inputString, re, envReplacement);
//      
//      return outputString;
//  }


}

