#include "FileReader.hh"
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip> //std::setw

void FileReader::registerIntParameter(const std::string &key, int init)
{
   IntParameter[key] = init;
}

void FileReader::registerRealParameter(const std::string &key, real init)
{
   RealParameter[key] = init;
}

void FileReader::registerStringParameter(const std::string &key, const std::string &init)
{
   StringParameter[key] = init;
}

void FileReader::setParameter(const std::string &key, const std::string &in)
{
   CHECK_MSG(StringParameter.find(key) != StringParameter.end(), "Parameter '" + key + "\"' not registered!");
   StringParameter[key] = in;
}

void FileReader::setParameter(const std::string &key, real in)
{
   CHECK_MSG(RealParameter.find(key) != RealParameter.end(), "Parameter '" + key + "\"' not registered!");
   RealParameter[key] = in;
}

void FileReader::setParameter(const std::string &key, int in)
{
   CHECK_MSG(IntParameter.find(key) != IntParameter.end(), "Parameter '" + key + "\"' not registered!");
   IntParameter[key] = in;
}


bool FileReader::readFile(const std::string &name)
{
   std::ifstream paramFile(name);
   CHECK_MSG(paramFile,"Could not open file '"+ name + "' which has to be in the current directory");
   std::string line, col1, col2;
   while (std::getline(paramFile, line)) {
       std::istringstream ss(line, std::istringstream::in);
       ss >> col1 >> col2;

       if ( (!col1.empty()) && col1.at(0) != '#' && (!col2.empty()) && col2.at(0) != '#') {
          if(string_only_digits(col2) )
             registerIntParameter(col1,std::stoi(col2));
          else if(string_only_reals(col2) )
             registerRealParameter(col1,std::stof(col2));
          else
             registerStringParameter(col1,col2); }

   }   return true;
}

const bool FileReader::find( const std::string & key ) const
{
   if (IntParameter.find(key) != IntParameter.end())
       return true;
   else if (RealParameter.find(key) != RealParameter.end())
       return true;
   else if (StringParameter.find(key) != StringParameter.end())
       return true;
   else
       return false;
}

bool FileReader::has_only_spaces(const std::string& str) {
	if(str.find_first_not_of (' ') == str.npos||str.find_first_not_of ('\t') == str.npos)
		return true;
	else
		return false;
}

bool FileReader::string_only_digits(const std::string& str){
    return str.find_first_not_of("0123456789-") == std::string::npos;
}

bool FileReader::string_only_reals(const std::string& str){
    return (   str.find_first_not_of("0123456789.e-") == std::string::npos 
            || str.find_first_not_of("0123456789.E-") == std::string::npos );
}

void FileReader::printParameters() const
{
    std::cout << "Printing out all integer parameters: \n";
    for (auto i = IntParameter.begin(); i != IntParameter.end(); ++i)
    {
	std::cout << i->first << " " << i->second << std::endl;
    }
    std::cout << std::endl;
    std::cout << "Printing out all real parameters: \n";
    for (auto i = RealParameter.begin(); i != RealParameter.end(); ++i)
    {
	std::cout << std::left<< std::setw(20)<< i->first << " " << i->second << std::endl;
    }
    std::cout << std::endl;
    std::cout << "Printing out all string parameters: \n";
    for (auto i = StringParameter.begin(); i != StringParameter.end(); ++i)
    {
	std::cout << std::left<< std::setw(20)<< i->first << " " << i->second << std::endl;
    }
    std::cout << std::endl;
    std::cout << "Finished printing all parameters \n\n";
}
