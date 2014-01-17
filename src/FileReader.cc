#include "FileReader.hh"
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip> //std::setw



void FileReader::registerIntParameter(const std::string &key, int init)
{
   std::stringstream ss;
   ss << init;
   parameters[key] = ss.str();
}

void FileReader::registerRealParameter(const std::string &key, real init)
{
   std::stringstream ss;
   ss << init;
   parameters[key] = ss.str();
}

void FileReader::registerStringParameter(const std::string &key, const std::string &init)
{
   parameters[key] = init;
}

void FileReader::setParameter(const std::string &key, const std::string &in)
{
   CHECK_MSG(parameters.find(key) != parameters.end(), "Parameter '" + key + "\"' not registered!");
   parameters[key] = in;
}

void FileReader::setParameter(const std::string &key, real in)
{
   std::stringstream ss;
   ss << in;
   CHECK_MSG(parameters.find(key) != parameters.end(), "Parameter '" + key + "\"' not registered!");
   parameters[key] = ss.str();
}

void FileReader::setParameter(const std::string &key, int in)
{
   std::stringstream ss;
   ss << in;
   CHECK_MSG(parameters.find(key) != parameters.end(), "Parameter '" + key + "\"' not registered!");
   parameters[key] = ss.str();
}


bool FileReader::readFile(const std::string &name)
{
   std::ifstream paramFile(name);
   CHECK_MSG(paramFile,"Could not open file '"+ name + "' which has to be in the current directory");
   std::string line, col1, col2;
   while (std::getline(paramFile, line)) {
       std::istringstream ss(line, std::istringstream::in);
       ss >> col1 >> col2;

       if ( !col1.empty() && col1.at(0) != '#' && !col2.empty() && col2.at(0) != '#') {
          if(parameters.find(col1) == parameters.end())
             registerStringParameter(col1,col2);
          else {
             try {stod(parameters[col1]);}
             catch(...){
                setParameter(col1,col2);
                continue;}
             try {stod(col2);}
             catch(...){
                CHECK_MSG(0, "Type mismatch for parameter '" + col1 +"' "); }
             setParameter(col1,real(stod(col2))); } }
   }   return true;
}

const bool FileReader::find( const std::string & key ) const
{
   if (parameters.find(key) != parameters.end())
       return true;
   else
       return false;
}

void FileReader::printParameters() const
{
   for(auto it = parameters.begin(); it != parameters.end(); ++it)
        std::cout << std::left<< std::setw(20)<< it->first << " " << it->second << "\n";
}
