#ifndef FILEREADER_HH
#define FILEREADER_HH

#include <string>
#include <map>
#include "Debug.hh"
#include "Types.hh"


//*******************************************************************************************************************
/*! Class for reading configuration from a file
*
* Configuration File Syntax:
*   - everything between a '#' character and the beginning of the next line is ignored
*   - lines can be empty
*   - line contain parameter key followed by white-spaces followed by their value
*
*  All possible keys (and the datatype of the value) have to be registered first:
*   - for example usage have a look at the FileReaderTest
*/
//*******************************************************************************************************************
class FileReader
{
private:

	std::map<std::string, std::string> parameters;

public:

	//register a new parameter with name key and initial int value
	void registerIntParameter( const std::string & key, int init = 0 );

	//register a new parameter with name key and initial double value
	void registerRealParameter( const std::string & key, real init = 0.0 );

	//register a new parameter with name key and initial string value
	void registerStringParameter( const std::string & key, const std::string & init = "" );

	//set a value for the key string with value in
	void setParameter( const std::string & key, const std::string & in );

	//set a value for the key string with value in
	void setParameter( const std::string & key, real in );

	//set a value for the key string with value in
	void setParameter( const std::string & key, int in );

	// get the int value of key 
	inline int getIntParameter( const std::string & key ) const;

	// get the double value of key 
	inline real getRealParameter( const std::string & key ) const;

	// get the string value of key 
	inline std::string getStringParameter( const std::string & key ) const;

	//try to read all registered parameters from file name
	bool readFile( const std::string & name );

	//print out all parameters to std:out
	void printParameters() const;
	
	//find if a parameter is registered
	const bool find( const std::string & key ) const;

};




inline int FileReader::getIntParameter(const std::string &key) const
{
   auto iter = parameters.find(key);
   int val=0;
   CHECK_MSG(iter != parameters.end(), "Parameter '" + key +"' does not exist!");
   try {val = std::stoi(iter->second);}
   catch (...)
       {CHECK_MSG(0, "Parameter type of '" + key +"' is not int!");}
   return val;
}

inline real FileReader::getRealParameter(const std::string &key) const
{
   auto iter = parameters.find(key);
   real val=0.0;
   CHECK_MSG(iter != parameters.end(), "Parameter '" + key +"' does not exist!");
   try {val = std::stof(iter->second);}
   catch (...)
       {CHECK_MSG(0, "Parameter type of \"" + key +"\" is not int!");}
   return val;
}

inline std::string FileReader::getStringParameter(const std::string &key) const
{
   auto iter = parameters.find(key);
   std::string val = "";
   CHECK_MSG(iter != parameters.end(), "Error: Parameter '" + key +"' does not exist!");
   val = iter->second;
   return val;
}





#endif //FILEREADER_HH

