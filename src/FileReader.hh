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

    bool has_only_spaces(const std::string& str);
    bool string_only_digits(const std::string& str);
    bool string_only_reals(const std::string& str);
    std::map<std::string,int> IntParameter;
    std::map<std::string,real> RealParameter;
    std::map<std::string,std::string> StringParameter;

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
   auto iter = IntParameter.find(key);
   CHECK_MSG(iter != IntParameter.end(), "Parameter '" + key +"' does not exist!");
   return iter->second;
}

inline real FileReader::getRealParameter(const std::string &key) const
{
   auto iter = RealParameter.find(key);
   CHECK_MSG(iter != RealParameter.end(), "Parameter '" + key +"' does not exist!");
   return iter->second;
}

inline std::string FileReader::getStringParameter(const std::string &key) const
{
   auto iter = StringParameter.find(key);
   CHECK_MSG(iter != StringParameter.end(), "Parameter '" + key +"' does not exist!");
   return iter->second;
}

#endif //FILEREADER_HH
