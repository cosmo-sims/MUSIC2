/*
 
 config_file.hh - This file is part of MUSIC -
 a code to generate multi-scale initial conditions 
 for cosmological simulations 
 
 Copyright (C) 2010  Oliver Hahn
 
 */

#ifndef __CONFIG_FILE_HH
#define __CONFIG_FILE_HH

#include <string>
#include <sstream>
#include <map>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <typeinfo>

#include "log.hh"

/*!
 * @class config_file
 * @brief provides read/write access to configuration options
 *
 * This class provides access to the configuration file. The 
 * configuration is stored in hash-pairs and can be queried and
 * validated by the responsible class/routine
 */
class config_file {
  
  //! current line number
  unsigned m_iLine;

  //! hash table for key/value pairs, stored as strings
  std::map<std::string, std::string> m_Items;

public:
	
  //! removes all white space from string source
  /*!
   * @param source the string to be trimmed
   * @param delims a string of delimiting characters
   * @return trimmed string
   */
  std::string trim(std::string const& source, char const* delims = " \t\r\n") const{
    std::string result(source);
    //... skip initial whitespace ...
    std::string::size_type index = result.find_last_not_of(delims);
    if(index != std::string::npos)
      result.erase(++index);
    //... find beginning of trailing whitespace ...
    index = result.find_first_not_of(delims);
    //... remove trailing whitespace ...
    if(index != std::string::npos)
      result.erase(0, index);
    else
      result.erase();
    return result;
  }
  
  //! converts between different variable types
  /*!
   *  The main purpose of this function is to parse and convert
   *  a string argument into numbers, booleans, etc...
   * @param ival the input value (typically a std::string)
   * @param oval the interpreted/converted value
   */
  template <class in_value, class out_value>
  void convert( const in_value & ival, out_value & oval) const
  {
		std::stringstream ss;
		ss << ival; //.. insert value into stream
		ss >> oval; //.. retrieve value from stream

		if (! ss.eof()) {
		  //.. conversion error
			std::cerr << "Error: conversion of \'" << ival << "\' failed." << std::endl;
		  throw ErrInvalidConversion(std::string("invalid conversion to ")+typeid(out_value).name()+'.');
		}
  }

  //! constructor of class config_file
  /*! @param FileName the path/name of the configuration file to be parsed
   */
  config_file( std::string const& FileName )
	: m_iLine(0), m_Items()	
  {
    std::ifstream file(FileName.c_str());
	  
	if( !file.is_open() )
		throw std::runtime_error(std::string("Error: Could not open config file \'")+FileName+std::string("\'"));

    std::string line;
    std::string name;
    std::string value;
    std::string inSection;
    int posEqual;
    m_iLine=0;
    //.. walk through all lines ..
    while (std::getline(file,line)) {
      ++m_iLine;
      //.. encounterd EOL ?
      if (! line.length()) continue;
      
      //.. encountered comment ?
      unsigned long idx;
      if( (idx=line.find_first_of("#;%")) != std::string::npos )
        line.erase(idx);

      //.. encountered section tag ?
      if (line[0] == '[') {
        inSection=trim(line.substr(1,line.find(']')-1));
        continue;
      }

      //.. seek end of entry name ..
      posEqual=line.find('=');
      name  = trim(line.substr(0,posEqual));
      value = trim(line.substr(posEqual+1));
	
	  if( (size_t)posEqual==std::string::npos && (name.size()!=0||value.size()!=0) )
	  {
		  LOGWARN("Ignoring non-assignment in %s:%d",FileName.c_str(),m_iLine);
		  continue;
	  }
	
	  if(name.length()==0&&value.size()!=0)
	  {  
		  LOGWARN("Ignoring assignment missing entry name in %s:%d",FileName.c_str(),m_iLine);
		  continue;
		  
	  }
		
	  if(value.length()==0&&name.size()!=0)
	  {	  
		  LOGWARN("Empty entry will be ignored in %s:%d",FileName.c_str(),m_iLine);
		  continue;
	  }
		
	  if( value.length()==0&&name.size()==0)
		  continue;
      
      //.. add key/value pair to hash table ..
	  if( m_Items.find(inSection+'/'+name) != m_Items.end() )
		  LOGWARN("Redeclaration overwrites previous value in %s:%d",FileName.c_str(),m_iLine);
	
	  m_Items[inSection+'/'+name] = value;
		
    }
  }
	
	//! inserts a key/value pair in the hash map
	/*! @param key the key value, usually "section/key"
	 *  @param value the value of the key, also a string
	 */
	void insertValue( std::string const& key, std::string const& value )
	{
		m_Items[key] = value;
	}
	
	//! inserts a key/value pair in the hash map
	/*! @param section section name. values are stored under "section/key"
	 *  @param key the key value usually "section/key"
	 *  @param value the value of the key, also a string
	 */
	void insertValue( std::string const& section, std::string const& key, std::string const& value )
	{
		m_Items[section+'/'+key] = value;
	}
	
	//! checks if a key is part of the hash map
	/*! @param section the section name of the key
	 *  @param key the key name to be checked
	 *  @return true if the key is present, false otherwise
	 */
	bool containsKey( std::string const& section, std::string const& key )
	{
		std::map<std::string,std::string>::const_iterator i = m_Items.find(section+'/'+key);
		if ( i == m_Items.end() ) 
			return false;
		return true;
	}
	
	//! checks if a key is part of the hash map
	/*! @param key the key name to be checked
	 *  @return true if the key is present, false otherwise
	 */
	bool containsKey( std::string const& key )
	{
		std::map<std::string,std::string>::const_iterator i = m_Items.find(key);
		if ( i == m_Items.end() ) 
			return false;
		return true;
	}
	
	
	//! return value of a key
	/*! returns the value of a given key, throws a ErrItemNotFound
	*  exception if the key is not available in the hash map.
	*  @param key the key name
	*  @return the value of the key
	*  @sa ErrItemNotFound
	*/
	template<class T> T getValue( std::string const& key ) const{
		return getValue<T>( "", key );
	}

	//! return value of a key
	/*! returns the value of a given key, throws a ErrItemNotFound
	 *  exception if the key is not available in the hash map.
	 *  @param section the section name for the key
	 *  @param key the key name
	 *  @return the value of the key
	 *  @sa ErrItemNotFound
	 */
	template<class T> T getValue( std::string const& section, std::string const& key ) const
	{
		T r;
		std::map<std::string,std::string>::const_iterator i = m_Items.find(section + '/' + key);
		if ( i == m_Items.end() ) 
		  throw ErrItemNotFound('\'' + section + '/' + key + std::string("\' not found."));
		  
		convert(i->second,r);
		return r;
	}

	//! exception safe version of getValue
	/*! returns the value of a given key, returns a default value rather
	 *  than a ErrItemNotFound exception if the key is not found.
	 *  @param section the section name for the key
	 *  @param key the key name
	 *  @param default_value the value that is returned if the key is not found
	 *  @return the key value (if key found) otherwise default_value
	 */
	template<class T> T getValueSafe( std::string const& section, std::string const& key, T default_value ) const
	{
		T r;
		try{
		  r = getValue<T>( section, key );
		} catch( ErrItemNotFound& ) {
		  r = default_value;
		}
		return r;
	}

	
	//! exception safe version of getValue
	/*! returns the value of a given key, returns a default value rather
	 *  than a ErrItemNotFound exception if the key is not found.
	 *  @param key the key name
	 *  @param default_value the value that is returned if the key is not found
	 *  @return the key value (if key found) otherwise default_value
	 */
	template<class T> T getValueSafe( std::string const& key, T default_value ) const
	{
		return getValueSafe( "", key, default_value );
	}

	
	//! dumps all key-value pairs to a std::ostream
	void dump( std::ostream& out )
	{
		std::map<std::string,std::string>::const_iterator i = m_Items.begin();
		while( i!=m_Items.end() )
		{
			if( i->second.length() > 0 )
				out << std::setw(24) << std::left << i->first << "  =  " << i->second  << std::endl;
			++i;
		}
	}
	
	void log_dump( void )
	{
		LOGUSER("List of all configuration options:");
		std::map<std::string,std::string>::const_iterator i = m_Items.begin();
		while( i!=m_Items.end() )
		{
			if( i->second.length() > 0 )
				LOGUSER("  %24s = %s",(i->first).c_str(),(i->second).c_str());//out << std::setw(24) << std::left << i->first << "  =  " << i->second  << std::endl;
			++i;
		}
	}

  //--- EXCEPTIONS ---

	//! runtime error that is thrown if key is not found in getValue
	class ErrItemNotFound : public std::runtime_error{
	public:
		ErrItemNotFound( std::string itemname )
		  : std::runtime_error( itemname.c_str() )
		{}
	};

	//! runtime error that is thrown if type conversion fails
	class ErrInvalidConversion : public std::runtime_error{
	public:
		ErrInvalidConversion( std::string errmsg )
		  : std::runtime_error( errmsg )
		{}
	};


	//! runtime error that is thrown if identifier is not found in keys
	class ErrIllegalIdentifier : public std::runtime_error{
	public:
		ErrIllegalIdentifier( std::string errmsg )
		  : std::runtime_error( errmsg )
		{}
	};
  
};

//==== below are template specialisations =======================================//

//... Function: getValue( strSection, strEntry ) ...
//... Descript: specialization of getValue for type boolean to interpret strings ...
//...           like "true" and "false" etc.
//...           converts the string to type bool, returns type bool ...
template<> 
inline bool config_file::getValue<bool>( std::string const& strSection, std::string const& strEntry ) const{
  std::string r1 = getValue<std::string>( strSection, strEntry );
    if( r1=="true" || r1=="yes" || r1=="on" || r1=="1" )
      return true;
    if( r1=="false" || r1=="no" || r1=="off" || r1=="0" )
      return false;
    throw ErrIllegalIdentifier(std::string("Illegal identifier \'")+r1+std::string("\' in \'")+strEntry+std::string("\'."));
    //return false;
  }

template<>
inline bool config_file::getValueSafe<bool>( std::string const& strSection, std::string const& strEntry, bool defaultValue ) const{
  std::string r1;
  try{
    r1 = getValue<std::string>( strSection, strEntry );
    if( r1=="true" || r1=="yes" || r1=="on" || r1=="1" )
      return true;
    if( r1=="false" || r1=="no" || r1=="off" || r1=="0" )
      return false;
  } catch( ErrItemNotFound& ) {
    return defaultValue;
  }
  return defaultValue;
}

template<>
inline void config_file::convert<std::string,std::string>( const std::string & ival, std::string & oval) const
{
	oval = ival;
}

#endif //__CONFIG_FILE_HH
