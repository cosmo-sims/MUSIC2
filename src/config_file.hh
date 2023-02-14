// This file is part of monofonIC (MUSIC2)
// A software package to generate ICs for cosmological simulations
// Copyright (C) 2020 by Oliver Hahn
// 
// monofonIC is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// monofonIC is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#pragma once

#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <typeinfo>

#include <logger.hh>

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
  unsigned iline_;

  //! hash table for key/value pairs, stored as strings
  std::map<std::string, std::string> items_;

public:
  //! removes all white space from string source
  /*!
   * @param source the string to be trimmed
   * @param delims a string of delimiting characters
   * @return trimmed string
   */
  std::string trim(std::string const &source,
                   char const *delims = " \t\r\n") const {
    std::string result(source);
    //... skip initial whitespace ...
    std::string::size_type index = result.find_last_not_of(delims);
    if (index != std::string::npos)
      result.erase(++index);
    //... find beginning of trailing whitespace ...
    index = result.find_first_not_of(delims);
    //... remove trailing whitespace ...
    if (index != std::string::npos)
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
  void convert(const in_value &ival, out_value &oval) const {
    std::stringstream ss;
    ss << ival; //.. insert value into stream
    ss >> oval; //.. retrieve value from stream

    if (!ss.eof()) {
      //.. conversion error
      music::elog << "Error: conversion of \'" << ival << "\' failed."
                << std::endl;
      throw except_invalid_conversion(std::string("invalid conversion to ") +
                                 typeid(out_value).name() + '.');
    }
  }

  //! constructor of class config_file
  /*! @param filename the path/name of the configuration file to be parsed
   */
  explicit config_file(std::string const &filename) : iline_(0), items_() {
    std::ifstream file(filename.c_str());

    if (!file.is_open()){
      music::elog << "Could not open config file \'" << filename << "\'." << std::endl;
      throw std::runtime_error(
          std::string("Error: Could not open config file \'") + filename +
          std::string("\'"));
    }

    std::string line;
    std::string name;
    std::string value;
    std::string in_section;
    int pos_equal;
    iline_ = 0;
    //.. walk through all lines ..
    while (std::getline(file, line)) {
      ++iline_;
      //.. encounterd EOL ?
      if (!line.length())
        continue;

      //.. encountered comment ?
      unsigned long idx;
      if ((idx = line.find_first_of("#;%")) != std::string::npos)
        line.erase(idx);

      //.. encountered section tag ?
      if (line[0] == '[') {
        in_section = trim(line.substr(1, line.find(']') - 1));
        continue;
      }

      //.. seek end of entry name ..
      pos_equal = line.find('=');
      name = trim(line.substr(0, pos_equal));
      value = trim(line.substr(pos_equal + 1));

      if ((size_t)pos_equal == std::string::npos &&
          (name.size() != 0 || value.size() != 0)) {
        music::wlog << "Ignoring non-assignment in " << filename << ":"
                  << iline_ << std::endl;
        continue;
      }

      if (name.length() == 0 && value.size() != 0) {
        music::wlog << "Ignoring assignment missing entry name in "
                  << filename << ":" << iline_ << std::endl;
        continue;
      }

      if (value.length() == 0 && name.size() != 0) {
        music::wlog << "Empty entry will be ignored in " << filename << ":"
                  << iline_ << std::endl;
        continue;
      }

      if (value.length() == 0 && name.size() == 0)
        continue;

      //.. add key/value pair to hash table ..
      if (items_.find(in_section + '/' + name) != items_.end()) {
        music::wlog << "Redeclaration overwrites previous value in "
                  << filename << ":" << iline_ << std::endl;
      }

      items_[in_section + '/' + name] = value;
    }
  }

  //! inserts a key/value pair in the hash map
  /*! @param key the key value, usually "section/key"
   *  @param value the value of the key, also a string
   */
  void insert_value(std::string const &key, std::string const &value) {
    items_[key] = value;
  }

  //! inserts a key/value pair in the hash map
  /*! @param section section name. values are stored under "section/key"
   *  @param key the key value usually "section/key"
   *  @param value the value of the key, also a string
   */
  void insert_value(std::string const &section, std::string const &key,
                   std::string const &value) {
    items_[section + '/' + key] = value;
  }

  //! checks if a key is part of the hash map
  /*! @param section the section name of the key
   *  @param key the key name to be checked
   *  @return true if the key is present, false otherwise
   */
  bool contains_key(std::string const &section, std::string const &key) {
    std::map<std::string, std::string>::const_iterator i =
        items_.find(section + '/' + key);
    if (i == items_.end())
      return false;
    return true;
  }

  //! checks if a key is part of the hash map
  /*! @param key the key name to be checked
   *  @return true if the key is present, false otherwise
   */
  bool contains_key(std::string const &key) {
    std::map<std::string, std::string>::const_iterator i = items_.find(key);
    if (i == items_.end())
      return false;
    return true;
  }

  //! return value of a key
  /*! returns the value of a given key, throws a except_item_not_found
   *  exception if the key is not available in the hash map.
   *  @param key the key name
   *  @return the value of the key
   *  @sa except_item_not_found
   */
  template <class T> T get_value(std::string const &key) const {
    return get_value<T>("", key);
  }

  //! return value of a key
  /*! returns the value of a given key, throws a except_item_not_found
   *  exception if the key is not available in the hash map.
   *  @param section the section name for the key
   *  @param key the key name
   *  @return the value of the key
   *  @sa except_item_not_found
   */
  template <class T>
  T get_value_basic(std::string const &section, std::string const &key) const {
    T r;
    std::map<std::string, std::string>::const_iterator i =
        items_.find(section + '/' + key);
    if (i == items_.end()){
      throw except_item_not_found('\'' + section + '/' + key +
                            std::string("\' not found."));
    }

    convert(i->second, r);
    return r;
  }

  template <class T>
  T get_value(std::string const &section, std::string const &key) const
  {
    T r;
    try
    {
      r = get_value_basic<T>(section, key);
    }
    catch (except_item_not_found& e)
    {
      music::elog << e.what() << std::endl;
      throw;
    }
    return r;
  }

  //! exception safe version of get_value
  /*! returns the value of a given key, returns a default value rather
   *  than a except_item_not_found exception if the key is not found.
   *  @param section the section name for the key
   *  @param key the key name
   *  @param default_value the value that is returned if the key is not found
   *  @return the key value (if key found) otherwise default_value
   */
  template <class T>
  T get_value_safe(std::string const &section, std::string const &key,
                 T default_value) const {
    T r;
    try {
      r = get_value_basic<T>(section, key);
    } catch (except_item_not_found&) {
      r = default_value;
      music::dlog << "Item \'" << section << "/" << key << " not found in config. Default = \'" << default_value << "\'" << std::endl;
    }
    return r;
  }

  //! exception safe version of get_value
  /*! returns the value of a given key, returns a default value rather
   *  than a except_item_not_found exception if the key is not found.
   *  @param key the key name
   *  @param default_value the value that is returned if the key is not found
   *  @return the key value (if key found) otherwise default_value
   */
  template <class T>
  T get_value_safe(std::string const &key, T default_value) const {
    return get_value_safe("", key, default_value);
  }

  //! dumps all key-value pairs to a std::ostream
  void dump(std::ostream &out) {
    std::map<std::string, std::string>::const_iterator i = items_.begin();
    while (i != items_.end()) {
      if (i->second.length() > 0)
        out << std::setw(24) << std::left << i->first << "  =  " << i->second
            << std::endl;
      ++i;
    }
  }

  void dump_to_log(void) {
    music::ulog << "List of all configuration options:" << std::endl;
    std::map<std::string, std::string>::const_iterator i = items_.begin();
    while (i != items_.end()) {
      if (i->second.length() > 0)
        music::ulog << std::setw(28) << i->first << " = " << i->second
                  << std::endl;
      ++i;
    }
  }

  //--- EXCEPTIONS ---

  //! runtime error that is thrown if key is not found in get_value
  class except_item_not_found : public std::runtime_error {
  public:
    except_item_not_found(std::string itemname)
        : std::runtime_error(itemname.c_str()) {}
  };

  //! runtime error that is thrown if type conversion fails
  class except_invalid_conversion : public std::runtime_error {
  public:
    except_invalid_conversion(std::string errmsg) : std::runtime_error(errmsg) {}
  };

  //! runtime error that is thrown if identifier is not found in keys
  class ErrIllegalIdentifier : public std::runtime_error {
  public:
    ErrIllegalIdentifier(std::string errmsg) : std::runtime_error(errmsg) {}
  };
};

//==== below are template specialisations
//=======================================//

//... Function: get_value( strSection, strEntry ) ...
//... Descript: specialization of get_value for type boolean to interpret strings
//...
//...           like "true" and "false" etc.
//...           converts the string to type bool, returns type bool ...
template <>
inline bool config_file::get_value<bool>(std::string const &strSection,
                                       std::string const &strEntry) const {
  std::string r1 = get_value<std::string>(strSection, strEntry);
  if (r1 == "true" || r1 == "yes" || r1 == "on" || r1 == "1")
    return true;
  if (r1 == "false" || r1 == "no" || r1 == "off" || r1 == "0")
    return false;
  music::elog << "Illegal identifier \'" << r1 << "\' in \'" << strEntry << "\'." << std::endl;
  throw ErrIllegalIdentifier(std::string("Illegal identifier \'") + r1 +
                             std::string("\' in \'") + strEntry +
                             std::string("\'."));
  // return false;
}

template <>
inline bool config_file::get_value_safe<bool>(std::string const &strSection,
                                           std::string const &strEntry,
                                           bool defaultValue) const {
  std::string r1;
  try {
    r1 = get_value_basic<std::string>(strSection, strEntry);
    std::transform(r1.begin(), r1.end(),r1.begin(), ::toupper);
    if (r1 == "YAY" || r1 == "TRUE" || r1 == "YES" || r1 == "ON" || r1 == "1")
      return true;
    if (r1 == "NAY" || r1 == "FALSE" || r1 == "NO" || r1 == "OFF" || r1 == "0")
      return false;
  } catch (except_item_not_found&) {
    return defaultValue;
  }
  return defaultValue;
}

template <>
inline void
config_file::convert<std::string, std::string>(const std::string &ival,
                                              std::string &oval) const {
  oval = ival;
}
