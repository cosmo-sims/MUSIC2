/*
 
 log.cc - This file is part of MUSIC -
 a code to generate multi-scale initial conditions 
 for cosmological simulations 
 
 Copyright (C) 2010  Oliver Hahn
 
 */

#include "log.hh"
#include <iostream>
#include <algorithm>


std::string MUSIC::log::outputFile_;
std::ofstream MUSIC::log::outputStream_;
std::list<MUSIC::log::message> MUSIC::log::messages_;
void (*MUSIC::log::receiver)(const message&) = NULL;
MUSIC::log::messageType MUSIC::log::logLevel_;


std::string RemoveMultipleWhiteSpaces( std::string s )
{ 
	std::string search = "  "; // this is 2 spaces
	size_t index;
	
	while( (index = s.find(search)) != std::string::npos )
	{ // remove 1 character from the string at index
		s.erase(index,1);
	}
	
	return s;
}

void MUSIC::log::send(messageType type, const std::string& text_)
//void MUSIC::log::send(messageType type, std::stringstream& textstr)
{
	std::string text(text_);// = textstr.str();
	// Skip logging if minimum level is higher
    if (logLevel_)
		if (type < logLevel_) return;
	// log message
	MUSIC::log::message m;
	m.type = type;
	m.text = text;
	time_t t = time(NULL);
	m.when = localtime(&t);
	messages_.push_back(m);
	
	if( type==Info||type==Warning||type==Error||type==FatalError )
	{	
		std::cout << " - "; 
		if(type==Warning)
			std::cout << "WARNING: ";
		if(type==Error)
			std::cout << "ERROR: ";
		if(type==FatalError)
			std::cout << "FATAL: ";
		std::cout << text << std::endl;
	}
	
	std::replace(text.begin(),text.end(),'\n',' ');
	RemoveMultipleWhiteSpaces(text);
	
	// if enabled logging to file
	if(outputStream_.is_open())
	{
		// print time
		char buffer[9];
		strftime(buffer, 9, "%X", m.when);
		outputStream_ << buffer;
		
		// print type
		switch(type)
		{
			case Info:		outputStream_ << " |  info   | "; break;
			case DebugInfo: outputStream_ << " |  debug  | "; break;
			case Warning:	outputStream_ << " | warning | "; break;
			case Error:		outputStream_ << " |  ERROR  | "; break;
			case FatalError:outputStream_ << " |  FATAL  | "; break;
			case User:		outputStream_ << " |  info   | "; break;
			default:		outputStream_ << " | ";
		}
		
		// print description
		outputStream_ << text << std::endl;
	}
	
	// if user wants to catch messages, send it to him
	if(receiver)
		receiver(m);
}


void MUSIC::log::setOutput(const std::string& filename)
{
	//logDebug("Setting output log file: " + filename);
	outputFile_ = filename;
	
	// close old one
	if(outputStream_.is_open())
		outputStream_.close();
	
	// create file
	outputStream_.open(filename.c_str());
	if(!outputStream_.is_open())
		LOGERR("Cannot create/open \'%s\' for logging",filename.c_str());
		//send(Error, "Cannot create/open '" + filename + "' for logging");
}

void MUSIC::log::setLevel(const MUSIC::log::messageType level) 
{
    logLevel_ = level;
}


MUSIC::log::~log()
{
	if(outputStream_.is_open())
		outputStream_.close();
}

