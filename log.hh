/*
 
 log.hh - This file is part of MUSIC -
 a code to generate multi-scale initial conditions 
 for cosmological simulations 
 
 Copyright (C) 2010  Oliver Hahn
 
 */

#ifndef __LOG_HH
#define __LOG_HH

#include <string>
#include <list>
#include <fstream>
#include <ctime>
#include <cstdarg>
#include <sstream>

/*!
 *	\brief	System for logging runtime library errors, warnings, etc.
 *
 *	This is the class that catches every (debug) info, warning, error, or user message and
 *	processes it. Messages can be written to files and/or forwarded to user function for
 *	processing messages.
 */
namespace MUSIC
{
	
class log
{
public:
	log(){}
	~log();
	
	/*!
	 *	\brief	Types of logged messages.
	 */
	enum messageType
	{
		Info,
		DebugInfo,
		Warning,
		Error,
		FatalError,
		User
	};
	
	/*!
	 *	\brief	Logged message of type MessageType with some info.
	 */
	struct message
	{
		messageType type;
		std::string text;
		tm* when;
	};
	
	/*!
	 *	\brief	Open file where to log the messages.
	 */
	static void setOutput(const std::string& filename);
	
	/*!
	 *	\brief	Get the filename of log.
	 */
	static const std::string& output() { return outputFile_; }
	
	/*!
	 *	\brief	Add a new message to log.
	 *	\param	type	Type of the new message.
	 *	\param	text	Message.
	 *	\remarks Message is directly passes to user reciever if one is set.
	 */
	static void send(messageType type, const std::string& text);
	//static void send(messageType type, std::string& text);
	
	/*!
	 *	\brief	Get the list of all of the logged messages.
	 */
	static const std::list<message>& messages() { return messages_; }
	
	/*!
	 *	\brief	Get the last logged message.
	 */
	static const message& lastMessage() { return messages_.back(); }
	
	/*!
	 *	\brief	Set user function to receive newly sent messages to logger.
	 */
	static void setUserReceiver(void (*userFunc)(const message&)) { receiver = userFunc; }
	
	/*!
	 *	\brief	Set minimum level of message to be logged.
	 */
	static void setLevel(const log::messageType level);
	
private:
	
	static std::string outputFile_;
	static std::ofstream outputStream_;
	static std::list<message> messages_;
	static messageType logLevel_;
	static void (*receiver)(const message&);
};

}


inline void LOGERR( const char* str, ... )
{
	char out[1024];
	va_list argptr;
	va_start(argptr,str);
	va_end(argptr);
	vsprintf(out,str,argptr);
	MUSIC::log::send(MUSIC::log::Error, std::string(out));
}

inline void LOGWARN( const char* str, ... )
{
	char out[1024];
	va_list argptr;
	va_start(argptr,str);
	va_end(argptr);
	vsprintf(out,str,argptr);
	MUSIC::log::send(MUSIC::log::Warning, std::string(out));
}

inline void LOGFATAL( const char* str, ... )
{
	char out[1024];
	va_list argptr;
	va_start(argptr,str);
	va_end(argptr);
	vsprintf(out,str,argptr);
	MUSIC::log::send(MUSIC::log::FatalError, std::string(out));
}

inline void LOGDEBUG( const char* str, ... )
{
	char out[1024];
	va_list argptr;
	va_start(argptr,str);
	va_end(argptr);
	vsprintf(out,str,argptr);
	MUSIC::log::send(MUSIC::log::DebugInfo, std::string(out));
}

inline void LOGUSER( const char* str, ... )
{
	char out[1024];
	va_list argptr;
	va_start(argptr,str);
	va_end(argptr);
	vsprintf(out,str,argptr);
	MUSIC::log::send(MUSIC::log::User, std::string(out));
}

inline void LOGINFO( const char* str, ... )
{
	char out[1024];
	va_list argptr;
	va_start(argptr,str);
	va_end(argptr);
	vsprintf(out,str,argptr);
	MUSIC::log::send(MUSIC::log::Info, std::string(out));
}



/*#ifndef LOGERR
#define LOGERR(x) { std::stringstream ss; ss<<(x); MUSIC::log::send(MUSIC::log::Error, ss); }
#endif

#ifndef LOGINFO
#define LOGINFO(x) { std::stringstream ss; ss<<(x); MUSIC::log::send(MUSIC::log::Info, ss); }
#endif

#ifndef LOGWARN
#define LOGWARN(x) { std::stringstream ss; ss<<(x); MUSIC::log::send(MUSIC::log::Warning, ss); }
#endif

#ifndef LOGFATAL
#define LOGFATAL(x) { std::stringstream ss; ss<<(x); MUSIC::log::send(MUSIC::log::FatalError, ss); }
#endif

#ifndef LOGDEBUG
#define LOGDEBUG(x) { std::stringstream ss; ss<<(x); MUSIC::log::send(MUSIC::log::DebugInfo, ss); }
#endif

#ifndef LOGUSER
#define LOGUSER(x) { std::stringstream ss; ss<<(x); MUSIC::log::send(MUSIC::log::User, ss); }
#endif*/

/*#define LOGINFO(x) MUSIC::log::send(MUSIC::log::Info, std::stringstream()<<(x));
#define LOGWARN(x) MUSIC::log::send(MUSIC::log::Warning, std::stringstream()<<(x));
#define LOGFATAL(x) MUSIC::log::send(MUSIC::log::FatalError, std::stringstream()<<(x));
*/


#endif //__LOG_HH


