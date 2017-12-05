#ifndef Logger_hh
#define Logger_hh

#include <sstream>
#include <string>
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <mutex>

#include "Destroyer.hh" 
#include "DateTime.hh" 

using namespace std;


//Exception class thrown when Logger::ERROR written
class CriticalError : public std::exception {
public:
  CriticalError() : std::exception() {}
  virtual const char* what() const throw() { return "CRITICAL ERROR"; }
};


//stdout color handling
//Stolen from http://stackoverflow.com/questions/2616906/how-do-i-output-coloured-text-to-a-linux-terminal
#include <ostream>
namespace Color {
  enum Code {
    FG_RED      = 31,
    FG_GREEN    = 32,
    FG_YELLOW   = 33,
    FG_BLUE     = 34,
    FG_DEFAULT  = 39,
    BG_DEFAULT  = 49
    };
  class Modifier {
    Code code;
  public:
     Modifier(Code pCode) : code(pCode) {}
     void setCode(Code pCode) { code = pCode; }
     friend std::ostream& operator<<(std::ostream& os, const Modifier& mod) {
       return os << "\033[" << mod.code << "m";
     }
  };
}


//Logger singleton class
class Logger {
public:
  static Logger* Instance();
  static Logger* Instance(const unsigned int logLevel_);

  //Define various levels of debugging/print-out (less to more)
  static unsigned int const ERROR                    = 0;
  static unsigned int const INFO                     = 1;
  static unsigned int const WARNING                  = 2;
  static unsigned int const NOTE                     = 3;
  static unsigned int const DEBUG                    = 4;
  
  void setLogLevel(unsigned int logLevel_);
  unsigned int getLogLevel();
  void setStylePlain();
  void setStyle(unsigned int style_);

  void LogToFile(string fileName_);
  void CloseLogFile();

  void enableCriticalErrorThrow() { _throwCriticalErrors = true; }
  void disableCriticalErrorThrow() { _throwCriticalErrors = false; }

  void setUseColor(bool use_) { _useColor = use_; }

  // Write a message with a given log level (can pass custom stream if desired)
  void write(unsigned int level_, const std::string & message_);

  std::string getMostRecentMsg() const;

  //External access to color modifiers
  static Color::Modifier& red() { return _red; }
  static Color::Modifier& yellow() { return _yel; }
  static Color::Modifier& green() { return _green; }
  static Color::Modifier& blue() { return _blue; }
  static Color::Modifier& def() { return _def; }

  //static void setOutstream(const std::ostream & o ) { _outstream = o; }

protected:
  Logger() { }
  static void constructor(const unsigned int logLevel_);
  void write(unsigned int level_, const std::string & message_, const std::lock_guard<std::mutex> & lock_);

  friend class Destroyer<Logger>;
  virtual ~Logger() { }

private:
  static Logger* _instance;
  static Destroyer<Logger> _destroyer;
  static unsigned int _logLevel;
  static unsigned int _style;
  static DateTime _time;

  static std::mutex _mutex; //For making logging thread-safe

  static string _opFileName; //Name and path of file to output to
  static ofstream _opFile;
  static bool _logToFile;

  static bool _throwCriticalErrors;

  static std::string _prevMsg; //Most recent message sent

  //Color handling
  static bool _useColor;
  static Color::Modifier _blue;
  static Color::Modifier _red;
  static Color::Modifier _yel;
  static Color::Modifier _green;
  static Color::Modifier _def;

};


#endif
