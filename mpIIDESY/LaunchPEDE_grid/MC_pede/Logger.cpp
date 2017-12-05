#include <stdlib.h> 

#include "Logger.hh"
#include "DateTime.hh"

// Must constructorialise the singletons/static variables here
Logger* Logger::_instance = 0;
unsigned int Logger::_logLevel = 0;
unsigned int Logger::_style = 0;
DateTime Logger::_time;
std::mutex Logger::_mutex;
string Logger::_opFileName;
ofstream Logger::_opFile;
bool Logger::_logToFile;
bool Logger::_throwCriticalErrors;
bool Logger::_useColor;
Color::Modifier Logger::_red(Color::FG_RED);
Color::Modifier Logger::_yel(Color::FG_YELLOW);
Color::Modifier Logger::_green(Color::FG_GREEN);
Color::Modifier Logger::_blue(Color::FG_BLUE);
Color::Modifier Logger::_def(Color::FG_DEFAULT);
std::string Logger::_prevMsg;

// Define template
Destroyer<Logger> Logger::_destroyer;

Logger* Logger::Instance () {
  std::lock_guard<std::mutex> lock(_mutex); //Lock mutex for duration of function
  constructor(Logger::INFO); //Default log level
  return _instance;
}


Logger* Logger::Instance (const unsigned int logLevel_) {
  std::lock_guard<std::mutex> lock(_mutex); //Lock mutex for duration of function
  constructor(logLevel_);
  return _instance;
}


//Common function called by the overloaded "Instance" functions
void Logger::constructor (const unsigned int logLevel_) {

  if (!_instance) {
    _instance = new Logger;
    _destroyer.SetDoomed(_instance);
    _logLevel = logLevel_;
    _logToFile = false;
    _style = 0;;
   _throwCriticalErrors = true;
   _useColor = true;
  }

  return;
}


void Logger::write(unsigned int level_, const std::string & message_) {

  //Write doesn't lock mutex itself by default in case it is called by
  //other functions that have already locked the mutex.
  //This public overload does lock it.
  std::lock_guard<std::mutex> lock(_mutex);
  this->write(level_,message_,lock);

}


void Logger::write(unsigned int level_, const std::string & message_, const std::lock_guard<std::mutex> & lock_) {

  // Only print if level is above the set log level threshold
  if ( level_ <= _logLevel ) {

    std::stringstream msg;

    string levelStr = "UNKNOWN";
    if(level_ == Logger::ERROR) levelStr = "ERROR";
    else if(level_ == Logger::INFO) levelStr = "INFO";
    else if(level_ == Logger::WARNING) levelStr = "WARNING";
    else if(level_ == Logger::NOTE) levelStr = "NOTE";
    else if(level_ == Logger::DEBUG) levelStr = "DEBUG";
    else {
      std::stringstream ss; ss << "[ Logger :: write ] Invalid log level " << level_ << " chosen, must be ERROR, INFO, WARNING, NOTE or DEBUG";
      this->write(Logger::ERROR,ss.str(),lock_);
      return;
    }

    //Store as new "most recent" message
    _prevMsg = message_;

    //Compose message based on style
    if (_style > 2) _style = 0; //Default to 0 if invalid style set
    if (_style == 0) {
      msg << message_; //Message only
    }
    else if (_style == 1)  {
      msg << levelStr << " : " << message_; //Add message level
    }
    else if (_style == 2)  {
      string date = _time.Format() ;
      msg << levelStr << " : " << date << " : " << message_; //Add message level and date
    }

    //Send msg to stdout
    if( _useColor ) { //Add color case
      std::stringstream colorMsg;
      if( level_ == Logger::ERROR ) colorMsg << _red;
      else if( level_ == Logger::WARNING ) colorMsg << _yel;
      else colorMsg << _def;
      colorMsg << msg.str() << _def << std::endl;
      std::cout << colorMsg.str();
    }
    else { //No color case
      std::cout << msg.str() << std::endl;
    }
    
    //Send same message to file (no color)
    if (_logToFile) {
      _opFile << msg.str() << std::endl;
      _opFile.flush(); //Write the stream to the file
    }

    //Throw CriticalError if required
    if(_throwCriticalErrors) {
      if( level_ == Logger::ERROR ) throw CriticalError();
    }

  }//if within debug level

}

unsigned int Logger::getLogLevel() {
  std::lock_guard<std::mutex> lock(_mutex); //Lock mutex for duration of function
  return _logLevel;
}

void Logger::setLogLevel(const unsigned int logLevel) {

  std::lock_guard<std::mutex> lock(_mutex); //Lock mutex for duration of function

  //Check for valid debug level and set if OK
  if( logLevel > 4 ) {
    std::stringstream ss; ss << "[ Logger :: setLogLevel ] Invalid log level " << logLevel << " chosen, must be ERROR, INFO, WARNING, NOTE or DEBUG";
    this->write(Logger::ERROR,ss.str(),lock);
    _logLevel = 0;
  }
  else {
    _logLevel = logLevel;
  }


}

void Logger::setStylePlain() {
  this->setStyle(0);
}

void Logger::setStyle(unsigned int style_) {

  std::lock_guard<std::mutex> lock(_mutex); //Lock mutex for duration of function

  //Check for valid style level and set if OK
  if( style_ > 2 ) {
    std::stringstream ss; ss << "[ Logger :: setStyle ] Invalid style " << style_ << " chosen, must be in range 0-2";
    this->write(Logger::ERROR,ss.str(),lock);
    _style = 0;
  }
  else {
    _style = style_;
  }

}

//Enable logging to file (as well as stdout)
void Logger::LogToFile(string fileName_) {

  std::lock_guard<std::mutex> lock(_mutex); //Lock mutex for duration of function

  //Check if already logging to file
  if(_logToFile) {
    stringstream ss; ss << "[ Logger :: LogToFile ] Already logging to file";
    this->write(Logger::ERROR,ss.str(),lock);
    return;
  }

  //Check if stream already exists
  //TODO

  //Open file
  _opFileName = fileName_;
  _opFile.open(_opFileName.c_str()); //TODO Catch exceptions

  //If open was successful, set flag
  if(_opFile) {
    _logToFile = true;
    stringstream ss; ss << "[ Logger :: LogToFile ] Logging to file: " << _opFileName;
    this->write(Logger::INFO,ss.str(),lock);
  }
  else {
    stringstream ss; ss << "[ Logger :: LogToFile ] Failed to open file: " << _opFileName << std::endl;
    this->write(Logger::ERROR,ss.str(),lock);
  }

}


//Write the file stream to the file
void Logger::CloseLogFile() {

  std::lock_guard<std::mutex> lock(_mutex); //Lock mutex for duration of function

  //Check that we are logging to file
  if(_logToFile) {

    _logToFile = false;

    stringstream ss; ss << "[ Logger :: CloseLogFile ] Closing log file: " << _opFileName;
    write(Logger::INFO,ss.str(),lock);

    _opFile.close();

  }
  else {
    stringstream ss; ss << "[ Logger :: CloseLogFile ] Cannot close log file as not currently logging to file";
    write(Logger::ERROR,ss.str(),lock);
    return;
  }

}


std::string Logger::getMostRecentMsg() const { 
  std::lock_guard<std::mutex> lock(_mutex); //Lock mutex for duration of function
  return _prevMsg;
}

