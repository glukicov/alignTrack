#ifndef DATETIME_HH
#define DATETIME_HH 1

// Simple utility class for manipulating, printing times and dates

#include <string>
#include <vector>
#include <time.h>
#include <iostream>
#include <cstring>

using std::istream;
using std::ostream;
using std::endl;
using std::cout;
using std::string;

class DateTime
{
public:
        DateTime() ;
        DateTime(const DateTime&);
        DateTime(int time);       
        DateTime(const char* timeStr);
        ~DateTime() {} ;

	int value() const { return _value; }

        int CurrentTime()  ;
        void SetTime (int time) { _value = time;}

        void SetDate(int year,int month,int day);

        int Year() ;
        int Year(int time) ;
        int Month() ;
        int Month(int time) ;
        int Day() ;
        int Day(int time) ;
        int Hour() ;
        int Hour(int time) ;
        int Minute() ;
        int Minute(int time) ;
        int Second() ;
        int Second(int time) ;

        std::string TimeFull(int time);
        std::string TimeFull();
        std::string Format(int option=0);

        operator int() const { return _value; }
        DateTime& operator=(const DateTime& d);

        friend ostream& operator << (std::ostream& ost, const DateTime& t) {
	  time_t    now;  
	  now = (time_t) t._value ;
	  ost << ctime(&now) ;
          return ost;
	}

        friend istream& operator >> (istream& ist, DateTime& t) {
	  time_t    now;  
	  now = (time_t) t._value ;
	  if (!(ist >> ctime(&now))) return ist;
	  return ist;
	}

private:
       int _value;


};

inline DateTime::DateTime(const DateTime& d){
  _value = d._value;
}

inline DateTime& DateTime::operator=(const DateTime& d){
  _value=d._value ;
  return *this;
}

inline DateTime::DateTime( ){
   _value = CurrentTime() ; 
}

inline DateTime::DateTime(int time): _value(time) { } 

inline void DateTime::SetDate(int year,int month,int day)  { 
       
  struct tm tmx; ;
  memset(&tmx, 0, sizeof(tmx));
  if (year >= 0 && year < 1900) year = 2000 + year ; 
  tmx.tm_year   = year - 1900;
  tmx.tm_mon    = month - 1 ;
  tmx.tm_mday   = day;
  
  _value  = mktime(&tmx);
  
}


inline int DateTime::CurrentTime()  { 

  time_t    now;  now = time((time_t *)NULL);
  return  (int) now; 
  
}
inline std::string DateTime::TimeFull()     { 
  time_t    now;  
  now = (time_t) _value ;
  return ctime(&now) ;
  }

inline std::string DateTime::TimeFull(int time)    { 
  time_t    now;  
  now = (time_t) time ;
  return ctime(&now) ;
}

inline int DateTime::Year()     { 
  time_t    now;  	struct tm *l_time;

  now    = (time_t) _value ;
  l_time = localtime(&now);           
  
  return 1900+l_time->tm_year ;
  
}

inline int DateTime::Year(int time)     { 
  time_t    now;  	struct tm *l_time;
  
  now    = (time_t) time ;
  l_time = localtime(&now);           
  
  return 1900+l_time->tm_year ;
  
}


inline int DateTime::Month()     { 
  time_t    now;  	struct tm *l_time;
  
  now    = (time_t) _value ;
  l_time = localtime(&now);           
  
  return l_time->tm_mon + 1 ;
  
}


inline int DateTime::Month(int time)     { 
  time_t    now;  	struct tm *l_time;
  
  now    = (time_t) time ;
  l_time = localtime(&now);           
  
  return l_time->tm_mon + 1;
  
}


inline int DateTime::Day()     { 
  time_t    now;  	struct tm *l_time;
  
  now    = (time_t) _value ;
  l_time = localtime(&now);           
  
  return l_time->tm_mday ;

}

inline int DateTime::Day(int time)     { 
  time_t    now;  	struct tm *l_time;

  now    = (time_t) time ;
  l_time = localtime(&now);           
                              
  return l_time->tm_mday ;

}

inline int DateTime::Hour()     { 
  time_t    now;  	struct tm *l_time;

  now    = (time_t) _value ;
  l_time = localtime(&now);           
                              
  return l_time->tm_hour ;

}

inline int DateTime::Hour(int time)     { 
  time_t    now;  	struct tm *l_time;
  
  now    = (time_t) time ;
  l_time = localtime(&now);           
  
  return l_time->tm_hour ;
  
}

inline int DateTime::Minute()     { 
  time_t    now;  	struct tm *l_time;
  
  now    = (time_t) _value ;
  l_time = localtime(&now);           
                              
  return l_time->tm_min ;
}

inline int DateTime::Minute(int time)     { 
  time_t    now;  	struct tm *l_time;

  now    = (time_t) time ;
  l_time = localtime(&now);           
                              
  return l_time->tm_min ;

}

inline int DateTime::Second()     { 
  time_t    now;  	struct tm *l_time;

  now    = (time_t) _value ;
  l_time = localtime(&now);           
  
  return l_time->tm_sec ;

}

inline int DateTime::Second(int time)     { 
  time_t    now;  	struct tm *l_time;

  now    = (time_t) time ;
  l_time = localtime(&now);           
  
  return l_time->tm_sec ;
	
}

inline string DateTime::Format(int option)     { 
  char test [50];
  std::string fmt;
  time_t now;       
  time(&now);       
  

  if (option == 0) {
    fmt = "-%Y-%m-%d_%H_%M_%S";
  }
  else if (option == 1) {
    fmt = "%Y-%m-%d_%H_%M_%S";
  }
  else if (option == 2) {
    fmt = "-%Y-%m-%d";
  }
  else if (option == 3) {
    fmt = "%Y-%m-%d";
  }
  else if (option == 4) {
    fmt = "-%H_%M_%S";
  }
  else if (option == 5) {
    fmt = "%H:%M:%S";
  }
  else if (option == 6) {
    fmt = "%m/%d/%Y %H:%M:%S";
  }

  strftime(test,sizeof test, fmt.c_str(), localtime(&now));
  return string(test);
  
}

#endif
