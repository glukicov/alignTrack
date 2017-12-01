#ifndef ALIGNTRACKER_H
#define ALIGNTRACKER_H

/* This header file contains definitions of auxiliary functionality 
*/

#include <iostream>
#include <fstream>
#include <string>

// overload the << operator to output to the terminal and log file simultaneously 
//https://stackoverflow.com/questions/24413744/using-operator-to-write-to-both-a-file-and-cout
struct MyStreamingHelper
{
    MyStreamingHelper(std::ostream& out1,
                      std::ostream& out2) : out1_(out1), out2_(out2) {}
    std::ostream& out1_;
    std::ostream& out2_;
};

template <typename T>
MyStreamingHelper& operator<<(MyStreamingHelper& h, T const& t)
{
   h.out1_ << t;
   h.out2_ << t;
   return h;
}

MyStreamingHelper& operator<<(MyStreamingHelper& h, std::ostream&(*f)(std::ostream&))
{
   h.out1_ << f;
   h.out2_ << f;
   return h;
}


#endif