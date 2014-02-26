#ifndef __LRVTRACK_DEBUG_HPP__
#define __LRVTRACK_DEBUG_HPP__
#include <string>
#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>

using namespace std;
namespace logging = boost::log;
//
//Debugging functions for output
void log_init()
{
      logging::core::get()->set_filter
            (
                     logging::trivial::severity >= logging::trivial::debug
                         );
}
void verbosePrint(stringstream &toPrint)
{
  if(LRVTRACK_VERBOSE_LEVEL>0)
    {
      cout << "LrvTrack DEBUG: " << toPrint.str() << endl;
      toPrint.str("");
      toPrint.clear();
    }
}

void verbosePrint(const char * toPrint)
{
  if(LRVTRACK_VERBOSE_LEVEL>0)
    {
      cout << "LrvTrack DEBUG: " << toPrint << endl;
    }
}
void verbosePrint(string &toPrint)
{
  if(LRVTRACK_VERBOSE_LEVEL>0)
    {
      cout << "LrvTrack DEBUG: " << toPrint << endl;
    }
}

string printUIMap(map<unsigned int, unsigned int> &A)
{
  stringstream F;
  map<unsigned int, unsigned int>::iterator Ait=A.begin();
  F << "[ ";
  while(Ait!=A.end())
  {
    F << Ait->first << " -> " << Ait->second << ",";
    ++Ait;
  }
  F << " ]";
  return F.str();
}

#endif
