#include "lrvTrackDebug.hpp"

namespace logging = boost::log;
void log_init()
{
  logging::core::get()->set_filter
    (
     logging::trivial::severity >= logging::trivial::debug
    );
}
void verbosePrint(std::stringstream &toPrint)
{
  if(LRVTRACK_VERBOSE_LEVEL>0)
    {
      std::cout << "LrvTrack DEBUG: " << toPrint.str() << std::endl;
      toPrint.str("");
      toPrint.clear();
    }
}

void verbosePrint(const char * toPrint)
{
  if(LRVTRACK_VERBOSE_LEVEL>0)
    {
      std::cout << "LrvTrack DEBUG: " << toPrint << std::endl;
    }
}
void verbosePrint(std::string &toPrint)
{
  if(LRVTRACK_VERBOSE_LEVEL>0)
    {
      std::cout << "LrvTrack DEBUG: " << toPrint << std::endl;
    }
}

std::string printUIMap(std::map<unsigned int, unsigned int> &A)
{
  std::stringstream F;
  std::map<unsigned int, unsigned int>::iterator Ait=A.begin();
  F << "[ ";
  while(Ait!=A.end())
  {
    F << Ait->first << " -> " << Ait->second << ",";
    ++Ait;
  }
  F << " ]";
  return F.str();
}

