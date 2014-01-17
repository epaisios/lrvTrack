#ifndef __LRVTRACK_DEBUG_HPP__
#define __LRVTRACK_DEBUG_HPP__
#include <string>

using namespace std;
//Debugging functions for output
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
