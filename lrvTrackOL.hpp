#ifndef __LRVTRACKOL_HPP__
#define __LRVTRACKOL_HPP__
#include <string>
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/timer/timer.hpp>
#include <iomanip>
#include <numeric>
#include <string>
#include "lrvTrack.hpp"
#include "cvblob.h"
#include "lrvTrackBase.hpp"
#include "blobUtils.hpp"
#include "larvaDistanceMap.hpp"
#include "larvaSkel.hpp"
#include "lrvTrackDebug.hpp"
#include <fstream>
#ifndef _WIN32
#include <sys/time.h>
#endif

namespace po = boost::program_options;
namespace fs = boost::filesystem;
using namespace cv;
using namespace std;

void updateOneLarva(cvb::CvBlobs &In,
                    cvb::CvBlobs &Prev,
                    cvb::CvBlobs::iterator it,
                    tbb::concurrent_hash_map<unsigned int, larvaObject> &NEW_LARVA);

double mh_dist(unsigned int N,unsigned int C);
double kn_dist(unsigned int N,unsigned int C);

/* 
 * Class to perform the computations for each larva.
 * The class is needed to perform these operations in parallel.
 * The format is based on this answer:
 * http://answers.opencv.org/question/3730/how-to-use-parallel_for/
 */
class larvaeUpdateBody : public ParallelLoopBody
{
private:
  cvb::CvBlobs &In;
  cvb::CvBlobs &Prev;
  tbb::concurrent_hash_map<unsigned int, larvaObject> &NEW_LARVA;
  cvb::CvBlobs::iterator it;

public:
  larvaeUpdateBody(cvb::CvBlobs &IIn, 
      cvb::CvBlobs &IPrev,
      tbb::concurrent_hash_map<unsigned int, larvaObject> &n
      ): 
    In(IIn),
    Prev(IPrev),
    NEW_LARVA(n)
  {}
    void operator ()(const Range& range) const
    {
        for (int i = range.start; i < range.end; ++i)
        {
          cvb::CvBlobs::iterator it=In.begin();
          advance(it,i);
            updateOneLarva(In,Prev,it,NEW_LARVA);
        }
    }
    void operator ()(const BlockedRange& range) const
    {
        for (int i = range.begin(); i < range.end(); ++i)
        {
          cvb::CvBlobs::iterator it=In.begin();
          advance(it,i);
            updateOneLarva(In,Prev,it,NEW_LARVA);
        }
    }
};


/*
 * class to encode the possible mappings for the larvae before and after a collision.
 */
class lrvMapping {

    double dst;
  
  public:
    //The actual mapping
    // newID -> detectedLarvaeID
    pair<unsigned int,unsigned int> mapping;
    //vector<unsigned int> candidates;
    unsigned int nlrv;
    unsigned int plrv;
    //Constructor function sets up the mapping
    lrvMapping(unsigned int a,unsigned int b)
    {
      mapping=make_pair<unsigned int,unsigned int>(a,b);
      nlrv=mapping.first;
      plrv=mapping.second;
      dst=-1;
    }
    void setDistance()
    {
       //dst=kn_dist(nlrv,plrv);
       dst=mh_dist(nlrv,plrv);
    }
    void print()
    {
      cerr << "M[" << nlrv << "->" << plrv << "]" << " #" << dst << endl;
    }
    double getDistance()
    {
      if(dst==-1)
        dst=mh_dist(nlrv,plrv);
        //dst=kn_dist(nlrv,plrv);
      
      return dst;
    }
};

// Vector of mappings: covers the possible assignments.
typedef vector<lrvMapping> ltAssignments;

void pair_powersets(vector<lrvMapping> &IN, 
    vector<ltAssignments > &OUT);


#endif
