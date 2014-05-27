#ifndef __LRVTRACK_FIT_HPP
#define __LRVTRACK_FIT_HPP
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/core/core.hpp>
#include <vector>
#include "cvblob.h"
#include "larvaObject.hpp"

//We direct our approximation initially 
//to the following spine points (p[0-11], 0:head 11:tail):
//  p[8]: main reference point
//  angleg: global angle of p[8] p[6]
//  angles: spine angles at p[2],p[4],p[6],p[8]
//
//  We need points p[0],p[2],p[4],p[6],p[8],p[11]
//  and their corresponding width points.
//
//  Angles are in degrees and anticlockwise, measured
//  from the old vector to the new.


using namespace cv;
using namespace cvb;
using namespace std;

extern size_t FRAME_COLS;
extern size_t FRAME_ROWS;
extern Mat greyFrame;

#define LRVFIT_ANGLES 4
#define IPOL 3

class larvaFit 
{
  private:
    size_t ID;
    size_t ipol=IPOL;
    double orig_width;
    double orig_length;
    std::vector<double> width_lookup;
    std::vector<Point2f> spine;
    std::vector<Point2f> intSpine;
    Point *cpoints;
    size_t cpoints_size;
    double degree_step=15*0.0174532925199432957692369076848;
    double wstep=0.02;
    vector<double> ga;
    vector<double> a1;
    vector<double> a2;
    vector<double> a3;
    vector<double> a4;
    vector<double> wl;
    void createContourFromFit(vector<Point2f> &F);
    void setupSpine();
    void pointToPoint(Point2f &p1, double angle, double d, Point2f &p2);
    //Mat completeFrame;
    void calculateContourPoints(Point2f &a,
        Point2f &b,
        Point2f &c,
        double b_index,
        double width,
        Point &cl,
        Point &cr);
    void calculateContourPoints(Point2f &a,
        Point2f &b,
        Point2f &c,
        Point2f &bp,
        double b_index,
        double width,
        Point &cl,
        Point &cr);
    void paintPoly(Mat &ROI, Point* f,size_t fsize);
  public:
    size_t PAD=10;
    unsigned long long ppdiffmax;
    unsigned long long createMatFromFitTime;
    size_t createMatFromFitCalls;
    unsigned long long paintPolyTime;
    size_t paintPolyCalls;
    unsigned long long setupSpineTime;
    size_t setupSpineCalls;
    unsigned long long calculateContourPointsTime;
    size_t calculateContourPointsCalls;
    void setloops(int ds,double ws);
    void setloops();
    struct fitData
    {
      Point2f midtail;
      double length;
      double width;
      vector<double> angles;
      double global_angle;


      fitData()
      {};
      fitData(Point2f mp,
              double l,
              double w,
              double a1,
              double a2,
              double a3,
              double a4,
              double ga):
                  midtail(mp),
                  length(l),
                  width(w),
                  global_angle(ga)
      {
        angles.push_back(a1);
        angles.push_back(a2);
        angles.push_back(a3);
        angles.push_back(a4);
      }
      fitData(Point2f mp,
              double l,
              double w,
              vector<double> a,
              double ga):
                  midtail(mp),
                  length(l),
                  width(w),
                  angles(a),
                  global_angle(ga){}
    };
    typedef struct fitData fitData;
    //std::vector<fitData> larvaFitData;
    fitData larvaFitData;
    void generate(std::vector<fitData> &l);
    ~larvaFit();
    larvaFit();
    larvaFit(const larvaFit &p);
    void setup(larvaObject &l,size_t FRAME=0);
    void setup(larvaObject &l,int dstep, double wstep, size_t FRAME=0);
    larvaFit(larvaObject &l,int dstep,double wstep, size_t FRAME=0);
    larvaFit(larvaObject &l,size_t FRAME=0);
    double errorFunc(Mat &b1C,
        size_t minx,
        size_t maxx,
        size_t miny,
        size_t maxy
        );
    double errorFunc(cvb::CvBlob &blob);
    void innerCountersReset();
    void showContour();
    double optimize2(cvb::CvBlob &blob);
    double optimize(cvb::CvBlob &blob);
    double optimize(vector<cvb::CvBlob> &vec, cvb::CvBlob &res);
    void createMatfromFit(Mat &larvaFitContour,
                                Mat &fitBase,
                                size_t minx,
                                size_t maxx,
                                size_t miny,
                                size_t maxy,
                                bool verbose
                                );
    void createMatfromFit(Mat &larvaFitContour,
                          size_t minx,
                          size_t maxx,
                          size_t miny,
                          size_t maxy,
                          bool vervose=false
                          );
    void filterAngle(std::vector<double> &a, 
        double &val,
        double lim,
        double add);
};


ostream& operator<<(ostream& cout,larvaFit &obj);

#endif
