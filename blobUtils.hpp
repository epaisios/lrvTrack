#ifndef __LRVTRACK_BLOBUTILS_HPP__
#define __LRVTRACK_BLOBUTILS_HPP__
#include <opencv2/core/core.hpp>
#include "cvblob.h"
#define ROI_PADDING 0

unsigned long long getmsofday();

void spline3(std::vector<cv::Point2f> &cp,
           std::vector<float> &d,
           std::vector<float> &w,
           int n,
           int RES,
           std::vector<cv::Point2f> &np,
           std::vector<size_t> &di);

void spline4(std::vector<cv::Point2f> &cp,
           std::vector<float> &d,
           std::vector<float> &w,
           int n,
           int RES,
           std::vector<cv::Point2f> &np,
           std::vector<size_t> &di,
           std::map<float,size_t> &cm,
           std::vector<float> &vcurv,
           double &cvariance);

void spline2(std::vector<cv::Point2f> &cp,
           std::vector<float> &d,
           std::vector<float> &w,
           int n,
           int RES,
           std::vector<cv::Point2f> &np,
           std::vector<size_t> &di,
           std::map<float,size_t> &cm,
           std::vector<float> &vcurv);

double p2fdist(double x1,double y1, double x2, double y2);

double p2fdist(cv::Point2f &a, cv::Point2f &b);

double diff(cv::Point2f &a, cv::Point2f &b);
void blobToPointVector(cvb::CvBlob &p,std::vector<cv::Point2f> &points,size_t PAD=0);

double vecVariance(std::vector<double> &c);

void pointsToContourVector(cvb::CvBlob &blob,
                         std::vector<cv::Point2f> &kpoints,
                         cv::Mat &frame, 
                         int PAD,
                         std::vector<cv::Point2f> &points);

void blobToContourVector(cvb::CvBlob &p,
                         cv::Mat &frame, 
                         int PAD,
                         std::vector<cv::Point2f> &points);

void createLarvaROI(cv::Mat &frame, cv::Mat &ROI, cvb::CvBlob &blob);

void createLarvaContour(cv::Mat &lrvROI,
                        cvb::CvBlob &blob,
                        int type=CV_8UC1,
                        int PADDING=0,
                        bool FILL=true,
                        cv::Scalar color=cv::Scalar(255),
                        int connectivity=4,
                        cv::Scalar bg=cv::Scalar(0)
                        );

void createLarvaContour_custom(cv::Mat &lrvROI,
                        cvb::CvBlob &blob,
                        int type=CV_8UC1,
                        int minx=0,
                        int maxx=0,
                        int miny=0,
                        int maxy=0,
                        int PAD=0,
                        bool FILL=true,
                        cv::Scalar color=cv::Scalar(255),
                        int connectivity=4,
                        cv::Scalar bg=cv::Scalar(0)
    );;

/*
void createLarvaContourCV(cv::Mat &lrvROI,
                        cvb::CvBlob &blob,
                        int type=CV_8UC1,
                        int PADDING=0,
                        bool FILL=true);

*/
cv::Point2f px3Smooth(cv::Mat &f, cv::Point2f &e, cv::Point2f &a, cv::Point2f &b, cv::Point2f &c, cv::Point2f &d);
cv::Point2f px5Smooth(cv::Mat &f, cv::Point2f &a, cv::Point2f &b, cv::Point2f &c, cv::Point2f &d, cv::Point2f &e);

void smoothPoints(std::vector<cv::Point2f> &origPoints, std::vector<cv::Point2f> &smoothened,cv::Mat &frame, size_t rep);

void spline(std::vector<cv::Point2f> &cp,
           std::vector<float> &d,
           std::vector<float> &w, 
           int n, 
           float xp1, 
           float yp1, 
           float xpn, 
           float ypn, 
           std::vector<float> &x2, 
           std::vector<float> &y2);

void splint(std::vector<cv::Point2f> &cp,
            std::vector<float> &x2a, 
            std::vector<float> &y2a, 
            std::vector<float> &d, 
            int n, 
            float t, 
            int khi,
            int klo,
            cv::Point2f &np);

/*void sspline(std::vector<float> &cp,
           std::vector<float> &t,
           std::vector<float> &w,
           int n, 
           float s,
           std::vector<float> &ax, 
           std::vector<float> &bx, 
           std::vector<float> &cx, 
           std::vector<float> &dx);

void ssplint(std::vector<float> &d,
           float t,
           int klo,
           std::vector<float> &ax, 
           std::vector<float> &bx, 
           std::vector<float> &cx, 
           std::vector<float> &dx,
           std::vector<float> &ay, 
           std::vector<float> &by, 
           std::vector<float> &cy, 
           std::vector<float> &dy,
           cv::Point2f &np);
*/

void createLarvaContourPacked(cv::Point &first, 
                              size_t &size,
                              std::string &STR,
                              cvb::CvBlob &blob);

void createLarvaContourPoints(cv::Mat &lrvROI,
                              cvb::CvBlob &blob,
                              int type=CV_8UC1,
                              int PADDING=0);

void lengthAreaPerimeter(double a,double b,double &length, double &width);

double getGreyValue(cv::Mat &larvaROI, cvb::CvBlob &blob,cv::Mat &grey_frame);

double getPerimeter(cvb::CvBlob &blob);

double getSurroundingSize(cv::Point2f &point, cvb::CvBlob &blob,cv::Mat &grey_frame);
double getSurroundingSize(cv::Point2f &point, cvb::CvBlob &blob,cv::Mat &grey_frame,cv::Mat &preFrame);
double plotAngle(cvb::CvBlob *blob,cv::Mat &ROIimg,int PAD=0);
double angle( cv::Point2f &pt1, cv::Point2f &pt0, cv::Point2f &pt2 );
double angleD( cv::Point2f &pt1, cv::Point2f &pt0, cv::Point2f &pt2 );
double angleC( cv::Point2f &pt1, cv::Point2f &pt0, cv::Point2f &pt2 );
cv::Point2f perp(cv::Point2f &p1, cv::Point2f &p2, 
                 cv::Point2f &p3, double mag=1.0, 
                 bool right=true);
bool intersection(cv::Point2f o1, cv::Point2f p1, 
                  cv::Point2f o2, cv::Point2f p2, cv::Point2f &r);

void derivVec(std::vector<cv::Point2f> &in,
              std::vector<float> &p,
    std::vector<cv::Point2f> &d);

void deriv2Vec(
    std::vector<cv::Point2f> &d1,
              std::vector<float> &p,
    std::vector<cv::Point2f> &d2);

void curvVec(std::vector<cv::Point2f> &in,
             std::vector<float> &p,
             std::vector<float> &c);

void getBestCurvature(std::vector<float> &c,
                      std::vector<size_t> &di,
                      std::vector<float> &dmax
                      );
void getBestCurvatureS(std::vector<float> &c,
                      std::map<float,size_t> &curvatures,
                      std::vector<size_t> &di,
                      std::vector<float> &dmax,
                      double &variance
                      );
namespace std{
template<typename data>
  data getNeighboursAvg(std::vector<data> &in,int idx, int range, data initval)
  {
    data sum=initval;
    for(int i=(idx-range/2);i<=(idx+range/2);i++)
    {
      int id=i;
      if(i<0)
        id=in.size()+i-1;
      if(i>(int) in.size()-1)
        id=i-in.size();
      sum+=in[id];
    }
    data retval=sum*(1.0/range);
    return retval;
  }

template<typename data>
  void smoothVecMap(std::vector<data> &in,
                 std::map<data,size_t> &out,
                 std::vector<data> &out2,
                 int range, 
                 data initval)
  {
    for(size_t i=0;i<in.size();i++)
    {
      data ret;
      ret=getNeighboursAvg(in,i,range,initval);
      out[ret]=i;
      out2.push_back(ret);
    }
  }

template<typename data>
  void smoothVec(std::vector<data> &in,
                 std::vector<data> &out, 
                 int range, 
                 data initval)
  {
    for(size_t i=0;i<in.size();i++)
    {
      data ret;
      ret=getNeighboursAvg(in,i,range,initval);
      out.push_back(ret);
    }
  }
}

#endif
