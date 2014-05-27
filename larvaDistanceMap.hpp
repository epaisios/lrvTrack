#ifndef __LRVTRACK_LARVADISTANCEMAP_HPP
#define __LRVTRACK_LARVADISTANCEMAP_HPP
#include <opencv2/core/core.hpp>
#include <vector>
#include "cvblob.h"

#define SPINE_SEGMENTS 12

extern size_t CURRENT_FRAME;
typedef std::pair<cv::Point2f,cv::Point2f > PointPair;

double p2fdist(cv::Point2f &a, cv::Point2f &b);
double p2fdist(double x1,double y1, double x2, double y2);

class larvaDistanceMap
{
private:
  std::vector<double> distances;
  cv::Mat px;
  cv::Mat py;
public:
  double MaxDist;
  double WidthDist;
  double WidthAvg;
  double WidthSD;
  size_t spineSegments;
  PointPair MaxDistPoints;
  cv::Point2f MidPoint;
  cv::Point2f p20;
  cv::Point2f p80;
  std::vector<cv::Point2f> points;
  std::vector<cv::Point2f> Spine;
  std::vector<double> Angles;
  std::vector<double> Widths;
  std::vector<PointPair> spinePairs;
  double maxAngle;
  double firstHalfWidthsSum;
  double secondHalfWidthsSum;
  double curvatureBias; // 0.0 <= Back <= 0.5 <= Front <= 1.0
  double curvatureVariance; // 0.0 <= Back <= 0.5 <= Front <= 1.0
  int maxAngleLocation;
  class my2ndPoint
  {
  private:
    larvaDistanceMap &parent;
    int p1;
  public:
    my2ndPoint(larvaDistanceMap& p, int fstPoint):parent(p),p1(fstPoint) {}

    double &operator[](int p2)
    {
      return parent.distances[p1*parent.points.size()+p2];
    }
    double &operator[](cv::Point2f p2)
    {
      int index_p2=std::distance(parent.points.begin(),
                                 std::find(parent.points.begin(),
                                           parent.points.end(),
                                           p2)
                                );
      return parent.distances[p1*parent.points.size()+index_p2];
    }

  };

  my2ndPoint operator[](int p1)
  {
    return my2ndPoint(*this,p1);
  }
  my2ndPoint operator[](cv::Point2f p1)
  {
    int index_p1=std::distance(points.begin(),
                               std::find(points.begin(),
                                         points.end(),
                                         p1)
                              );
    return my2ndPoint(*this,index_p1);
  }
  friend class my2ndPoint;

  larvaDistanceMap(std::vector<cv::Point2f> ps):spineSegments(SPINE_SEGMENTS),
                                                points(ps),
                                                maxAngle(DBL_MIN),
                                                maxAngleLocation(-1)
  {
    //distances.reserve((points.size())*(points.size()));
    Spine.resize(spineSegments);
  }

  void getPxPy(cv::Mat &x,cv::Mat &y)
  {
    x=px;
    y=py;
  }
  void getDistances(cv::Point2f p1)
  {
  }

};
void lBFS(int p1,
          std::vector<cv::Point2f> &Points ,
          larvaDistanceMap &Distances
         );

void computeInnerDistances(cvb::CvBlob &blob,
                           larvaDistanceMap &Distances,
                           cv::Point2f &MidPoint);

double distanceBetweenSpines(larvaDistanceMap &a,larvaDistanceMap &b);

void fixContour(
    cvb::CvBlob &blob,
    larvaDistanceMap &Distances,
    size_t RES,
    cv::Mat &frame,
    cv::Mat &previousFrame,
    std::vector<cv::Point2f> *heads=NULL,
    std::vector<cv::Point2f> *tails=NULL,
    std::vector<cvb::CvBlob> *blobs=NULL);

void computeSpine(
    cvb::CvBlob &blob,
    larvaDistanceMap &Distances,
    cv::Mat &frame
    );

#endif
