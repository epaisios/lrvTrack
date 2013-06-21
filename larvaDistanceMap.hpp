#ifndef __LRVTRACK_LARVADISTANCEMAP_HPP
#define __LRVTRACK_LARVADISTANCEMAP_HPP
#include <opencv2/core/core.hpp>
#include <vector>
#include "lrvTrackBase.hpp"
#include "cvblob.h"

//typedef std::unordered_map<PointPair, double> DistanceMap;

class larvaDistanceMap{
  private:
    std::vector<double> distances; 
    cv::Mat px;
    cv::Mat py;
  public:
    double MaxDist;
    double WidthDist;
    PointPair MaxDistPoints;
    PointPair WidthDistPoints;
    std::vector<cv::Point> &points;
    class my2ndPoint 
    {
      private:
        larvaDistanceMap &parent;
        int p1;
      public:
      my2ndPoint(larvaDistanceMap& p, int fstPoint):parent(p),p1(fstPoint){}

      double &operator[](int p2)
      {
        return parent.distances[p1*parent.points.size()+p2];
      }
      double &operator[](cv::Point p2)
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
    my2ndPoint operator[](cv::Point p1)
    {
      int index_p1=std::distance(points.begin(),
          std::find(points.begin(),
            points.end(),
            p1)
          );
      return my2ndPoint(*this,index_p1);
    }
    friend class my2ndPoint;

    larvaDistanceMap(std::vector<cv::Point> &ps):points(ps){
      distances.reserve((points.size())*(points.size()));
    }
    void getPxPy(cv::Mat &x,cv::Mat &y){
      x=px;
      y=py;
    }
    void getDistances(cv::Point p1)
    {
    }

};
void lBFS(int p1, 
          std::vector<cv::Point> &Points ,
          larvaDistanceMap &Distances
          );

void computeInnerDistances(cvb::CvBlob &blob,
                           larvaDistanceMap &Distances,
                           cv::Point &MidPoint);


#endif
