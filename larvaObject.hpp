#ifndef __LRVTRACK_LARVAOBJECT_HPP__
#define __LRVTRACK_LARVAOBJECT_HPP__
#include "cvblob.h"
#include <vector>
#include "larvaSkel.hpp"

/*
 * Class containing information about each larva and its history on the plate.
 */
class larvaObject
{
public:
  unsigned int start_frame;
  unsigned int cur_frame;
  unsigned int lifetimeWithStats;
  unsigned int lastIdxWithStats;
  unsigned int larva_ID;
  unsigned int parentBlobID;
  bool isCluster;
  std::vector<cvb::CvBlob> blobs; //Blob for each frame for a given larva
  std::vector<float> area;
  float area_mean;
  float area_sum;
  float area_max;
  float area_min;

  std::vector<float> grey_value;
  float grey_value_mean;
  float grey_value_sum;
  float grey_value_max;
  float grey_value_min;

  std::vector<float> length;
  float length_mean;
  float length_sum;
  float length_max;
  float length_min;

  std::vector<float> perimeter;
  float perimeter_mean;
  float perimeter_sum;
  float perimeter_max;
  float perimeter_min;

  std::vector<float> width;
  float width_mean;
  float width_sum;
  float width_max;
  float width_min;

  std::vector<float> roundness;
  std::vector<float> angular_speed;

  std::vector<cv::Point> centroids;
  std::vector<float> centroid_speed_x;
  std::vector<float> centroid_speed_y;

  std::vector<float> centroid_distance_x;
  std::vector<float> centroid_distance_y;
  float centroid_distance_x_sum;
  float centroid_distance_y_sum;

  std::vector<unsigned int> inCluster;

  std::vector<larvaSkel> lrvskels;

  std::vector<cv::Point> heads;
  std::vector<cv::Point> tails;

  larvaObject():
    start_frame(0),
    cur_frame(0),
    lifetimeWithStats(0),
    lastIdxWithStats(0),
    larva_ID(0),
    parentBlobID(0),
    isCluster(false),
    area_mean(0),
    area_sum(0),
    area_max(0),
    area_min(0),
    grey_value_mean(0),
    grey_value_sum(0),
    grey_value_max(0),
    grey_value_min(0),
    length_mean(0),
    length_sum(0),
    length_max(0),
    length_min(0),
    perimeter_mean(0),
    perimeter_sum(0),
    perimeter_max(0),
    perimeter_min(0),
    width_mean(0),
    width_sum(0),
    width_max(0),
    width_min(0),
    centroid_distance_x_sum(0),
    centroid_distance_y_sum(0)
  {}

  int switchFaultyAssignment(
    std::map<unsigned int,std::vector<unsigned int> > &detected_clusters,
    std::map<unsigned int,larvaObject> &detected_larvae
  );
};
#endif
