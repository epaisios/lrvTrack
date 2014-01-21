#ifndef __LRVTRACK_LARVAOBJECT_HPP__
#define __LRVTRACK_LARVAOBJECT_HPP__
#include "cvblob.h"
#include <vector>
#include "larvaDistanceMap.hpp"

/*
 * Class containing information about each blob
 * The class contains a bunch of synchronized vectors.
 */

class blobObject_OL{
  private:
    // ID of the larva
    unsigned int ID;

    // ID of the originating blob of the larva.
    //  if the blob was there from the begining then
    //    originating_ID=ID
    //  if the blob is comming after dissapearing
    //    originating_ID=0
    unsigned int originating_ID;

    // ID where the blob entered before it dissappeared
    //  if the blob ends without a collision or dissapearing
    //  ending_ID=ID
    //  if the blob dissapears
    //  ending_ID=0
    
    unsigned int ending_ID;

    // Frame number of Nth data instance
    std::vector<unsigned int> FRAME;

    // Frames of existence until the Nth data instance
    std::vector<unsigned int> FRAMES_COUNT;

    //Vectors of length, width,size, perimeter, grey_value,various speeds etc.
    std::vector<double> length;
    std::vector<double> width;
    std::vector<double> size;
    std::vector<double> perimeter;
    std::vector<double> greyvalue;
    std::vector<double> spine0_speedx;
    std::vector<double> spine0_speedy;
    std::vector<double> spineM_speedx;
    std::vector<double> spineM_speedy;
    std::vector<double> spineL_speedx;
    std::vector<double> spineL_speedy;
    std::vector<double> centroid_speedx;
    std::vector<double> centroid_speedy;

    //Vector of centroids
    std::vector<double> centroids;

    //Vector of detailed Contour, spine reconstructions
    std::vector<larvaDistanceMap> larvaDistances;

    //Vector of original blobs
    std::vector<cvb::CvBlob> blobs;

    bool collisionObject;
};

/*
 * Class containing information about each larva and its history on the plate.
 */
class larvaObject
{
public:
  unsigned int start_frame;
  unsigned int lifetimeWithStats;
  unsigned int lastBlobWithStats;
  unsigned int lastFrameWithStats;
  unsigned int larva_ID;
  unsigned int old_ID;
  unsigned int parentBlobID;
  std::vector<unsigned int> childrenIDs;
  bool isCluster;
  std::vector<double> capture_times;
  std::vector<cvb::CvBlob> blobs; //Blob for each frame for a given larva
  std::vector<double> area;
  double area_mean;
  double area_sum;
  double area_max;
  double area_min;

  std::vector<double> grey_value;
  double grey_value_mean;
  double grey_value_sum;
  double grey_value_max;
  double grey_value_min;

  std::vector<double> length;
  double length_mean;
  double length_sum;
  double length_max;
  double length_min;

  std::vector<double> perimeter;
  double perimeter_mean;
  double perimeter_sum;
  double perimeter_max;
  double perimeter_min;

  std::vector<double> width;
  double width_mean;
  double width_sum;
  double width_max;
  double width_min;

  std::vector<double> headBodyAngle;
  std::vector<double> orientationAngle;

  std::vector<double> roundness;
  std::vector<double> angular_speed;

  std::vector<cv::Point2f> centroids;
  std::vector<cv::Point2f> centroids_full;
  std::vector<double> midpoint_speed_x;
  std::vector<double> midpoint_speed_y;
  double max_midpoint_speed;
  double min_midpoint_speed;

  std::vector<double> centroid_speed_x;
  std::vector<double> centroid_speed_y;
  double max_centroid_speed;
  double min_centroid_speed;

  std::vector<double> centroid_distance_x;
  std::vector<double> centroid_distance_y;
  double centroid_distance_x_sum;
  double centroid_distance_y_sum;

  std::vector<unsigned int> inCluster;

  std::vector<larvaDistanceMap> lrvDistances;

  std::vector<cv::Point2f> heads;
  std::vector<cv::Point2f> tails;

  larvaObject():
    start_frame(0),
    lifetimeWithStats(0),
    lastBlobWithStats(0),
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
  
  void dump() const;
  int switchFaultyAssignment(
    std::map<unsigned int,std::vector<unsigned int> > &detected_clusters,
    std::map<unsigned int,larvaObject> &detected_larvae
  );
};
#endif
