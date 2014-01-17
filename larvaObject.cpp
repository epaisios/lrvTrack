/*
 * Class containing information about each larva and its history on the plate.
 */
#include "larvaObject.hpp"
#include <iostream>

using namespace std;

void larvaObject::dump() const
{
  cerr << endl;
  cerr << "============== Larva ID:" << larva_ID << "===============" <<  endl;
  cerr << "start_frame: " << start_frame << endl;
  cerr << "lifetimeWithStats: " << lifetimeWithStats << endl;
  cerr << "lastBlobWithStats: " << lastBlobWithStats << endl;
  cerr << "lastFrameWithStats: " << lastFrameWithStats << endl;
  cerr << "larva_ID: " << larva_ID << endl;
  cerr << "old_ID: " << old_ID << endl;
  cerr << "parentBlobID: " << parentBlobID << endl;
  cerr << "isCluster: " << isCluster << endl;
  cerr << "vector capture_times: " << endl;
  cerr << "  ";
  cerr << "size: " << capture_times.size() << endl;
  cerr << "  ";
  cerr << "content: " << printVector(capture_times) << endl;
  cerr << "vector blobs: " << endl;
  cerr << "  ";
  cerr << "size: " << capture_times.size() << endl;

  cerr << "vector area: " << endl;
  cerr << "  ";
  cerr << "size: " << area.size() << endl;
  cerr << "  ";
  cerr << "content: " << printVector(area) << endl;
  cerr << "area_mean: " << area_mean << endl;
  cerr << "area_sum: " << area_sum << endl;
  cerr << "area_max: " << area_max << endl;
  cerr << "area_min: " << area_min << endl;

  cerr << "vector grey_value: " << endl;
  cerr << "  ";
  cerr << "size: " << grey_value.size() << endl;
  cerr << "  ";
  cerr << "content: " << printVector(grey_value) << endl;
  cerr << "grey_value_mean: " << grey_value_mean << endl;
  cerr << "grey_value_sum: " << grey_value_sum << endl;
  cerr << "grey_value_max: " << grey_value_max << endl;
  cerr << "grey_value_min: " << grey_value_min << endl;

  cerr << "vector length: " << endl;
  cerr << "  ";
  cerr << "size: " << length.size() << endl;
  cerr << "  ";
  cerr << "content: " << printVector(length) << endl;
  cerr << "length_mean: " << length_mean << endl;
  cerr << "length_sum: " << length_sum << endl;
  cerr << "length_max: " << length_max << endl;
  cerr << "length_min: " << length_min << endl;

  cerr << "vector perimeter: " << endl;
  cerr << "  ";
  cerr << "size: " << perimeter.size() << endl;
  cerr << "  ";
  cerr << "content: " << printVector(perimeter) << endl;
  cerr << "perimeter_mean: " << perimeter_mean << endl;
  cerr << "perimeter_sum: " << perimeter_sum << endl;
  cerr << "perimeter_max: " << perimeter_max << endl;
  cerr << "perimeter_min: " << perimeter_min << endl;

  cerr << "vector width: " << endl;
  cerr << "  ";
  cerr << "size: " << width.size() << endl;
  cerr << "  ";
  cerr << "content: " << printVector(width) << endl;
  cerr << "width_mean: " << width_mean << endl;
  cerr << "width_sum: " << width_sum << endl;
  cerr << "width_max: " << width_max << endl;
  cerr << "width_min: " << width_min << endl;

  cerr << "vector headBodyAngle: " << endl;
  cerr << "  ";
  cerr << "size: " << headBodyAngle.size() << endl;
  cerr << "  ";
  cerr << "content: " << printVector(headBodyAngle) << endl;

  cerr << "vector orientationAngle: " << endl;
  cerr << "  ";
  cerr << "size: " << orientationAngle.size() << endl;
  cerr << "  ";
  cerr << "content: " << printVector(orientationAngle) << endl;

  cerr << "vector roundness: " << endl;
  cerr << "  ";
  cerr << "size: " << roundness.size() << endl;
  cerr << "  ";
  cerr << "content: " << printVector(roundness) << endl;

  cerr << "vector angular_speed: " << endl;
  cerr << "  ";
  cerr << "size: " << angular_speed.size() << endl;
  cerr << "  ";
  cerr << "content: " << printVector(angular_speed) << endl;

  cerr << "vector midpoint_speed_x: " << endl;
  cerr << "  ";
  cerr << "size: " << midpoint_speed_x.size() << endl;
  cerr << "  ";
  cerr << "content: " << printVector(midpoint_speed_x) << endl;

  cerr << "vector midpoint_speed_y: " << endl;
  cerr << "  ";
  cerr << "size: " << midpoint_speed_y.size() << endl;
  cerr << "  ";
  cerr << "content: " << printVector(midpoint_speed_y) << endl;

  cerr << "max_midpoint_speed: " << max_midpoint_speed << endl;
  cerr << "min_midpoint_speed: " << min_midpoint_speed << endl;

  cerr << "vector centroid_speed_x: " << endl;
  cerr << "  ";
  cerr << "size: " << centroid_speed_x.size() << endl;
  cerr << "  ";
  cerr << "content: " << printVector(centroid_speed_x) << endl;

  cerr << "vector centroid_speed_y: " << endl;
  cerr << "  ";
  cerr << "size: " << centroid_speed_y.size() << endl;
  cerr << "  ";
  cerr << "content: " << printVector(centroid_speed_y) << endl;

  cerr << "max_centroid_speed: " << max_centroid_speed << endl;
  cerr << "min_centroid_speed: " << min_centroid_speed << endl;

}

int larvaObject::switchFaultyAssignment(
  std::map<unsigned int,std::vector<unsigned int> > &detected_clusters,
  std::map<unsigned int,larvaObject> &detected_larvae)
{
  int i=inCluster.size()-1;
  int collisionIndex;
  int exchangeDuration; //in Frames up to current moment
  unsigned int blob;

// Look for last collision
  for (; i>=0; i--)
    {
      if (inCluster[i]>0)
        {
          blob=inCluster[i];
          collisionIndex=i;
          break;
        }
    }
  if(i==0)
    {
      std::cerr << "No Collision found, do nothing" << std::endl;
      return(0);
    }

// Find range of values that need to be switched

  // exchangeDuration:
  // The duration of the exchange. In this simple case,
  // it also provides us with the range.
  // For all values:
  //    Range: array.size()-1-exchangeDuration  up to array.size()-1
  // Here it's -2 since the index is always shifted by
  // 1 (starts at 0 and we do not want to change the
  // information while the blob was still there (which is
  // indicated by the collision index)
  exchangeDuration=inCluster.size()-collisionIndex-2;

// Find the other larva involved
  unsigned int otherLarvaID=0;
  std::vector<unsigned int>::iterator otherLarvaIT=detected_clusters[blob].begin();
  for (; otherLarvaIT!=detected_clusters[blob].begin(); ++otherLarvaIT)
    {
      if (*otherLarvaIT!=larva_ID)
        {
          otherLarvaID=*otherLarvaIT;
          break;
        }
    }

// Exchange the values

  //larvaObject &otherLarva=detected_larvae[otherLarvaID];
  //int exchangeIdx=0;

  return 0;
}
