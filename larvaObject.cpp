/*
 * Class containing information about each larva and its history on the plate.
 */
#include "larvaObject.hpp"
#include <iostream>

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

  larvaObject &otherLarva=detected_larvae[otherLarvaID];
  int exchangeIdx=0;

}
