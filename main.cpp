#include <lrvTrack.hpp>
#include <fstream>
#ifndef _WIN32
#include <sys/time.h>
#endif
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/timer/timer.hpp>
#include <iomanip>
#include <string>

#include "cvblob.h"
#include "lrvTrackBase.hpp"
#include "blobUtils.hpp"
#include "larvaDistanceMap.hpp"

namespace po = boost::program_options;
namespace fs = boost::filesystem;

void findHeadTail(std::vector<cv::Point> &startPoints,
                  larvaObject &lrv,
                  cv::Point &Head,
                  cv::Point &Tail,
                  bool force_SurroundingValSearch=false)
{
  int max=0,min=65535;
  int i=0;
  if (startPoints.size()<2)
    {
      Head.x=0.0;
      Head.y=0.0;
      Tail.x=0.0;
      Tail.y=0.0;
    }

  if( lrv.start_frame==CURRENT_FRAME ||
      lrv.inCluster.back()>0 ||
      force_SurroundingValSearch ||
      lrv.roundness.back()<=3.0 ||
      (lrv.heads.back().x == 0 && lrv.heads.back().y==0 &&
       lrv.tails.back().x == 0 && lrv.tails.back().y==0)
    )
    {
      for (i=0; i<startPoints.size(); ++i)
        {
          cv::Point p=startPoints[i];
          cvb::CvBlob blob=lrv.blobs.back();
          double tmpsz=getSurroundingSize(p,blob,grey_frame);
          if (tmpsz>max)
            {
              Tail=startPoints[i];
              max=tmpsz;
            }
          if (tmpsz<min)
            {
              Head=startPoints[i];
              min=tmpsz;
            }
        }
      if (max==min)
        std::cerr << "PROBLEM AREAS AROUND HEAD AND TAIL ARE THE SAME" << std::endl;
      if (Head==Tail && i!=2)
        {
          std::cerr << "PROBLEM HEAD AND TAIL ARE THE SAME" << std::endl;
        }
    }
  else
    {
      double hmin=65535;
      double tmin=65535;
      for (i=0; i<startPoints.size(); ++i)
        {
          double hdiff=diff(lrv.heads.back(),startPoints[i]);
          double tdiff=diff(lrv.tails.back(),startPoints[i]);
          if (hdiff<hmin)
            {
              Head=startPoints[i];
              hmin=hdiff;
            }
          if (tdiff<tmin)
            {
              Tail=startPoints[i];
              tmin=tdiff;
            }

        }
    }
}

void updateLarvae(cvb::CvBlobs &In, cvb::CvBlobs &Prev)
{
  cvb::CvBlobs::iterator it=In.begin();
  std::vector<unsigned int> larvaeInClusters;

  while (it!=In.end())
    {
      unsigned int ID=(*it).first;
      cvb::CvBlob blob=*((*it).second);
      cv::Mat larvaROI,cntPoints;
      createLarvaContour(larvaROI,blob);
      createLarvaContourPoints(cntPoints,blob);

      //cv::imshow("Current Larva",larvaROI);
      //cv::imshow("Contour Points",cntPoints);
      //cv::waitKey(1);
      std::map<unsigned int,larvaObject>::iterator curLarva;
      // NEW LARVA OBJECT!
      if ((curLarva=detected_larvae.find(ID))==detected_larvae.end())
        {
          // Create and allocate the new object
          larvaObject newLarva;
          //
          // Initialize the speed to 0
          newLarva.centroid_speed_x.push_back(0);
          newLarva.centroid_speed_y.push_back(0);

          // State that the larva is not in a blob
          newLarva.inCluster.push_back(false);

          cv::Point centroid=cv::Point(
                               static_cast<int>(blob.centroid.x-blob.minx+ROI_PADDING+0.5),
                               static_cast<int>(blob.centroid.y-blob.miny+ROI_PADDING+0.5)
                             );
          newLarva.centroids.push_back(centroid);
          ++newLarva.lifetimeWithStats;
          newLarva.lastIdxWithStats=CURRENT_FRAME;
          // Set the frame of it's existence
          newLarva.start_frame=CURRENT_FRAME;

          // Add the blob of the larva to its blob history
          newLarva.blobs.push_back(blob);
          // Give the larva the necessary ID
          newLarva.larva_ID=ID;

          if (detected_clusters.find(ID)==detected_clusters.end())
            {
              newLarva.isCluster=false;

              //Initialize the area values
              newLarva.area.push_back(blob.area);
              newLarva.area_mean=blob.area;
              newLarva.area_sum=blob.area;
              newLarva.area_max=newLarva.area_min=blob.area;
              newLarva.area_min=newLarva.area_min=blob.area;

              double greyVal=getGreyValue(larvaROI,blob,grey_frame);
              newLarva.grey_value.push_back(greyVal);
              newLarva.grey_value_mean = greyVal;
              newLarva.grey_value_sum= greyVal;
              newLarva.grey_value_max = greyVal;
              newLarva.grey_value_min = greyVal;

              double perimeter=getPerimeter(blob);
              newLarva.perimeter.push_back(perimeter);
              newLarva.perimeter_mean=perimeter;
              newLarva.perimeter_sum=perimeter;
              newLarva.perimeter_max=perimeter;
              newLarva.perimeter_min=perimeter;

              // Initialize the speed to 0
              newLarva.centroid_speed_x.push_back(0);
              newLarva.centroid_speed_y.push_back(0);

              newLarva.roundness.push_back((perimeter*perimeter)/(2*CV_PI*blob.area));

              // State that the larva is not in a blob
              newLarva.inCluster.push_back(false);

              cv::Point centroid=cv::Point(
                                   static_cast<int>(blob.centroid.x-blob.minx+ROI_PADDING+0.5),
                                   static_cast<int>(blob.centroid.y-blob.miny+ROI_PADDING+0.5)
                                 );
              newLarva.centroids.push_back(centroid);
              larvaSkel newLarvaSkel(larvaROI,centroid);
              newLarva.lrvskels.push_back(newLarvaSkel);

              // In this block we compute the inner distances for the larva
              std::vector<cv::Point> cntPoints;
              blobToPointVector(blob,cntPoints);
              larvaDistanceMap Distances(cntPoints);
              computeInnerDistances(blob,Distances,newLarvaSkel.MidPoint);
              newLarva.length.push_back(Distances.MaxDist);
              newLarva.length_mean = Distances.MaxDist;
              newLarva.length_sum= Distances.MaxDist;
              newLarva.length_max = Distances.MaxDist;
              newLarva.length_min = Distances.MaxDist;
              PointPair MAXPair=Distances.MaxDistPoints;

              cv::Point Head,Tail;
              std::vector<cv::Point> startPoints;
              startPoints.push_back(cv::Point(
                                      MAXPair.first.x-blob.minx+ROI_PADDING,
                                      MAXPair.first.y-blob.miny+ROI_PADDING)
                                   );
              startPoints.push_back(cv::Point(
                                      MAXPair.second.x-blob.minx+ROI_PADDING,
                                      MAXPair.second.y-blob.miny+ROI_PADDING)
                                   );
              newLarva.angular_speed.push_back(0);

              findHeadTail(startPoints,newLarva,Head,Tail);
              newLarva.heads.push_back(Head);
              newLarva.tails.push_back(Tail);

              newLarva.width.push_back(Distances.WidthDist);
              newLarva.width_mean = Distances.WidthDist;
              newLarva.width_sum= Distances.WidthDist;
              newLarva.width_max = Distances.WidthDist;
              newLarva.width_min = Distances.WidthDist;


              if(DEBUG_INFO!=0)
                {
                  std::cout << CURRENT_FRAME <<
                            " , " << newLarva.larva_ID <<
                            " , " << blob.area <<
                            " , " << Distances.MaxDist <<
                            " , " << greyVal <<
                            " , " << perimeter <<
                            " , " << Distances.WidthDist <<
                            " , " << newLarva.roundness.back() <<
                            " , " << newLarva.inCluster.back() <<
                            std::endl;
                }
            }
          else
            {
              // Larva is a blob
              // IMPORTANT: When the larva IS a blob NOTHING is updated
              //            The only fields that get updated are:
              //             - isCluster
              //             - centroid
              //             - centroid speeds
              //
              // IMPORTANT: When the larva is INSIDE a blob NOTHING is updated
              //            The only field that gets updated is:
              //             - inCluster
              //             - blob
              //
              newLarva.isCluster=true;
              std::vector<unsigned int>::iterator cl;
              for ( cl=detected_clusters[ID].begin() ; cl!=detected_clusters[ID].end() ; ++cl)
                {
                  larvaObject &clusterLarva=detected_larvae[*cl];
                  clusterLarva.inCluster.push_back(ID);
                  clusterLarva.blobs.push_back(blob);
                  if(DEBUG_INFO!=0)
                    {
                      std::cout << CURRENT_FRAME <<
                                " , " << clusterLarva.larva_ID <<
                                " , " << blob.area <<
                                " , " << clusterLarva.length.back() <<
                                " , " << clusterLarva.grey_value.back() <<
                                " , " << clusterLarva.perimeter.back() <<
                                " , " << clusterLarva.width.back() <<
                                " , " << clusterLarva.roundness.back() <<
                                " , " << clusterLarva.inCluster.back() <<
                                std::endl;
                    }
                }

            }
          detected_larvae[ID]=newLarva;
        }
      // UPDATED LARVA OBJECT!
      else
        {
          //Reference to current larva
          larvaObject &cur_larva=(*curLarva).second;
          //Pointer for the previous blob
          cvb::CvBlob preBlob = cur_larva.blobs.back();

          // Set the ID of the larvaObject to the ID found TODO:Probably unnecessary
          cur_larva.larva_ID=ID;

          // Add the current blob to the blobs history of the larva
          cur_larva.blobs.push_back(blob);

          // Create the skeleton of the larva and add it to the skeletons history
          cv::Point centroid=cv::Point(
                               static_cast<int>(blob.centroid.x-blob.minx+ROI_PADDING+0.5),
                               static_cast<int>(blob.centroid.y-blob.miny+ROI_PADDING+0.5)
                             );

          cur_larva.centroids.push_back(centroid);

          // Update centroid_speeds (in pixel per second per axis)
          double FrameEllapsedSeconds=FrameEllapsedTime.wall/1000000000.0;
          cur_larva.centroid_speed_x.push_back(
            (blob.centroid.x - preBlob.centroid.x)/FrameEllapsedSeconds);
          cur_larva.centroid_speed_y.push_back(
            (blob.centroid.y - preBlob.centroid.y)/FrameEllapsedSeconds);

          double curAngle=cvb::cvAngle(&blob);
          double preAngle=cvb::cvAngle(&preBlob);

          cur_larva.angular_speed.push_back(cv::fast_abs(curAngle-preAngle)/FrameEllapsedSeconds);

          // Look if larva is a blob
          if ( (!cur_larva.isCluster) &&
               (detected_clusters.find(ID)==detected_clusters.end()))
            {

              ++cur_larva.lifetimeWithStats;
              cur_larva.lastIdxWithStats=cur_larva.area.size();

              larvaSkel newLarvaSkel(larvaROI,centroid);
              cur_larva.lrvskels.push_back(newLarvaSkel);

              // If not then:
              //  Update area values for larva.
              cur_larva.area.push_back(blob.area);

              cur_larva.area_mean=(cur_larva.area_mean+blob.area)/2;
              cur_larva.area_sum = cur_larva.area_sum + blob.area;
              if (cur_larva.area_max < blob.area)
                {
                  cur_larva.area_max=blob.area;
                }
              if (cur_larva.area_min > blob.area)
                {
                  cur_larva.area_min=blob.area;
                }


              // Try to find Head and Tail

              // Point coordinates for head/tail
              cv::Point Head,Tail;

              // Map to keep the distances of each point to all the others
              // Pair of points to keep the points with the Maximum distance (i.e. head and tail :) )
              std::vector<cv::Point> cntPoints;
              blobToPointVector(blob,cntPoints);
              larvaDistanceMap Distances(cntPoints);
              // Compute all the inner distances for the larva
              computeInnerDistances(blob,Distances,newLarvaSkel.MidPoint);
              cur_larva.length.push_back(Distances.MaxDist);
              cur_larva.length_mean=(cur_larva.length_mean+Distances.MaxDist)/2;
              cur_larva.length_sum=cur_larva.length_sum+Distances.MaxDist;
              if (cur_larva.area_max < Distances.MaxDist)
                {
                  cur_larva.area_max=Distances.MaxDist;
                }
              if (cur_larva.area_min > Distances.MaxDist)
                {
                  cur_larva.area_min=Distances.MaxDist;
                }
              PointPair MAXPair=Distances.MaxDistPoints;

              cur_larva.width.push_back(Distances.WidthDist);
              cur_larva.width_mean=(cur_larva.width_mean+Distances.WidthDist)/2;
              cur_larva.width_sum=cur_larva.width_sum+Distances.WidthDist;
              if (cur_larva.area_max < Distances.WidthDist)
                {
                  cur_larva.area_max=Distances.WidthDist;
                }
              if (cur_larva.area_min > Distances.WidthDist)
                {
                  cur_larva.area_min=Distances.WidthDist;
                }

              double greyVal=getGreyValue(larvaROI,blob,grey_frame);
              cur_larva.grey_value.push_back(greyVal);
              cur_larva.grey_value_mean=(cur_larva.grey_value_mean+greyVal)/2;
              cur_larva.grey_value_sum=cur_larva.grey_value_sum+greyVal;
              if (cur_larva.area_max < greyVal)
                {
                  cur_larva.grey_value_max=greyVal;
                }
              if (cur_larva.area_min > greyVal)
                {
                  cur_larva.grey_value_min=greyVal;
                }

              double perimeter=getPerimeter(blob);
              cur_larva.perimeter.push_back(perimeter);
              cur_larva.perimeter_mean=(cur_larva.perimeter_mean+perimeter)/2;
              cur_larva.perimeter_sum=cur_larva.perimeter_sum+perimeter;
              if (cur_larva.area_max < perimeter)
                {
                  cur_larva.perimeter_max=perimeter;
                }
              if (cur_larva.area_min > perimeter)
                {
                  cur_larva.perimeter_min=perimeter;
                }

              cur_larva.roundness.push_back((perimeter*perimeter)/(2*CV_PI*blob.area));

              // Construct a vector of points including both
              // head and tail to decide which is which
              std::vector<cv::Point> startPoints;

              // Points must be corrected to match the ROI including the padding set
              startPoints.push_back(cv::Point(
                                      MAXPair.first.x-blob.minx+ROI_PADDING,
                                      MAXPair.first.y-blob.miny+ROI_PADDING)
                                   );
              startPoints.push_back(cv::Point(
                                      MAXPair.second.x-blob.minx+ROI_PADDING,
                                      MAXPair.second.y-blob.miny+ROI_PADDING)
                                   );
              //execute the function to find which is which and assign appropriately
              //to Head/Tail.
              //findHeadTail(startPoints,cur_larva,Head,Tail,true);
              findHeadTail(startPoints,cur_larva,Head,Tail);

              //Add head and tail to the history
              cur_larva.heads.push_back(Head);
              cur_larva.tails.push_back(Tail);

              //state that larva is not detected as part of a blob of larvae
              //NOTE: It is important to perform this value setting !AFTER!
              //      searching for the head tail because, the head tail
              //      search behaves differently if the larva was previously
              //      part of a blob.
              cur_larva.inCluster.push_back(0);

              if(DEBUG_INFO!=0)
                {
                  std::cout << CURRENT_FRAME <<
                            " , " << cur_larva.larva_ID <<
                            " , " << blob.area <<
                            " , " << cur_larva.length.back() <<
                            " , " << cur_larva.grey_value.back() <<
                            " , " << cur_larva.perimeter.back() <<
                            " , " << cur_larva.width.back() <<
                            " , " << cur_larva.roundness.back() <<
                            " , " << cur_larva.inCluster.back() <<
                            std::endl;
                }

            }
          else
            {
              // Larva is a blob
              // IMPORTANT: When the larva IS a blob NOTHING is updated
              //            The only fields that get updated are:
              //             - isCluster
              //             - blob
              //             - centroid
              //             - centroid speeds
              //
              // IMPORTANT: When the larva is INSIDE a blob NOTHING is updated
              //            The only field that gets updated is:
              //             - inCluster
              //             - blob
              //
              cur_larva.isCluster=true;
              std::vector<unsigned int>::iterator cl;
              for ( cl=detected_clusters[ID].begin() ; cl!=detected_clusters[ID].end() ; ++cl)
                {
                  larvaObject &clusterLarva=detected_larvae[*cl];
                  clusterLarva.inCluster.push_back(ID);
                  clusterLarva.blobs.push_back(blob);
                  if(DEBUG_INFO!=0)
                    {
                      std::cout << CURRENT_FRAME <<
                                " , " << clusterLarva.larva_ID <<
                                " , " << blob.area <<
                                " , " << clusterLarva.length.back() <<
                                " , " << clusterLarva.grey_value.back() <<
                                " , " << clusterLarva.perimeter.back() <<
                                " , " << clusterLarva.width.back() <<
                                " , " << clusterLarva.roundness.back() <<
                                " , " << clusterLarva.inCluster.back() <<
                                std::endl;
                    }
                }
            }
        }
      ++it;
    }
}

inline double SQUARE(double n)
{
  return n*n;
}

//Currently only for the 4 sized vector
inline void vmin(std::vector<double>& vals, std::vector<int>& seq)
{
  for (int i=0 ; i<vals.size() ; ++i )
    {
      int sum=0;
      for (int j=0; j<vals.size() ; ++j)
        {
          if (vals[i]>vals[j])
            {
              ++sum;
            }
        }
      seq[sum]=i;
    }
}

namespace std
{
template <typename number >
std::string printVector(std::vector<number> vec)
{
  std::stringstream sstm;
  //bool const is_number= std::is_arithmetic<number>::value;
  //static_assert( is_number, "Provided type is not an arithmetic type");
  sstm << "[" ;
  typename std::vector<number>::const_iterator i=vec.begin();
  sstm << *i ;
  ++i;
  for( ; i != vec.end(); ++i)
    sstm << ","<< *i;
  sstm << "]";
  return sstm.str();
}
}

void diverge_match_new(
  std::vector<unsigned int> &candidateLarvae,
  std::vector<unsigned int> &newLarvae,
  std::map<unsigned int, unsigned int> &newAssignments,
  cvb::CvBlobs &NEW)
{
  std::map<unsigned int, cv::Mat> candidateCovarMat;
  std::map<unsigned int, cv::Mat> candidateMeanMat;

  std::map<unsigned int, cv::Mat> newMeanMat;

  std::vector<unsigned int>::iterator cIT=candidateLarvae.begin();
  while(cIT!=candidateLarvae.end())
  {
    double size_avg=detected_larvae[*cIT].area_sum/
      detected_larvae[*cIT].area.size();
    double grey_value_avg=detected_larvae[*cIT].grey_value_sum/
      detected_larvae[*cIT].grey_value.size();
    double length_avg=detected_larvae[*cIT].length_sum/
      detected_larvae[*cIT].length.size();
    double perimeter_avg=detected_larvae[*cIT].perimeter_sum/
      detected_larvae[*cIT].perimeter.size();
    double width_avg=detected_larvae[*cIT].width_sum/
      detected_larvae[*cIT].width.size();

    cv::Mat InputArray;
    cv::hconcat(cv::Mat(detected_larvae[*cIT].area),
        cv::Mat(detected_larvae[*cIT].grey_value),
        InputArray);

    cv::hconcat(InputArray,
        cv::Mat(detected_larvae[*cIT].length),
        InputArray);

    cv::hconcat(InputArray,
        cv::Mat(detected_larvae[*cIT].perimeter),
        InputArray);

    cv::hconcat(InputArray,
        cv::Mat(detected_larvae[*cIT].width),
        InputArray);

    std::vector<double> meanVec;
    meanVec.push_back(size_avg);
    meanVec.push_back(grey_value_avg);
    meanVec.push_back(length_avg);
    meanVec.push_back(perimeter_avg);
    meanVec.push_back(width_avg);

    cv::Mat meanMat(meanVec);
    cv::Mat meanTMat;
    cv::transpose(meanMat,meanTMat);

    cv::Mat covarMat;

    cv::calcCovarMatrix(InputArray, covarMat,meanTMat,CV_COVAR_ROWS|CV_COVAR_NORMAL|CV_COVAR_USE_AVG);
    cv::invert(covarMat,covarMat,cv::DECOMP_SVD);

    candidateCovarMat[*cIT]=covarMat;
    candidateMeanMat[*cIT]=meanMat;

    ++cIT;
  }
  
  cIT=newLarvae.begin();
  while(cIT!=newLarvae.end())
  {
    //Setup of each new larva
    cv::Mat larvaROI;
    std::vector<cv::Point> newLarvaPoints;
    blobToPointVector(*NEW[*cIT],newLarvaPoints);
    createLarvaContour(larvaROI,(*NEW[*cIT]));
    larvaDistanceMap dstLarva(newLarvaPoints);
    cv::Point centroid;
    centroid.x=NEW[*cIT]->centroid.x;
    centroid.y=NEW[*cIT]->centroid.y;
    larvaSkel newLarvaSkel(larvaROI,centroid);
    computeInnerDistances(*NEW[*cIT],dstLarva,newLarvaSkel.MidPoint);

    double newSize=NEW[*cIT]->area;
    double newGreyval=getGreyValue(larvaROI,*NEW[*cIT],grey_frame);
    double newLength=dstLarva.MaxDist;
    double newPerimeter=getPerimeter(*NEW[*cIT]);
    double newWidth=dstLarva.WidthDist;

    std::vector<double> meanVec;
    meanVec.push_back(newSize);
    meanVec.push_back(newGreyval);
    meanVec.push_back(newLength);
    meanVec.push_back(newPerimeter);
    meanVec.push_back(newWidth);

    cv::Mat meanMat(meanVec);
    newMeanMat[*cIT]=meanMat;

  }

    // For each new Larva we look at the MH distances with the candidates
    //   if the min distance is less than MH_THRESHOLD assign it
    //   if the assignment is already done assign the smaller of the two
    //      and free the other.

  //std::map <unsigned int, unsigned int> newAssignments;
  std::map <unsigned int, unsigned int> oldAssignments;
  std::map <unsigned int, double> oldAssignmentDistances;

  cIT=newLarvae.begin();
  while(cIT!=newLarvae.end())
  {
    newAssignments[*cIT]=0;
    ++cIT;
  }

  bool changed=true;
  while (changed)
  {
    changed=false;
    cIT=newLarvae.begin();
    while(cIT!=newLarvae.end())
    {
      double minDist=9999.0;
      unsigned int minCand=0;
      if(newAssignments[*cIT]==0)
      {
        std::vector<unsigned int>::iterator oIT=candidateLarvae.begin();
        while(oIT!=candidateLarvae.end())
        {
          double Dist=cv::Mahalanobis(newMeanMat[*cIT],
              candidateMeanMat[*oIT],
              candidateCovarMat[*oIT]);
          if(Dist<LARVA_MAHALANOBIS_THRESHOLD && 
              Dist < minDist)
          {
            if(oldAssignmentDistances.find(*oIT)==oldAssignmentDistances.end())
            {
              minDist=Dist;
              minCand=*oIT;
            }
            else if(oldAssignmentDistances[*oIT]>Dist)
            {
              minDist=Dist;
              minCand=*oIT;
            }
          }
          ++oIT;
        }
        if(oldAssignmentDistances.find(minDist)==oldAssignmentDistances.end())
        {
          newAssignments[*cIT]=minCand;
          oldAssignments[minCand]=*cIT;
          oldAssignmentDistances[minCand]=minDist;
        }
        else if(oldAssignmentDistances[*oIT]>minDist)
        {
          newAssignments[*cIT]=minCand;
          newAssignments[oldAssignments[minCand]]=0;
          oldAssignments[minCand]=*cIT;
          changed=true;
        }
      }
      ++cIT;
    }
  }
}

void diverge_match(
  unsigned int &candidate_larva_a,
  unsigned int &candidate_larva_b,
  cvb::CvBlob  *newLarva1,
  cvb::CvBlob *newLarva2)
{
  larvaObject &LarvaA=detected_larvae[candidate_larva_a];
  larvaObject &LarvaB=detected_larvae[candidate_larva_b];
  cv::Point mdpLarva1, mdpLarva2;
  mdpLarva1.x=newLarva1->centroid.x;
  mdpLarva1.y=newLarva1->centroid.y;
  mdpLarva2.x=newLarva2->centroid.x;
  mdpLarva2.y=newLarva2->centroid.y;

  double size_a=LarvaA.area_sum/LarvaA.area.size();
  double size_b=LarvaB.area_sum/LarvaB.area.size();;
  double size_1=newLarva1->area;
  double size_2=newLarva2->area;

  cv::Mat larvaROI1,larvaROI2;
  createLarvaContour(larvaROI1,*newLarva1);
  createLarvaContour(larvaROI2,*newLarva2);
  double grey_value_a=LarvaA.grey_value_sum/LarvaA.grey_value.size();
  double grey_value_b=LarvaB.grey_value_sum/LarvaB.grey_value.size();
  double grey_value_1=getGreyValue(larvaROI1,*newLarva1,grey_frame);
  double grey_value_2=getGreyValue(larvaROI2,*newLarva2,grey_frame);

  // Length
  std::vector<cv::Point> newLarva1Points;
  std::vector<cv::Point> newLarva2Points;
  blobToPointVector(*newLarva1,newLarva1Points);
  blobToPointVector(*newLarva2,newLarva2Points);
  larvaDistanceMap dstLarva1(newLarva1Points), dstLarva2(newLarva2Points);

  cv::Point centroid1;
  centroid1.x=newLarva1->centroid.x;
  centroid1.y=newLarva1->centroid.y;
  cv::Point centroid2;
  centroid2.x=newLarva2->centroid.x;
  centroid2.y=newLarva2->centroid.y;

  larvaSkel newLarvaSkel1(larvaROI1,centroid1);
  larvaSkel newLarvaSkel2(larvaROI2,centroid2);

  computeInnerDistances(*newLarva1,dstLarva1,newLarvaSkel1.MidPoint);
  computeInnerDistances(*newLarva2,dstLarva2,newLarvaSkel2.MidPoint);


  double length_a=LarvaA.length_sum/LarvaA.length.size();
  double length_b=LarvaB.length_sum/LarvaB.length.size();
  double length_1=dstLarva1.MaxDist;
  double length_2=dstLarva2.MaxDist;

  double perimeter_a=LarvaA.perimeter_sum/LarvaA.perimeter.size();
  double perimeter_b=LarvaB.perimeter_sum/LarvaB.perimeter.size();
  double perimeter_1=getPerimeter(*newLarva1);
  double perimeter_2=getPerimeter(*newLarva2);

  double width_a=LarvaA.width_sum/LarvaA.width.size();
  double width_b=LarvaB.width_sum/LarvaB.width.size();
  double width_1=dstLarva1.WidthDist;
  double width_2=dstLarva2.WidthDist;

  cv::Mat InputArrayA;
  cv::Mat InputArrayB;
  cv::hconcat(cv::Mat(LarvaA.area),cv::Mat(LarvaA.grey_value),InputArrayA);
  cv::hconcat(InputArrayA,cv::Mat(LarvaA.length),InputArrayA);
  cv::hconcat(InputArrayA,cv::Mat(LarvaA.perimeter),InputArrayA);
  cv::hconcat(InputArrayA,cv::Mat(LarvaA.width),InputArrayA);

  cv::hconcat(cv::Mat(LarvaB.area),cv::Mat(LarvaB.grey_value),InputArrayB);
  cv::hconcat(InputArrayB,cv::Mat(LarvaB.length),InputArrayB);
  cv::hconcat(InputArrayB,cv::Mat(LarvaB.perimeter),InputArrayB);
  cv::hconcat(InputArrayB,cv::Mat(LarvaB.width),InputArrayB);

  /*
     std::cerr << InputArrayA << std::endl;
     std::cerr << InputArrayB << std::endl;
     std::cerr << "===========================================" << std::endl;
     */
  std::vector<double> meanVecA;
  meanVecA.push_back(size_a);
  meanVecA.push_back(grey_value_a);
  meanVecA.push_back(length_a);
  meanVecA.push_back(perimeter_a);
  meanVecA.push_back(width_a);

  std::vector<double> meanVecB;
  meanVecB.push_back(size_b);
  meanVecB.push_back(grey_value_b);
  meanVecB.push_back(length_b);
  meanVecB.push_back(perimeter_b);
  meanVecB.push_back(width_b);

  cv::Mat meanMatA(meanVecA);
  cv::Mat meanMatB(meanVecB);
  cv::Mat meanTMatA, meanTMatB;
  cv::transpose(meanMatA,meanTMatA);
  cv::transpose(meanMatB,meanTMatB);
  cv::Mat covarMatA;
  cv::Mat covarMatB;
  /*
     printVector(LarvaA.area);
     printVector(LarvaA.grey_value);
     printVector(LarvaA.perimeter);

     printVector(LarvaB.area);
     printVector(LarvaB.grey_value);
     printVector(LarvaB.perimeter);
     */
  cv::calcCovarMatrix(InputArrayA, covarMatA,meanTMatA,CV_COVAR_ROWS|CV_COVAR_NORMAL|CV_COVAR_USE_AVG);
  cv::calcCovarMatrix(InputArrayB, covarMatB,meanTMatB,CV_COVAR_ROWS|CV_COVAR_NORMAL|CV_COVAR_USE_AVG);

  /*
     std::cerr << "CovarMatA" << covarMatA << std::endl;
     std::cerr << "CovarMatB" << covarMatB << std::endl;
     */
  cv::invert(covarMatA,covarMatA,cv::DECOMP_SVD);
  cv::invert(covarMatB,covarMatB,cv::DECOMP_SVD);

  std::vector<double> vec1;
  vec1.push_back(size_1);
  vec1.push_back(grey_value_1);
  vec1.push_back(length_1);
  vec1.push_back(perimeter_1);
  vec1.push_back(width_1);

  std::vector<double> vec2;
  vec2.push_back(size_2);
  vec2.push_back(grey_value_2);
  vec2.push_back(length_2);
  vec2.push_back(perimeter_2);
  vec2.push_back(width_2);

  cv::Mat mat1(vec1);
  cv::Mat mat2(vec2);

  double DistA1 = cv::Mahalanobis(mat1, meanMatA, covarMatA);
  double DistA2 = cv::Mahalanobis(mat2, meanMatA, covarMatA);
  double DistB1 = cv::Mahalanobis(mat1, meanMatB, covarMatB);
  double DistB2 = cv::Mahalanobis(mat2, meanMatB, covarMatB);

  double MHDiff = (DistA1+DistB2) - (DistA2+DistB1);
  std::vector<double> vals;
  vals.push_back(DistA1);
  vals.push_back(DistB2);
  vals.push_back(DistA2);
  vals.push_back(DistB1);

  std::vector<int> minv;
  minv.reserve(4);
  vmin(vals,minv);
  int mini=minv[0];
  if (DistA1+DistB2 > DistA2+DistB1)
    /*
       double miniDiff = cv::fast_abs(vals[minv[0]] - vals[minv[1]]);
       if (mini >=2 && MHDiff < 0 )
       {
       if (cv::fast_abs(MHDiff) > miniDiff)
       {
       mini=0;
       }
       }
       else if (mini < 2 && MHDiff > 0 )
       {
       if (cv::fast_abs(MHDiff) > miniDiff)
       {
       mini=2;
       }
       }

       if (mini >=2 )
       */
    {
      candidate_larva_a=newLarva2->label;
      candidate_larva_b=newLarva1->label;
      std::cerr << "MH<" << CURRENT_FRAME << "> : Assigning " << LarvaA.larva_ID << " -> 2 [" << newLarva2->centroid.x << ", " << newLarva2->centroid.y << "] MHDiff: " << MHDiff << " min i: " << mini <<  std::endl;
      std::cerr << "MH<" << CURRENT_FRAME << "> : Assigning " << LarvaB.larva_ID << " -> 1 [" << newLarva1->centroid.x << ", " << newLarva1->centroid.y << "] Dists: A1: " << DistA1 << " B2: " << DistB2 << ", A2: " << DistA2 << " B1: " << DistB1 << std::endl;
    }
  else
    {
      candidate_larva_a=newLarva1->label;
      candidate_larva_b=newLarva2->label;
      std::cerr << "MH<" << CURRENT_FRAME << "> : Assigning " << LarvaA.larva_ID << " -> 1 [" << newLarva1->centroid.x << ", " << newLarva1->centroid.y << "] MHDiff: " << MHDiff << " min i: " << mini << std::endl;
      std::cerr << "MH<" << CURRENT_FRAME << "> : Assigning " << LarvaB.larva_ID << " -> 2 [" << newLarva2->centroid.x << ", " << newLarva2->centroid.y << "] Dists: A1: " << DistA1 << " B2: " << DistB2 << ", A2: " << DistA2 << " B1: " << DistB1 << std::endl;
    }
}

//Returns all the larvae in the set Blobs within an area around "Blob".
// The area is determined by the size of the Blob (the longest side 
// of the bounding box multiplied by PADRatio).
//
// The vector nearbyLarvae is filled by those found sorted from closest
// to furthest.
void getNearbyLarvae(cvb::CvBlobs &Blobs, cvb::CvBlob *Blob, 
		            std::vector<unsigned int> &nearbyLarvae,double PADRatio=2)
{
  std::vector<double> distances;
	double MaxDist = std::max(Blob->maxx-Blob->minx,Blob->maxy-Blob->miny);
	MaxDist=PADRatio*MaxDist;

  cvb::CvBlobs::iterator It=Blobs.begin();
  while (It!=Blobs.end())
	{
		cvb::CvBlob *cBlob=It->second;
    if (cBlob->centroid.x < (Blob->maxx + MaxDist/2) &&
        cBlob->centroid.x > (Blob->minx - MaxDist/2) &&
        cBlob->centroid.y < (Blob->maxy + MaxDist/2) &&
        cBlob->centroid.y > (Blob->miny - MaxDist/2))
    {
      double DIST=cv::fast_abs(Blob->centroid.x - cBlob->centroid.x) +
        cv::fast_abs(Blob->centroid.y - cBlob->centroid.y);
      std::vector<double>::iterator dIt=distances.begin();
      std::vector<unsigned int>::iterator nIt=nearbyLarvae.begin();
      if(nearbyLarvae.size()>0)
      {
        while(dIt!=distances.end())
        {
          if(DIST < *dIt)
          {
            distances.insert(dIt,DIST);
            nearbyLarvae.insert(nIt,It->first);
            break;
          }
          ++dIt;
          ++nIt;
        }
        if(dIt==distances.end())
        {
          nearbyLarvae.push_back(It->first);
          distances.push_back(DIST);
        }
      }
      else
      {
        distances.push_back(DIST);
        nearbyLarvae.push_back(It->first);
      }
    }
		++It;
	}
}

void findCombination(cvb::CvBlobs &In, cvb::CvBlobs &Prev, std::vector<unsigned int> &matching, std::vector<unsigned int> &withinDIST)
{
}

void findClosest(cvb::CvBlobs &In, cvb::CvBlob *BLOB, unsigned int &closestBLOB, double &dist)
{
  cvb::CvBlobs::iterator It=In.begin();
  double MIN=65000;
  unsigned int MINID=0;
  while (It!=In.end())
  {
    cvb::CvBlob *cBlob=It->second;
    double DIST=cv::fast_abs(BLOB->centroid.x - cBlob->centroid.x) +
      cv::fast_abs(BLOB->centroid.y - cBlob->centroid.y);
    if(DIST < MIN)
    {
      dist=DIST;
      closestBLOB=It->first;
    }
    ++It;
  }
}

bool blobSizeIsRelevant(
    cvb::CvBlob *BLOB1,
    cvb::CvBlob *BLOB2,
    double ratio=LARVA_SIZE_COMPARISON_FACTOR)
{
		return (ratio*BLOB1->area > BLOB2->area &&
						((2-ratio)*BLOB1->area < BLOB2->area));
}

bool blobSizeIsRelevant(
    cvb::CvBlobs &In,
    cvb::CvBlob *BLOB,
    std::vector<unsigned int> &larvae,
    double ratio=LARVA_SIZE_COMPARISON_FACTOR)
{
  std::vector<unsigned int>::iterator IT=larvae.begin();
  double areaSUM=0;
  while (IT!=larvae.end())
  {
    areaSUM+=In[*IT]->area;
    ++IT;
  }
		return (ratio*BLOB->area > areaSUM &&
						((2-ratio)*BLOB->area < areaSUM ));
}

void findInDistance(cvb::CvBlobs &In, cvb::CvBlob *BLOB, double dist, std::map <unsigned int, double> neighbours)
{
}

void verbosePrint(std::stringstream &toPrint)
{
  if(LRVTRACK_VERBOSE_LEVEL>0)
    {
      std::cout << "LrvTrack DEBUG: " << toPrint.str() << std::endl;;
      toPrint.clear();
    }
}

void verbosePrint(const char * toPrint)
{
  if(LRVTRACK_VERBOSE_LEVEL>0)
    {
      std::cout << "LrvTrack DEBUG: " << toPrint << std::endl;;
    }
}
void verbosePrint(std::string &toPrint)
{
  if(LRVTRACK_VERBOSE_LEVEL>0)
    {
      std::cout << "LrvTrack DEBUG: " << toPrint << std::endl;;
    }
}

double getManhattanDistance(cvb::CvBlob *blob1, cvb::CvBlob *blob2)
{
  return (cv::fast_abs(blob1->centroid.x - blob2->centroid.x) +
          cv::fast_abs(blob1->centroid.y - blob2->centroid.y));
}

// Is blob1 closer to blob than blob2
// 0 equal
// 1 yes
// -1 no
int closerTo(cvb::CvBlob *blob1, cvb::CvBlob *blob2,cvb::CvBlob* blob)
{
  double diff = getManhattanDistance(blob1,blob) -
                getManhattanDistance(blob2,blob);
  if(diff > 0)
    return -1;
  else if (diff < 0)
    return 1;
  else
    return 0;
}

bool centresMatch(
    cvb::CvBlob *blob1,
    cvb::CvBlob *blob2,
    double factor=1-LARVA_SIZE_COMPARISON_FACTOR)
{
  double objectLength=
      std::max(blob1->maxx-blob1->minx,blob1->maxy-blob1->miny);
  if ((blob1->centroid.x - blob2->centroid.x)< factor*objectLength &&
      (blob1->centroid.y - blob2->centroid.y)< factor*objectLength )
    {
      return true;
    }
  else
    {
      return false;
    }
}

void assign_one(unsigned int preID,unsigned int postID)
{
  assignedPrevious[preID].push_back(postID);
  assignedNew[postID].push_back(preID);
  assignedPreMap[preID]=1;
}

void assign_one(unsigned int preID,
                std::vector<unsigned int> postID)
{
  std::vector<unsigned int>::iterator postIT=postID.begin();
  while(postIT!=postID.end())
  {
    assignedPrevious[preID].push_back(*postIT);
    assignedNew[*postIT].push_back(preID);
    ++postIT;
  }
  assignedPreMap[preID]=postID.size();
}

void assign_one(std::vector<unsigned int> preID,
                unsigned int postID)
{
  std::vector<unsigned int>::iterator preIT=preID.begin();
  while(preIT!=preID.end())
  {
    assignedNew[postID].push_back(*preIT);
    assignedNew[*preIT].push_back(postID);
    assignedPreMap[*preIT]=1;
    ++preIT;
  }
}

void assign_diverging(cvb::CvBlobs &New,
                      unsigned int CLUSTER_ID,
                      std::vector<unsigned int> &IDs
                      )
{
  // We have the following cases here:
  //  1) CLUSTER_ID Corresponds to no cluster: This means
  //     the cluster is newly found (it started as a cluster)
  //  2) CLUSTER_ID Corresponds to an existing cluster:
  //      - We need to understand how it split.
  //        M -> M?
  //        M -> N<M?
  //        This will give us a clue about which ones are left.
  std::map<unsigned int,std::vector<unsigned int> >::iterator dcIT;
  if((dcIT=detected_clusters.find(CLUSTER_ID))==detected_clusters.end())
  {
    //Not found new cluster NEW IDs to be given to the vector
    //elements
    assign_one(CLUSTER_ID,IDs);
    detected_larvae[CLUSTER_ID].isCluster=true;
  }
  else
  {
    std::vector<unsigned int> candidateLarvae(dcIT->second.begin()+1,
        dcIT->second.end());
    std::map<unsigned int, unsigned int> newAssignments;
    diverge_match_new(candidateLarvae,
                      IDs,
                      newAssignments,
                      New);
    std::vector<unsigned int>::iterator nIT=IDs.begin();
    std::vector<unsigned int> newCluster(dcIT->second.begin()+1,
        dcIT->second.end());
    std::vector<unsigned int> unassigned;
    while(nIT!=IDs.end())
    {
      if(newAssignments[*nIT]!=0)
      {
        assign_one(newAssignments[*nIT],*nIT);
        std::vector<unsigned int>::iterator eIT=newCluster.begin();
        while(eIT!=newCluster.end())
        {
          if (*eIT==newAssignments[*nIT])
            newCluster.erase(eIT);
          eIT++;
        }
      }
      else
      {
        unassigned.push_back(*nIT);
      }
      ++nIT;
    }
    if(unassigned.size()!=1)
    {
      verbosePrint("CASE OF TWO CLUSTERS EMERGING OUT OF ONE!!! Unhandled!!");
    }
    else
    {
      detected_clusters[unassigned[0]].push_back(newCluster.size());
      detected_clusters[unassigned[0]].insert(
          detected_clusters[unassigned[0]].end(),
          newCluster.begin(),
          newCluster.end());
    }
  }
}

void assign_clustering(bool NEW,
                      unsigned int CLUSTER_ID,
                      std::vector<unsigned int> &IDs
                      )
{
  std::vector<unsigned int>::iterator IT=IDs.begin();
  detected_clusters[CLUSTER_ID][0]=IDs.size();
  current_clusters[CLUSTER_ID][0]=IDs.size();
  std::map<unsigned int,std::vector<unsigned int> >::iterator dcIT;
  assign_one(IDs,CLUSTER_ID);
  while (IT!=IDs.end())
  {
    if((dcIT=detected_clusters.find(*IT))!=detected_clusters.end())
    {
      detected_clusters[CLUSTER_ID].insert(
          detected_clusters[CLUSTER_ID].end(),
          dcIT->second.begin(),
          dcIT->second.end());
      // -2 because the cluster ID should be removed and replaced
      // by it's elements (-1) and another -1 because the 
      // detected_clusters contains as first element the size 
      // estimate of the cluster (the candidates can be more than that)
      detected_clusters[CLUSTER_ID][0]+=dcIT->second.size()-2;
      detected_clusters.erase(dcIT);
    }
    else
    {
      detected_clusters[CLUSTER_ID].push_back(*IT);
      current_clusters[CLUSTER_ID].push_back(*IT);
    }
    ++IT;
  }
  current_clusters[CLUSTER_ID]=detected_clusters[CLUSTER_ID];
}

bool centresMatch(
    cvb::CvBlobs &In, 
    cvb::CvBlob *blob,
    std::vector<unsigned int> &larvae, 
    double factor=1-LARVA_SIZE_COMPARISON_FACTOR)
{
  double xcomb=0, ycomb=0;
  double objectLength=
      std::max(blob->maxx-blob->minx,blob->maxy-blob->miny);
  double lrvAreaSum=0;
  if(larvae.size()==1)
  {
    xcomb=In[larvae[0]]->centroid.x;
    ycomb=In[larvae[0]]->centroid.y;
  }
  else
  {
    std::vector<unsigned int>::iterator it=larvae.begin();
    while(it!=larvae.end())
    {
      xcomb+=In[*it]->centroid.x*In[*it]->area;
      ycomb+=In[*it]->centroid.y*In[*it]->area;
      lrvAreaSum+=In[*it]->area;
      ++it;
    }
    xcomb=xcomb/lrvAreaSum;
    ycomb=ycomb/lrvAreaSum;
  }
  if ((blob->centroid.x - xcomb)< factor*objectLength &&
      (blob->centroid.y - ycomb)< factor*objectLength )
    {
      return true;
    }
  else
    {
      return false;
    }
}

unsigned int factorial(unsigned int num)
{
  unsigned int val=1;
  while(num>1)
    val=(num--)*val;

  return val;
}

unsigned int kofn(unsigned int k, unsigned int n)
{
  return factorial(n)/(factorial(k)*factorial(n-k));
}

//This is not exactly the powerset.
//We do not need the sets of 1 element and the complete SET in there
void powersets(std::vector<unsigned int> &IN, std::vector<std::vector<unsigned int> > &OUT){
  for (unsigned int i=2 ; i<IN.size();i++)
  {
    std::vector<unsigned int> pointers;
    for(unsigned int k=0;k<i;k++)
    {
      pointers.push_back(k);
    }
    for (unsigned int j=0 ; j<kofn(i,IN.size());j++)
    {
      std::vector<unsigned int> cvec;
      for(unsigned int idx=0;idx<i;idx++)
      {
        cvec.push_back(IN[pointers[idx]]);
      }
      OUT.push_back(cvec);
      for(unsigned int inc=i;inc>0;inc--)
      {
        if(pointers[inc-1]<IN.size()-1-(i-inc))
        {
          pointers[inc-1]++;
          unsigned int add=0;
          for(unsigned int res=inc;res<i;res++)
          {
            add++;
            pointers[res]=pointers[inc-1]+add;
          }
          break;
        }
      }
    }
  }
}

int detect_diverging(std::vector<unsigned int> &preLarvaeNearby,
                       std::vector<unsigned int> &newLarvaeNearby,
                       cvb::CvBlobs &Pre,
                       cvb::CvBlobs &New)
{
  if(newLarvaeNearby.size()<=1)
    return -1; // No diverging is possible
  // Check the complete set first
  std::vector<unsigned int>::iterator pIT=preLarvaeNearby.begin();
  while(pIT!=preLarvaeNearby.end())
  {
    if(centresMatch(New,Pre[*pIT],newLarvaeNearby))
    {
      // Centres of all candidates match with new blob
      // cluster contains all. We can return
      assign_diverging(New,*pIT,newLarvaeNearby);
      continue;
    }
    std::vector<std::vector<unsigned int> > pSETS;
    powersets(newLarvaeNearby,pSETS);
    std::vector<std::vector<unsigned int> >::iterator pSIT=pSETS.begin();
    while(pSIT!=pSETS.end())
    {
      if(centresMatch(New,Pre[*pIT],*pSIT))
      {
        // Centres of all candidates match with new blob
        // cluster contains all. We can return
        assign_diverging(New,*pIT,*pSIT);
        break;
      }
      ++pSIT;
    }
    ++pIT;
  }
}

int detect_clustering(std::vector<unsigned int> &preLarvaeNearby,
                       std::vector<unsigned int> &newLarvaeNearby,
                       cvb::CvBlobs &Pre,
                       cvb::CvBlobs &New)
{
  if(preLarvaeNearby.size()<=1)
    return -1; // No clustering is possible

  // Check the complete set first
  std::vector<unsigned int>::iterator nIT=newLarvaeNearby.begin();
  while(nIT!=newLarvaeNearby.end())
  {
    if(centresMatch(Pre,New[*nIT],preLarvaeNearby))
    {
      // Centres of all candidates match with new blob
      // cluster contains all. We can return
      assign_clustering(true,++LARVAE_COUNT,preLarvaeNearby);
      continue;
    }
    std::vector<std::vector<unsigned int> > pSETS;
    powersets(preLarvaeNearby,pSETS);
    std::vector<std::vector<unsigned int> >::iterator pSIT=pSETS.begin();
    while(pSIT!=pSETS.end())
    {
      if(centresMatch(Pre,New[*nIT],*pSIT))
      {
        // Centres of all candidates match with new blob
        // cluster contains all. We can return
        assign_clustering(true,++LARVAE_COUNT,*pSIT);
        break;
      }
      ++pSIT;
    }
    ++nIT;
  }

}

// TODO:
// IMPORTANT NOTE!!
//  We have essentially three mappings:
//   1) The previous frame blobs -> new frame blobs
//   2) The new frame blobs -> previous frame blobs
//   3) The new frame blobs -> updated numbers
//
//   The updated numbers are:
//    a) the numbers that remained from the previous frame.
//    b) the numbers that diverged from frames before the previous
//    c) new numbers (increasing LARVAE_COUNT)
//
void newLarvaeTrack(cvb::CvBlobs &In, cvb::CvBlobs &Prev, cvb::CvBlobs &out)
{

  std::stringstream DEBUG;

  assignedNew.clear();
  assignedPrevious.clear();
  assignedPreMap.clear();

  cvb::CvBlobs::iterator prevIt=Prev.begin();
  while (prevIt!=Prev.end())
    {
			assignedPreMap[prevIt->first]=0;
			++prevIt;
		}

  verbosePrint("Starting Tracking loop");

	prevIt=Prev.begin();
  while (prevIt!=Prev.end())
  {
    //If the larvae from the previous frame has not been assigned
    //  -- It may be that it is assigned during a previous larva
    //     inspection.
    if(assignedPreMap[prevIt->first]==0)
    {
      cvb::CvBlob *preBlob=prevIt->second;
      unsigned int preID=prevIt->first;

      // Get larva around it before and after
      std::vector<unsigned int> preLarvaeNearby;
      std::vector<unsigned int> postLarvaeNearby;
      getNearbyLarvae(Prev,preBlob,preLarvaeNearby);
      getNearbyLarvae(In,preBlob,postLarvaeNearby);

      // Great case collection now:

      // CASE 1 -> 1
      if(preLarvaeNearby.size()==1 && postLarvaeNearby.size()==1)
      {
        if(blobSizeIsRelevant(In[postLarvaeNearby[0]],preBlob))
        {
          if(centresMatch(In,Prev[preLarvaeNearby[0]],postLarvaeNearby,1.5))
          {
            assign_one(preID,postLarvaeNearby[0]);
          }
          else
          {
            // Centres do NOT match! Maybe it disappeared.
            // TODO: Check
            verbosePrint("1-1: with size similar\
                and centres not matching. !UNHANDLED!");
          }
        }
        else
        {
          // Sizes are too different:
          //  Cases:
          //    1) Object dissapeared and another came into sight.
          //       Case will be handled when the other one is handled
          //    2) Merge with a large object.
          //       Case will be handled when the large object is handled.
          verbosePrint("1-1: with size too different. No assignment.");
        }
      }
      // END OF 1-1 CASE
      // GENERAL CASE:
      else
      {
        // Try to assign obvious ones:
        std::vector<unsigned int>::iterator preNearIt=preLarvaeNearby.begin();
        while(preNearIt!=preLarvaeNearby.end())
        {
          std::vector<unsigned int>::iterator postNearIt=postLarvaeNearby.begin();
          while(postNearIt!=postLarvaeNearby.end())
          {
            if (blobSizeIsRelevant(Prev[*preNearIt],In[*postNearIt]) &&
                centresMatch(Prev[*preNearIt],In[*postNearIt]))
            {
              //1-1 and one extra in sight
              //New Larva matches preID larva therefore assign and ignore the
              //other.
              assign_one(*preNearIt,*postNearIt);
              // erase the obvious assignments
              preLarvaeNearby.erase(preNearIt);
              postLarvaeNearby.erase(postNearIt);
              break;
            }
            ++postNearIt;
          }
          ++preNearIt;
        }
        // The rest are either appearing/disappearing/clustering/diverging
        // BUT if the main one was indeed assigned then we are in luck 
        // and we move on
        if(assignedPreMap[preID]>0)
          break;

        detect_clustering(preLarvaeNearby,postLarvaeNearby,Prev, In);

        if(assignedPreMap[preID]>0)
          break;

        detect_diverging(preLarvaeNearby,postLarvaeNearby,Prev, In);
        if(assignedPreMap[preID]>0)
          break;
        else
          verbosePrint("FOUND TESTCASE FOR MIXED DIVERGING/CONVERGING");
      }
    }
  }

  std::map<unsigned int, std::vector<unsigned int> >::iterator pIT=
    assignedPrevious.begin();
  while(pIT!=assignedPrevious.end())
  {
    pIT++;
  }
}

void findLarvaeInDistance(cvb::CvBlob *centerBlob, std::vector<unsigned int> &out, double distance, cvb::CvBlobs &set, unsigned int &minLabel, double &min)
{

  cvb::CvBlobs::iterator it=set.begin();
  double XVAL,YVAL,cur=0.0;
  while (it!=set.end())
    {
      cvb::CvBlob *blob;
      blob=((*it).second);
      // MIN_DISCARD_DISTANCE is used to quickly exclude larvae that are
      // too far to be considered
      if (((XVAL=cv::fast_abs(blob->centroid.x - centerBlob->centroid.x)) < MIN_DISCARD_DISTANCE ) &&
          ((YVAL=cv::fast_abs(blob->centroid.y - centerBlob->centroid.y)) < MIN_DISCARD_DISTANCE ))
        {
          cur=XVAL+YVAL;
          if (cur < distance)
            {
              out.push_back((*it).first);
            }
          if (cur < min)
            {
              min=cur;
              minLabel=(*it).first;
            }
        }
      ++it;
    }
}

void larvae_track(cvb::CvBlobs &In,cvb::CvBlobs &Prev,cvb::CvBlobs &out)
{

  // Map to keep track of the larvae that were tracked succesfully so far
  std::map<unsigned int,std::vector<unsigned int> > used_map;
  // Resetting current_clusters and current_diverged to empty
  current_clusters.clear();
  current_diverged.clear();
  current_new.clear();
  current_gone.clear();

  cvb::CvBlobs::iterator prevIt=Prev.begin();
  while (prevIt!=Prev.end())
    {
      // preBlob is the larva we are trying to match at this point.
      // preBlob belongs to the previous frame and we are trying to find
      // which larva of the current frame is more likely to be the one
      cvb::CvBlob *preBlob=(*prevIt).second;

      // Used to store the label of the larvae which is closest to the preBlob
      unsigned int minLabel;
      // Helper variables to store and compute distances
      double XVAL,YVAL,cur=0;
      // Stores the minimal distance found (actually the manhattan distance).
      double min=65535;

      // Here we check if the *prevIt was already assigned a match.
      // This can occur when diverging where the two diverging worms are
      // both assigned at the same time.

      if(out.find((*prevIt).first)!=out.end())
        {
          ++prevIt;
          continue;
        }
      // Here we start looking for the blob in "In" (i.e. the latest larvae set)
      // that is closest to the preBlob. Typically (in 15fps and up ) distances
      // of matching larvae should be quite small (0.0X for the centroids).

      std::vector<unsigned int> nearbyBlobs;
      findLarvaeInDistance(preBlob,nearbyBlobs,5.0,In,minLabel,min);
      /*    if(nearbyBlobs.size()==0)
            {

            findLarvaeInDistance(preBlob,nearbyBlobs,10.0,In,minLabel,min);
            if(nearbyBlobs.size()==1)
            {
      //Even for long jumps we have only one candidate so trick the
      //program to think it's closer.
      min=5;
      }
      }
      */
      // We went through all the current batch of larvae and found minimum distance min
      // for the larva with label minLabel
      if (min<=5)
        {
          // Short distance. Indicates:
          // 1) Most likely normal motion
          // 2) Joining of larvae without much change in the centroid (TODO)
          // 3) Spliting of larvae without much change in the centoid (TODO)
          out[(*prevIt).first]=In[minLabel];
          out[(*prevIt).first]->label=(*prevIt).first;
          used_map[minLabel].push_back((*prevIt).first);
          //inactive[(*prevIt).first]=0;
        }
      else if (min<100)
        {
          // Rather large jump. Indicates:
          // 1) larvae converged
          // 2) larvae diverged
          // 3) Big jump due to capture problems.
          // 4) TODO: Larva dissapeared and matched wrongly with another...

          // Check for converging
          // We check if more than one larvae match the new centroid and
          //  whether the middle of their centroids in the previous frame
          //  matches the centroid we have.

          // This is used only in the diverging case to point to the matching cluster.
          std::map< unsigned int, std::vector<unsigned int> >::iterator cluster;

          // If the larva was assigned previously to another one it may be that
          // we have a cluster! (used_map for the minLabel larvae has been assigned to
          // at least one larvae).
          if (used_map[minLabel].size()>0)
            {
              cvb::CvBlob *blob;
              blob=In[minLabel];
              // One or more larvae have matched to the new blob with ID minLabel
              // We will now test whether the centroids of the other matches give a
              //   good aproximation of the new centroid.
              //   TODO: Take care of case of more than 2 merging at one point.
              double newx = ((*prevIt).second->centroid.x +
                             Prev[used_map[minLabel][0]]->centroid.x);
              double newy = ((*prevIt).second->centroid.y +
                             Prev[used_map[minLabel][0]]->centroid.y);
              double XDIF,YDIF;
              // TODO: Use correct centroid calculation based on area not just middle
              if (((XDIF=cv::fast_abs(2*blob->centroid.x - newx)) < 10 ) &&
                  ((YDIF=cv::fast_abs(2*blob->centroid.y - newy)) < 10))
                {
                  // We're in luck! The center's match :) We have a cluster.
                  used_map[minLabel].push_back((*prevIt).first);
                  // Add both larvae to the detected clusters
                  // if any of the two larvae is already a cluster unite them in a new cluster
                  if(detected_larvae[used_map[minLabel][0]].isCluster)
                    {
                      if(detected_larvae[used_map[minLabel][1]].isCluster)
                        {
                          //if both larvae are clusters then:
                          //construct a composite cluster
                          detected_clusters[blob->label]=detected_clusters[used_map[minLabel][0]];
                          detected_clusters[blob->label].insert(
                            detected_clusters[blob->label].end(),
                            detected_clusters[used_map[minLabel][1]].begin(),
                            detected_clusters[used_map[minLabel][1]].end());
                          //delete the previous clusters
                          detected_clusters.erase(used_map[minLabel][0]);
                          detected_clusters.erase(used_map[minLabel][1]);
                        }
                      else
                        {
                          detected_clusters[blob->label]=detected_clusters[used_map[minLabel][0]];
                          detected_clusters[blob->label].push_back(used_map[minLabel][1]);
                        }
                    }
                  else if (detected_larvae[used_map[minLabel][1]].isCluster)
                    {
                      detected_clusters[blob->label]=detected_clusters[used_map[minLabel][1]];
                      detected_clusters[blob->label].push_back(used_map[minLabel][0]);
                    }
                  else
                    {
                      detected_clusters[blob->label].push_back(used_map[minLabel][0]);
                      detected_clusters[blob->label].push_back(used_map[minLabel][1]);
                    }
                  current_clusters[blob->label]=detected_clusters[blob->label];
                }
              else
                {
                  // May be that one larva dissapeared and matched wrongly with another nearby
                  // In this case we should only keep the one that is closest!
                  std::cerr << "PROBLEM!!! : Check how to resolve these cases where centers after do not match" << std::endl;
                }

            }
          // Check for diverging:
          // If the detected_clusters map contains a cluster for the
          // previous blob with the one we are examining we might be
          // looking into larvae that have separated.
          // TODO: nested clusters would also be nice...
          else if ((cluster=detected_clusters.find((*prevIt).first))!=detected_clusters.end())
            {
              cvb::CvBlobs::iterator it;
              it=In.begin();
              unsigned int secondDivergent=0;
              double XVAL,YVAL,cur=0;
              double min=40;
              while (it!=In.end())
                {
                  if ((*it).first!=minLabel)
                    {
                      cvb::CvBlob *blob;
                      blob=(*it).second;
                      if (((XVAL=cv::fast_abs(blob->centroid.x - preBlob->centroid.x)) < 90) &&
                          ((YVAL=cv::fast_abs(blob->centroid.y - preBlob->centroid.y)) < 90))
                        {
                          cur=XVAL+YVAL;
                          if (cur < min)
                            {
                              min=cur;
                              secondDivergent=(*it).first;
                            }
                        }
                    }
                  ++it;
                }
              if(secondDivergent==0)
                {
                  std::cerr << "PROBLEM!!! : Did not find second divergent." << std::endl;
                }
              else
                {
                  double newx = (In[minLabel]->centroid.x + In[secondDivergent]->centroid.x);
                  double newy = (In[minLabel]->centroid.y + In[secondDivergent]->centroid.y);
                  // TODO: Use correct centroid calculation based on area not just middle
                  if (((XVAL=cv::fast_abs(2*(*prevIt).second->centroid.x - newx)) < 20 ) &&
                      ((YVAL=cv::fast_abs(2*(*prevIt).second->centroid.y - newy)) < 20))
                    {
                      // We're in luck! The center's match :) We have a diverging cluster.
                      used_map[minLabel].push_back((*prevIt).first);
                      used_map[secondDivergent].push_back((*prevIt).first);
                      // Matching function here:
                      unsigned int Larva1=(*cluster).second[0];
                      unsigned int Larva2=(*cluster).second[1];
                      diverge_match(Larva1,Larva2,In[minLabel],In[secondDivergent]);
                      out[(*cluster).second[0]]=In[Larva1];
                      out[(*cluster).second[1]]=In[Larva2];
                      out[(*cluster).second[0]]->label=Larva1;
                      out[(*cluster).second[1]]->label=Larva2;
                      current_diverged[cluster->first].push_back(minLabel);
                      current_diverged[cluster->first].push_back(secondDivergent);
                      detected_clusters.erase(cluster);
                    }
                  else
                    {
                      std::cerr << "PROBLEM!!! : Check how to resolve these cases where centers do not match" << std::endl;
                    }
                }
            }
          else
            {
              // Here we have the following cases:
              // 1) First Larva Spotted belonging to a too large cluster
              //    -- Not the case since we have no unseen clusters
              // 2) Divergence of larvae that were clustered from the start <TODO>
              // 3) Big Jump due to dropped frame rate <TODO>
              //used_map[minLabel].push_back((*prevIt).first);
              //out[(*prevIt).first]=In[minLabel];
              //out[(*prevIt).first]->label=(*prevIt).first;

              used_map[minLabel].push_back((*prevIt).first);
              out[++LARVAE_COUNT]=In[minLabel];
              out[LARVAE_COUNT]->label=LARVAE_COUNT;
              current_new.push_back(LARVAE_COUNT);
            }
          //out[(*prevIt).first]=In[minLabel];
          //out[(*prevIt).first]->label=(*prevIt).first;
          //used_map[minLabel]++;
          //detected_clusters[minLabel].push_back((*prevIt).first);
        }
      ++prevIt;
    }

  cvb::CvBlobs::iterator it=In.begin();
  while (it!=In.end())
    {
      int ID=(*it).first;
      if (used_map.find(ID)==used_map.end())
        {
          out[++LARVAE_COUNT]=(*it).second;
          out[LARVAE_COUNT]->label=LARVAE_COUNT;
          current_new.push_back(LARVAE_COUNT);
        }
      ++it;
    }
  it=Prev.begin();
  while (it!=Prev.end())
    {
      int ID=(*it).first;
      if (out.find(ID)==out.end())
        {
          current_gone.push_back(ID);
        }
      ++it;
    }
  updateLarvae(out,Prev);

}

int handleArgs(int argc, char* argv[])
{
  char cDate[16];
  time_t curtime;
  struct tm *loctime;

  /* Get the current time. */
  curtime = time (NULL);

  /* Convert it to local time representation. */
  loctime = localtime (&curtime);

  strftime(cDate,16,"%Y%m%d_%H%M%S",loctime);

  std::string DATE=std::string(cDate);
  LRVTRACK_DATE=DATE;
  std::string VIDEO_NAME=DATE + VIDEO_TYPE;
  std::string PROCESSED_VIDEO_NAME=DATE + "_PROC"+VIDEO_TYPE;

  try
    {

      // allowed only on command line
      po::options_description generic("Generic options");
      generic.add_options()
      ("version,V", "Print version string")
      ("help,h", "Produce help message")
      ("list-cameras,l","List available cameras.")
      ("verbose,v",
       po::value<int>(&LRVTRACK_VERBOSE_LEVEL),
       "Verbosity (0-3).")
      ;

      po::options_description IOopts("Input Output Options");

      IOopts.add_options()

      ("results-folder,r",
       po::value<std::string>(&LRVTRACK_RESULTS_FOLDER)->implicit_value("."),
       "Main folder where the results will be recorded. The experiment folder for each experiment (DATE-SUFFIX) will be created under this folder.(Default: Current folder)"
      )

      ("file-suffix,f",
       po::value<std::string>(&LRVTRACK_NAME)->implicit_value(""),
       "Suffix to use for the experiment output folder and files. The experiment folder name will be DATE-<file-suffix>.(Default: "")"
      )

      ("save-video,s",
       po::value<std::string>
       (&LRVTRACK_SAVE_VIDEO)->implicit_value(DATE+VIDEO_TYPE),
       "Save original input as video. \
                         The option is disabled when the input is from a video \
                         source.(Default: <date-filesuffix>)"
      )

      ("save-tracked-video,p",
       po::value<std::string>
       (&LRVTRACK_SAVE_PROCESSED_VIDEO)->implicit_value(
         DATE+"_PROC"+VIDEO_TYPE),
       "Save resulting output as video with the filename \
                         specified.(Default: <date-filesuffix>_PROC)"
      )
      ("file-input,i",
       po::value<std::string>(&LRVTRACK_FILE_INPUT)->implicit_value(""),
       "Filename of video to be used as input for tracker."
      )

      ("camera-input,c",
       po::value<int>(&LRVTRACK_CAMERA_INPUT)->implicit_value(0),
       "Camera feed to be used as input for tracker. For now unfortunately \
                         one can only use numbers and no listing is possible."
      );

      po::options_description setupOpts("Setup Options");
      setupOpts.add_options()

      ("odour-cups,o",
       po::value<int> (&LRVTRACK_ODOUR_CUPS)->implicit_value(10),
       "When the flag is specified lrvTrack will look for odour cups. If specified\
                         with no arguments (e.g. -o) lrvTrack will look for as many as it can find.\
                         The number of odour cups can be specified as a value for this option."
      );

      po::options_description procOpts("Processing Options");
      procOpts.add_options()

      ("normalize-image,n",
       po::value<bool> (&LRVTRACK_NORMALIZE),
       "Set to true (default) we normalize the brightness values of \
                         each frame. This improves the accuracy of results but causes \
                         higher cpu utilization. (Default: True)"
      );

      po::options_description displayOpts("Display Options");
      displayOpts.add_options()

      ("show-skeleton,S",
       po::value<bool> (&LRVTRACK_SHOW_SKELETON),
       "Show skeleton for detected larvae (Default: False)"
      )

      ("show-contour,C",
       po::value<bool> (&LRVTRACK_SHOW_CONTOUR),
       "Show contour for detected larvae (Default: False)"
      )

      ("show-orientation,O",
       po::value<bool> (&LRVTRACK_SHOW_ORIENTATION),
       "Show orientation for detected larvae (Default: False)"
      )

      ("show-centroid,Z",
       po::value<bool> (&LRVTRACK_SHOW_CENTROID),
       "Show centroid for detected larvae (Default: True)"
      )

      ("show-head-tail,H",
       po::value<bool> (&LRVTRACK_SHOW_HEAD_TAIL),
       "Show head and tail for detected larvae (Default: True)"
      )

      ("show-tags,T",
       po::value<bool> (&LRVTRACK_SHOW_TAGS),
       "Show larvae tags (Default: True)"
      );

      po::options_description cmdline_options;
      cmdline_options.add(generic).add(IOopts).add(setupOpts).add(procOpts).add(displayOpts);

      po::variables_map vm;
      po::store(po::parse_command_line(argc, argv, cmdline_options), vm);
      po::notify(vm);

      if (vm.count("help"))
        {
          std::cout << cmdline_options << "\n";
          exit(1);
        }
      if (vm.count("results-folder"))
        {
          LRVTRACK_RESULTS_FOLDER=vm["results-folder"].as<std::string>();
        }
      else
        {
          LRVTRACK_RESULTS_FOLDER=".";
        }

      if(vm.count("odour-cups"))
        {
          LRVTRACK_ODOUR_CUPS=vm["odour-cups"].as<int>();
        }

      if (vm.count("file-suffix"))
        {
          LRVTRACK_NAME=DATE + "@" + LRVTRACK_NAME;
        }
      else
        {
          LRVTRACK_NAME=DATE;
        }

      if (vm.count("save-video"))
        {
          LRVTRACK_SAVE_VIDEO=LRVTRACK_NAME+VIDEO_TYPE;
        }
      else
        {
          LRVTRACK_SAVE_VIDEO="";
        }

      if (vm.count("save-tracked-video"))
        {
          LRVTRACK_SAVE_PROCESSED_VIDEO=LRVTRACK_NAME+"_PROC"+VIDEO_TYPE;
        }
      else
        {
          LRVTRACK_SAVE_PROCESSED_VIDEO="";
        }

      if (vm.count("list-cameras"))
        {

        }

      if (vm.count("file-input")<=0 && vm.count("camera-input")<=0)
        {
          std::cerr << "Error: No Input Specified. Exiting..." <<  std::endl;
          exit(2);
        }
      else if (vm.count("file-input")>0 && vm.count("camera-input")>0)
        {
          std::cerr << "Error: Ambiguous Input Specified. Exiting..." <<  std::endl;
          exit(2);
        }

      if (vm.count("file-input")>0)
        {
          LRVTRACK_FILE_INPUT=vm["file-input"].as<std::string>();
          if (LRVTRACK_FILE_INPUT=="")
            {
              std::cerr << "Error: Input file flag given but no file specified. Exiting..." << std::endl;
              exit(3);
            }
        }

      if (vm.count("camera-input")>0)
        {
          LRVTRACK_CAMERA_INPUT=vm["camera-input"].as<int>();
        }
      else
        {
          LRVTRACK_CAMERA_INPUT=-2;
        }
    }
  catch(...)
    {
      std::cerr << "Problem parsing options" << std::endl;
      exit(1);
    }

}

int setupTracker()
{
  if(!fs::is_directory(LRVTRACK_RESULTS_FOLDER))
    {
      try
        {
          fs::create_directory( LRVTRACK_RESULTS_FOLDER );
        }
      catch ( const std::exception & ex )
        {
          std::cerr << "LrvTrack warning: Could not create folder <"
                    << LRVTRACK_RESULTS_FOLDER
                    << ">. Got: "
                    << LRVTRACK_RESULTS_FOLDER
                    << " " << ex.what() << std::endl;
          exit(4);
        }
    }


  fs::path experimentFolder(LRVTRACK_RESULTS_FOLDER);
  experimentFolder=experimentFolder / LRVTRACK_DATE;
  if(!fs::is_directory(experimentFolder))
    {
      try
        {
          fs::create_directory(experimentFolder);
        }
      catch ( const std::exception & ex )
        {
          std::cerr << "LrvTrack warning: Could not create folder <"
                    << experimentFolder.string()
                    << ">. Got: "
                    << experimentFolder
                    << " "
                    << ex.what() << std::endl;
          exit(4);
        }
    }

}

void printSummary(cvb::CvBlobs &preBlobs,cvb::CvBlobs &blobs, bool first)
{

  if(first==true)
    {
      cpu_times elapsed(tS.elapsed());
      summary << "1 ";

      summary << std::left
              << std::fixed
              << std::setfill('0')
              << std::setprecision(3)
              << (double) elapsed.wall/1000000000.0
              << "  ";

      summary << blobs.size() << " ";

      summary << blobs.size() << " ";

      summary << "0.0" << "  ";
      summary << "0.00" << " ";
      summary << "0.000" << "  ";
      summary << "0.0" << " ";
      summary << "0.000" << "  ";
      summary << "0.0" << " ";
      summary << "0.000" << "  ";
      summary << "0.000" << " ";
      summary << "0.000" << "  ";
      summary << "0.000" << " ";
      summary << "0.000" << "  ";
      summary << "%% " ;
      cvb::CvBlobs::iterator it=blobs.begin();
      int i=1;
      while (it!=blobs.end())
        {
          summary << "0 " << i << " ";
          ++it;
          ++i;
        }
      summary << std::endl;
    }
  else
    {
      // output summary
      // Collect info

      //gettimeofday(&tC,0);

      //double elapsed=(tC.tv_sec - tS.tv_sec) + ((tC.tv_usec - tS.tv_usec)/1000000.0);
      cpu_times elapsed(tS.elapsed());
      double lifespanSUM=0;
      double speedSUM=0;
      double angularSpeedSUM=0;
      double lengthSUM=0;
      double lengthRelSUM=0;
      double widthSUM=0;
      double widthRelSUM=0;
      double aspectRatioSUM=0;
      double relAspectRatioSUM=0;
      double avgSizeSUM=0;
      double avgWiggleSUM=0;
      int larvaeToConsider=0;
      std::map<unsigned int,larvaObject>::iterator i=detected_larvae.begin();
      for (; i!=detected_larvae.end(); ++i)
        {
          larvaObject &cl=(*i).second;
          // avg duration
          if (cl.isCluster==true)
            break;

          ++larvaeToConsider;
          lifespanSUM+=(CURRENT_FRAME-cl.start_frame)/VIDEO_FPS;
          // avg speed
          float xvel,yvel;
          if(cl.centroid_speed_x.size()>0)
            xvel=cl.centroid_speed_x.back();
          else
            xvel=0.0;
          if(cl.centroid_speed_y.size()>0)
            yvel=cl.centroid_speed_y.back();
          else
            yvel=0.0;
          float sqvel=xvel*xvel+yvel*yvel;
          float vel[1];
          if(sqvel==0)
            vel[0]=0;
          else
            {
              ltsqrt(vel,&sqvel);
            }
          speedSUM+=vel[0];

          //avg angular speed
          angularSpeedSUM+=cl.angular_speed.back();

          //average length
          float length=cl.length.back();
          if(length!=0)
            {
              ltsqrt(vel,&length);
            }
          double realLength=vel[0];
          lengthSUM+=vel[0];
          lengthRelSUM+=cl.length.back()*cl.length.size()/cl.length_sum;

          //average width
          widthSUM+=cl.width.back();
          widthRelSUM+=cl.width.back()*cl.width.size()/cl.width_sum;

          //average aspect ratio
          aspectRatioSUM+=cl.width.back()/realLength;

          //average relative aspect ratio
          relAspectRatioSUM+=(cl.width.back()*cl.length_mean) /
                             (cl.width_mean*cl.length.back());
          //average wiggle
          double a1,a2;
          a1=angle(cl.heads.back(),
                   cl.lrvskels.back().Point20,
                   cl.tails.back());
          a2=angle(cl.heads.back(),
                   cl.lrvskels.back().Point80,
                   cl.tails.back());

          avgWiggleSUM += MAX(a1,a2);
          //average size
          avgSizeSUM+=cl.area.back();
        }

      std::setiosflags(std::ios::fixed);

      summary << CURRENT_FRAME-START_FRAME+1 << " " ;
      summary << std::left
              << std::fixed
              << std::setfill('0')
              << std::setprecision(3)
              << (double) elapsed.wall/1000000000.0
              << "  ";

      summary << detected_larvae.size() << " ";

      summary << detected_larvae.size() << " ";

      summary << std::left
              << std::fixed
              << std::setfill('0')
              << std::setprecision(1)
              << lifespanSUM/larvaeToConsider
              << "  ";

      //TODO angular speed
      summary << std::left
              << std::fixed
              << std::setfill('0')
              << std::setprecision(2)
              << speedSUM/larvaeToConsider
              << " ";

      summary << std::left
              << std::fixed
              << std::setfill('0')
              << std::setprecision(3)
              << angularSpeedSUM/larvaeToConsider
              << "  ";
      summary << std::left
              << std::fixed
              << std::setfill('0')
              << std::setprecision(1)
              << lengthSUM/larvaeToConsider
              << " ";

      summary << std::left
              << std::fixed
              << std::setfill('0')
              << std::setprecision(3)
              << lengthRelSUM/larvaeToConsider
              << "  ";

      summary << std::left
              << std::fixed
              << std::setfill('0')
              << std::setprecision(1)
              << widthSUM/larvaeToConsider
              << " ";

      summary << std::left
              << std::fixed
              << std::setfill('0')
              << std::setprecision(3)
              << widthRelSUM/larvaeToConsider
              << "  ";

      summary << std::left
              << std::fixed
              << std::setfill('0')
              << std::setprecision(3)
              << aspectRatioSUM/larvaeToConsider
              << " ";

      summary << std::left
              << std::fixed
              << std::setfill('0')
              << std::setprecision(3)
              << relAspectRatioSUM/larvaeToConsider
              << "  ";

      //TODO wiggle
      summary << std::left
              << std::fixed
              << std::setfill('0')
              << std::setprecision(3)
              << avgWiggleSUM/larvaeToConsider
              << " ";

      summary << std::left
              << std::fixed
              << std::setfill('0')
              << std::setprecision(3)
              << avgSizeSUM/larvaeToConsider
              << " ";

      if(current_new.size()>0 ||
          current_gone.size()>0 ||
          current_clusters.size()>0 ||
          current_diverged.size()>0)
        {
          summary << " %% ";
        }
      if(current_new.size()>0)
        {
          for (int i=0; i<current_new.size(); ++i)
            {
              summary << "0 " << current_new[i] << " ";
            }
        }
      if(current_gone.size()>0)
        {
          for (int i=0; i<current_gone.size(); ++i)
            {
              summary << current_gone[i] << " 0 ";
            }
        }
      if(current_clusters.size()>0)
        {
          std::map<unsigned int,std::vector<unsigned int> >::iterator ci;
          for (ci=current_clusters.begin(); ci!=current_clusters.end(); ++ci)
            {
              summary << ci->second[0] << " " << ci->first << " " ;
              summary << ci->second[1] << " " << ci->first << " " ;
            }
        }
      if(current_diverged.size()>0)
        {
          std::map<unsigned int,std::vector<unsigned int> >::iterator ci;
          for (ci=current_diverged.begin(); ci!=current_diverged.end(); ++ci)
            {
              summary << ci->first << " " << ci->second[0] << " " ;
              summary << ci->first << " " << ci->second[1] << " " ;
            }
        }
      summary << std::endl;
    }
}

void printBlobFile(larvaObject lrv)
{

  std::ostringstream BLOBFILENAME;
  std::ofstream blobFile;

  //double elapsed=(tC.tv_sec - tS.tv_sec) + ((tC.tv_usec - tS.tv_usec)/1000000.0);
  cpu_times elapsed(tS.elapsed());
  BLOBFILENAME <<
               LRVTRACK_NAME <<
               "." <<
               std::fixed <<
               std::setfill('0') <<
               std::setw(5) <<
               lrv.larva_ID <<
               ".dat";

  fs::path experimentFolderPath(LRVTRACK_RESULTS_FOLDER);
  experimentFolderPath = experimentFolderPath / LRVTRACK_DATE;

  fs::path blobFilePath(
    experimentFolderPath / BLOBFILENAME.str()
  );

  blobFile.open(blobFilePath.c_str());

  blobFile << CURRENT_FRAME-START_FRAME+1 << " " ;
  BLOBFILENAME << std::left
               << std::fixed
               << std::setfill('0')
               << std::setprecision(3)
               << elapsed.wall/1000000000LL
               << "  ";

  BLOBFILENAME << detected_larvae.size() << " ";
}

int main(int argc, char* argv[])
{
  bool SHOWTAGS=true;
  bool TRACK=false;
  bool STEP=true;
  bool showimg=true;

#ifdef LRVTRACK_WITH_CUDA
  cv::gpu::printShortCudaDeviceInfo(cv::gpu::getDevice());
#elif defined(LRVTRACK_WITH_OPENCL)
  std::vector<cv::ocl::Info> oclinfo;
  cv::ocl::getDevice(oclinfo);
  std::cout << "Found " << oclinfo.size() << " OpenCL devices. Using: ";
  std::cout << oclinfo[0].DeviceName[0] << std::endl;
#endif
  cv::VideoCapture capture;
  handleArgs(argc,argv);
  if(LRVTRACK_CAMERA_INPUT != -2)
    {
      //capture.open(CV_CAP_DC1394);
      capture.open(300);
      capture.set(CV_CAP_PROP_FRAME_WIDTH,1280);
      capture.set(CV_CAP_PROP_FRAME_HEIGHT,1024);
      capture.set(CV_CAP_PROP_FPS,24);
    }
  else if (LRVTRACK_FILE_INPUT != "")
    {
      capture.open(LRVTRACK_FILE_INPUT);
    }
  capture.set(CV_CAP_PROP_FORMAT,CV_8U);

  if (capture.get(CV_CAP_PROP_FPS)==0)
    {
      VIDEO_FPS=24.1;
    }
  else
    {
      if ( capture.get(CV_CAP_PROP_FPS) > 30 )
        {
          VIDEO_FPS=capture.get(CV_CAP_PROP_FPS)/2;
          std::cout << VIDEO_FPS << std::endl;
        }
      else
        {
          VIDEO_FPS=capture.get(CV_CAP_PROP_FPS);
          std::cout << VIDEO_FPS << std::endl;
        }
    }

  cv::Mat grey_bgFrame;
  unsigned int fg_threshold=10;
  unsigned int thresholdlow=90;
  unsigned int thresholdhigh=255;
  LARVAE_COUNT=0;

  if (!capture.isOpened())
    return 1;

  bool stop(false);
  cv::namedWindow("Extracted Frame");

  capture.read(bgFrame);

  cv::VideoWriter vidOut,vidPOut;

  if(bgFrame.channels()>1)
    {
      cvtColor(bgFrame,grey_bgFrame,CV_BGR2GRAY);
    }
  else
    {
      bgFrame.copyTo(grey_bgFrame);
    }

  //std::cout << bgFrame.rows << "," << bgFrame.cols << std::endl;
  std::vector<cv::Vec3f> circles;
  cv::HoughCircles(grey_bgFrame, circles, CV_HOUGH_GRADIENT,
                   1,   // accumulator resolution (size of the image / 2)
                   850,  // minimum distance between two circles
                   50, // Canny high threshold
                   180, // minimum number of votes
                   bgFrame.rows/3, bgFrame.rows/2); // min and max radiusV

  std::vector<cv::Vec3f> cups;
  if (LRVTRACK_ODOUR_CUPS>0 || LRVTRACK_ODOUR_CUPS==-1)
    {
      cv::Mat fgROI;
      fgROI=cv::Mat::zeros(grey_bgFrame.rows , grey_bgFrame.cols,grey_bgFrame.depth());
      if(circles.size()>0)
        {
          cv::circle(fgROI, cv::Point(circles[0][0],circles[0][1]),int(circles[0][2]/1.1),cv::Scalar(255),-1);
        }

      cv::Mat thr;

      thr=grey_bgFrame&fgROI;

      cv::morphologyEx(thr,thr,cv::MORPH_OPEN,cv::Mat(),cv::Point(-1,-1),5);
      cv::threshold(thr,
                    thr,
                    thresholdlow,
                    thresholdhigh,
                    cv::THRESH_BINARY|cv::THRESH_OTSU);

      cv::HoughCircles(thr, cups, CV_HOUGH_GRADIENT,
                       2,   // accumulator resolution (size of the image / 2)
                       1,  // minimum distance between two circles
                       240, // Canny high threshold
                       30, // minimum number of votes
                       3, 40); // min and max radiusV
    }

  std::vector<cv::Vec3f>::
  const_iterator itc= circles.begin();
  //std::cout << circles.size() << std::endl;
  if(circles.size()>0)
    {
      while (itc!=circles.end())
        {
          cv::circle(grey_bgFrame,
                     cv::Point((*itc)[0], (*itc)[1]), // circle centre
                     (*itc)[2]/1.1,       // circle radius
                     cv::Scalar(), // color
                     -1);              // thickness
          ++itc;
          break;
        }
    }

  fs::path experimentFolderPath(LRVTRACK_RESULTS_FOLDER);
  experimentFolderPath=experimentFolderPath / LRVTRACK_DATE;
  fs::path summaryFilePath(
    experimentFolderPath / (LRVTRACK_NAME + ".summary")
  );

  cvb::CvBlobs preBlobs;
  while (!stop)
    {
      // read next frame if any
      if (!capture.read(frame))
        break;

      if(TRACK)
        {
          if(START_FRAME==0)
            {
              setupTracker();
              if(!summary.is_open())
                {
                  summary.open(summaryFilePath.c_str());
                }
              tS.start();
              START_FRAME=CURRENT_FRAME;
            }

          if (LRVTRACK_SAVE_VIDEO!="" && !vidOut.isOpened())
            {
              vidOut.open((experimentFolderPath / LRVTRACK_SAVE_VIDEO).string(),
                          VIDEO_CODEC,
                          VIDEO_FPS,
                          bgFrame.size());
              if (!vidOut.isOpened())
                {
                  std::cerr << "Error opening video output file. Problems with video codec. Exiting..." << std::endl;
                  exit(5);
                }
            }

          if (LRVTRACK_SAVE_VIDEO!="")
            {
              vidOut << frame;
            }

          if (LRVTRACK_SAVE_PROCESSED_VIDEO!="" && !vidPOut.isOpened())
            {
              vidPOut.open((experimentFolderPath / LRVTRACK_SAVE_PROCESSED_VIDEO).string(),
                           VIDEO_CODEC,
                           VIDEO_FPS,
                           frame.size());
              if (!vidPOut.isOpened())
                {
                  std::cerr << "Error opening video output file. Problems with video codec. Exiting..." << std::endl;
                  exit(5);
                }
            }

          CURRENT_FRAME=capture.get(CV_CAP_PROP_POS_FRAMES);
          cv::Mat fg_frame;
          cv::Mat image;
          cv::Mat processed_frame;
          cv::Mat working_frame;
          cv::Mat fgROI;
          cv::Mat fg_image;
          cv::Mat mask;

          // Captured frame to BW (necessary for various filters and speed)
          if(frame.channels()>1)
            {
              cvtColor(frame,grey_frame,CV_BGR2GRAY);
            }
          else
            {
              frame.copyTo(grey_frame);
              cvtColor(grey_frame,frame,CV_GRAY2BGR);
            }

          //Background subtraction is done by eliminating areas outside of the
          //detected petri-dish
          // Here we subtract from the BW captured frame the background frame.
          //cv::addWeighted(grey_frame, 1.0, grey_bgFrame, -3.0, 0.0, fg_frame);
          cv::absdiff(grey_frame,grey_bgFrame,fg_frame);
          //cv::imshow("normWin",fg_frame);
          //fg_frame=grey_frame;
          fgROI=cv::Mat::zeros(fg_frame.rows , fg_frame.cols,fg_frame.depth());

          if(circles.size()>0)
            {
              cv::circle(fgROI, cv::Point(circles[0][0],circles[0][1]),int(circles[0][2]/0.99),cv::Scalar(255),-1);
              cv::circle(frame, cv::Point(circles[0][0],circles[0][1]),int(circles[0][2]/0.99),cv::Scalar(0,255,0),1);
            }

          if(cups.size()>0)
            {
              for( int i=0; i<cups.size(); ++i)
                {
                  cv::circle(fgROI, cv::Point(cups[i][0],cups[i][1]),int(cups[i][2]),cv::Scalar(0),-1);
                  cv::circle(frame, cv::Point(cups[i][0],cups[i][1]),int(cups[i][2]),cv::Scalar(0,0,255),1);
                }
            }

          fg_image=fg_frame&fgROI;

          cv::Mat fg_image_norm;

          lrvTrackNormalize(fg_image,fg_image_norm,0,255,cv::NORM_MINMAX);
          cv::threshold(fg_image,
                        thresholded_frame,
                        thresholdlow,
                        thresholdhigh,
                        cv::THRESH_BINARY|cv::THRESH_OTSU);
          //cv::blur(thresholded_frame,thresholded_frame,cv::Size(2,2));
          //cv::imshow("normWin",thresholded_frame);
          cvb::CvBlobs blobs;
          IplImage ipl_thresholded = thresholded_frame;
          labelImg=cvCreateImage(cvGetSize(&ipl_thresholded), IPL_DEPTH_LABEL, 1);

          unsigned int result=cvLabel(&ipl_thresholded, labelImg, blobs);
          cvb::cvFilterByArea(blobs, 32, 900);

          cvb::CvBlobs tracked_blobs;
          cvb::CvBlobs blob1;
          cvb::CvBlobs::iterator it=blobs.begin();

          if(preBlobs.size()>0)
            {
              FrameEllapsedTime = tP.elapsed();
              //larvae_track(blobs,preBlobs,tracked_blobs);
              newLarvaeTrack(blobs,preBlobs,tracked_blobs);
              tP.start();

              printSummary(preBlobs,blobs,false);

              //printBlobFile(detected_larvae[1]);

              preBlobs=tracked_blobs;
            }
          else
            {
              preBlobs.clear();
              //normalize blobs for next use
              cvb::CvBlobs::iterator it=blobs.begin();
              int i=1;
              while (it!=blobs.end())
                {
                  preBlobs[i]=(*it).second;
                  preBlobs[i]->label=i;
                  it++;
                  i++;
                }

              tP.start();

              printSummary(preBlobs,blobs,true);

              LARVAE_COUNT=preBlobs.size();
            }

          if (SHOWTAGS!=0)
            {
              it=tracked_blobs.begin();
              while (it!=tracked_blobs.end())
                {
                  std::stringstream sstm;
                  cvb::CvBlob *blob=(*it).second;
                  try
                    {
                      std::vector<unsigned int> &clusterMBs=detected_clusters.at((*it).first);
                      sstm << printVector(clusterMBs);
                      cv::putText(frame,
                                  sstm.str(),
                                  cv::Point2d(blob->centroid.x+12,blob->centroid.y+12),
                                  cv::FONT_HERSHEY_PLAIN,
                                  0.7,
                                  cv::Scalar(255,255,255),
                                  1,
                                  CV_AA);
                    }
                  catch (const std::out_of_range& oor)
                    {
                      //out of range. i.e. it's not there :)
                      sstm << (*it).first;
                      cv::putText(frame,
                                  sstm.str(),
                                  cv::Point2d(blob->centroid.x+12,blob->centroid.y+12),
                                  cv::FONT_HERSHEY_PLAIN,
                                  0.8,
                                  cv::Scalar(255,255,255),
                                  1,
                                  CV_AA);
                    }

                  // Here we define an extra PADDING
                  int PAD=2;
                  cv::Mat larvaROI(frame,
                                   cv::Rect(blob->minx-ROI_PADDING-PAD,
                                            blob->miny-ROI_PADDING-PAD,
                                            blob->maxx-blob->minx+2*(ROI_PADDING+PAD),
                                            blob->maxy-blob->miny+2*(ROI_PADDING+PAD))
                                  );

                  cv::circle(frame,
                             cv::Point2d(blob->centroid.x,blob->centroid.y),
                             1,
                             cv::Scalar(255,0,0),
                             -1);

                  if(detected_larvae[it->first].isCluster==false)
                    {
                      cv::circle(larvaROI,
                                 cv::Point2d(
                                   detected_larvae[it->first].lrvskels.back().Point20.x+PAD,
                                   detected_larvae[it->first].lrvskels.back().Point20.y+PAD),
                                 1,
                                 cv::Scalar(255,255,0),
                                 -1);

                      cv::circle(larvaROI,
                                 cv::Point2d(
                                   detected_larvae[it->first].heads.back().x+PAD,
                                   detected_larvae[it->first].heads.back().y+PAD),
                                 1,
                                 cv::Scalar(0,255,0),
                                 -1);

                      cv::circle(larvaROI,
                                 cv::Point2d(
                                   detected_larvae[it->first].tails.back().x+PAD,
                                   detected_larvae[it->first].tails.back().y+PAD),
                                 1,
                                 cv::Scalar(0,0,255),
                                 -1);

                      plotAngle(blob,larvaROI,PAD);
                    }
                  it++;
                }
            }



          cvReleaseImage(&labelImg);

        }
      //if (showimg)
      //{
      //cv::Mat flipframe;
      //cv::flip(frame,flipframe,0);
      //cv::imshow("Extracted Frame",flipframe);
      cv::imshow("Extracted Frame",frame);
      if (LRVTRACK_SAVE_PROCESSED_VIDEO!="")
        {
          vidPOut << frame;
        }

      //cv::displayStatusBar("Extracted Frame","FOO",0);
      // showimg=false;
      //}
      //else
      //  showimg=true;

      int k;

      if (STEP==true)
        k=46;
      else
        k=cv::waitKey(1);
      if (k>=0)
        {
          if (k==32)
            {
              while (cv::waitKey(1)!=32)
                {
                  //No OP
                }
            }
          if (k==46)
            {
              STEP=true;
              while ((k=cv::waitKey(1))>=127 || k<0)
                {
                }
              if (k!=46)
                STEP=false;
            }
          if (k==63234)
            {
              thresholdlow-=10;
              std::cerr << (int) thresholdlow << std::endl;
            }
          if (k==63235)
            {
              thresholdlow+=10;
              std::cerr << (int) thresholdlow << std::endl;
            }
          if (k==63232)
            {
              thresholdhigh+=10;
              std::cerr << (int) thresholdhigh << std::endl;
            }
          if (k==63233)
            {
              thresholdhigh-=10;
              std::cerr << (int) thresholdhigh << std::endl;
            }
          if (k==115)
            {
              SHOWTAGS=!SHOWTAGS;
            }
          if (k==116)
            {
              TRACK=!TRACK;
            }
          if (k==120)
            {
              exit(0);
            }
        }
    }

  summary.close();
  return 0;
}
