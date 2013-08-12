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

void verbosePrint(std::stringstream &toPrint)
{
  if(LRVTRACK_VERBOSE_LEVEL>0)
    {
      std::cout << "LrvTrack DEBUG: " << toPrint.str() << std::endl;;
      toPrint.str("");
      toPrint.clear();
    }
}

double avgVec(std::vector<double> vec)
{
  double SUM=0;
  std::vector<double>::iterator it=vec.begin();
  while(it!=vec.end())
  {
    SUM+=*it;
    ++it;
  }
  return SUM/vec.size();
}

double avgNVec(std::vector<int> vec)
{
  double SUM=0;
  unsigned int range;
  std::vector<int>::reverse_iterator it=vec.rbegin();
  if (vec.size()>=HISTORY_SIZE)
    range=HISTORY_SIZE;
  else
    range=vec.size();

  while(it!=vec.rbegin()+range)
  {
    SUM+=*it;
    ++it;
  }
  return SUM/vec.size();
}

bool centresMatch(
    cvb::CvBlob *blob1,
    cvb::CvBlob *blob2,
    double factor=LARVA_CENTRE_COMPARISON_FACTOR-1.0)
{
  std::stringstream DEBUG;
  double objectLength=
      std::min(std::max(blob1->maxx-blob1->minx,blob1->maxy-blob1->miny),
               std::max(blob2->maxx-blob2->minx,blob2->maxy-blob2->miny));

  DEBUG << "CentresMatch ["<< blob1->label << ", " << blob2->label << "]: Length: " << objectLength << " Difx: " << fabs(blob1->centroid.x - blob2->centroid.x) << " Dify: " << fabs(blob1->centroid.y - blob2->centroid.y) << " Threshold: " << factor*objectLength;
  verbosePrint(DEBUG);

  if (fabs(blob1->centroid.x - blob2->centroid.x)< factor*objectLength &&
      fabs(blob1->centroid.y - blob2->centroid.y)< factor*objectLength )
    {
      return true;
    }
  else
    {
      return false;
    }
}

void principalAxes(cvb::CvBlob &blob, 
                   cv::Point2f &major,
                   cv::Point2f &minor,
                   double &long_axis,
                   double &short_axis)
{
  // I is the smallest distance of the value of the pixel to
  // the thresholds high and low. J is a helper variable.
  short I;
  double a,b,f,L1,L2;
  double Sxx,Sxy,Syy,Si;
  cv::Mat cROI,ROI,lrvROI;

  try{
    cROI=grey_frame(cv::Rect(blob.minx,
          blob.miny,
          blob.maxx-blob.minx+1,
          blob.maxy-blob.miny+1)
        );
    cROI.copyTo(ROI);
    createLarvaContour(lrvROI, blob,CV_8UC1);
  }
  catch(...)
  {
    std::cerr << "principalAxes: Error creating contour" << std::endl;
    return;
  }
  //cv::Mat element = cv::getStructuringElement(cv::MORPH_CROSS, cv::Size(3, 3));
  //cv::dilate(lrvROI,lrvROI,element);
  //lrvTrackNormalize(lrvROI, lrvROI, 0, 255, CV_MINMAX );
  ROI=ROI&lrvROI;
  
  // Least squares fit (standard formula)
  Si = 0.0;
  Sxx = Sxy = Syy = 0.0;
  
  int nc=ROI.cols;
  int nr=ROI.rows;

  for(int i=0;i<nc;++i)
  {
    a=i-(blob.centroid.x-blob.minx);
    for(int j=0;j<nr;j++)
    {
      cv::Point p(i,j);
      uchar gdata= ROI.at<uchar>(p);
      uchar bwdata= lrvROI.at<uchar>(p);
      if (bwdata != 0)
      {
        I=gdata;
        if (255-I < I)
          I=255-I;
      }
      else
        continue;
      b=j-(blob.centroid.y-blob.miny);
      Si += I;
      Sxx += I*a*a;
      Sxy += I*a*b;
      Syy += I*b*b;
    }
  }

  Sxx /= Si;
  Sxy /= Si;
  Syy /= Si;
  f = sqrt( (Sxx-Syy)*(Sxx-Syy)+4*Sxy*Sxy )/2;
  a = 0.5*(Sxx+Syy);
  b = 0.5*(Sxx-Syy);
  L1 = a + f;
  L2 = a - f;
  // If major axis poorly defined for the shape, arbitrarily choose X direction.
  if (fabs(f-b) < 1e-6)
  {
    major = cv::Point2f( 1.0 , 0.0 );  
  }
  else
  {
    major = cv::Point2f(Sxy/(f-b) , 1.0);
    if (fabs(major.x) < 1e-6) 
      major.x = 0.0;

    //major = major*(1/(sqrt(major.x*major.x + major.y*major.y)));
    major.x=major.x/sqrt(major.x*major.x + major.y*major.y);
    major.y=major.y/sqrt(major.x*major.x + major.y*major.y);
  }
  // Likewise for minor axis
  if (fabs(b+f) < 1e-7)
  {
    minor = cv::Point2f( 1.0 , 0.0 );
  }
  else
  {    
    minor = cv::Point2f(Sxy/(b+f) , -1.0);
    if (fabs(minor.x) < 1e-6) 
      minor.x = 0.0;
    minor.x = minor.x/sqrt(minor.x*minor.x + minor.y*minor.y);
    minor.y = minor.y/sqrt(minor.x*minor.x + minor.y*minor.y);
  }

  // Find length along the least squares axes
  cv::Point2f p;
  double maj_max,maj_min,min_max,min_min;
  maj_max = maj_min = min_max = min_min = 0.0;

  for(int i=0;i<nc;++i)
  {
    int j=0;
    cv::Point p;
    double gval=0.0;
    while(gval==0 && j<nr)
    {
      p=cv::Point(i,j);
      gval=lrvROI.at<uchar>(p);
      j++;
    }
    if(j==nr)
      continue;
    
    p.x = i - (blob.centroid.x - blob.minx);
    p.y = j - (blob.centroid.y - blob.miny);

    a = p.x*major.x+p.y*major.y;
    b = p.x*minor.x+p.y*minor.y;
    if (a > maj_max) maj_max = a;
    else if (a < maj_min) maj_min = a;
    if (b > min_max) min_max = b;
    else if (b < min_min) min_min = b;
  }
  
  // Set data: axes (of standard deviation length) and length of object along two axes
  major *= sqrt(L1);
  minor *= sqrt(L2);
  long_axis = maj_max - maj_min;
  short_axis = min_max - min_min;
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

double avgNVec(std::vector<double> vec)
{
  double SUM=0;
  unsigned int range;
  std::vector<double>::reverse_iterator it=vec.rbegin();
  if (vec.size()>=HISTORY_SIZE)
    range=HISTORY_SIZE;
  else
    range=vec.size();

  while(it!=vec.rbegin()+range)
  {
    SUM+=*it;
    ++it;
  }
  return SUM/vec.size();
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


std::string printUIMap(std::map<unsigned int, unsigned int> &A)
{
  std::stringstream F;
  std::map<unsigned int, unsigned int>::iterator Ait=A.begin();
  F << "[ ";
  while(Ait!=A.end())
  {
    F << Ait->first << " -> " << Ait->second << ",";
    ++Ait;
  }
  F << " ]";
  return F.str();
}


void findHeadTail(std::vector<cv::Point2f> &startPoints,
                  larvaObject &lrv,
                  cv::Point2f &Head,
                  cv::Point2f &Tail,
                  bool force_SurroundingValSearch=false)
{
  int max=0,min=65535;
  unsigned int i=0;
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
      lrv.roundness.back()<=2.0 ||
      (lrv.heads.back().x == 0 && lrv.heads.back().y==0 &&
       lrv.tails.back().x == 0 && lrv.tails.back().y==0)
    )
    {
      for (i=0; i<startPoints.size(); ++i)
        {
          cv::Point2f p=startPoints[i];
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

void updateOneLarva(cvb::CvBlobs &In,
                    cvb::CvBlobs &Prev,
                    cvb::CvBlobs::iterator it,
                    tbb::concurrent_hash_map<unsigned int, larvaObject> &NEW)
{
  unsigned int ID=(*it).first;
  cvb::CvBlob blob=*((*it).second);
  cv::Mat larvaROI,cntPoints;
  createLarvaContour(larvaROI,blob);
  createLarvaContourPoints(cntPoints,blob);

  std::map<unsigned int,larvaObject>::iterator curLarva;
  // NEW LARVA OBJECT!
  if ((curLarva=detected_larvae.find(ID))==detected_larvae.end())
  {
    // Create and allocate the new object
    larvaObject newLarva;

    // Set the frame of it's existence
    newLarva.start_frame=CURRENT_FRAME;

    // Give the larva the necessary ID
    newLarva.larva_ID=ID;

    // Add the blob of the larva to its blob history
    newLarva.blobs.push_back(blob);

    // State that the larva is not in a blob
    newLarva.inCluster.push_back(false);

    cv::Point2f centroid=cv::Point2f(
        (blob.centroid.x-blob.minx+ROI_PADDING),
        (blob.centroid.y-blob.miny+ROI_PADDING)
        );

    cv::Point2f centroidf=cv::Point2f(
        (blob.centroid.x),
        (blob.centroid.y)
        );

    double FrameEllapsedSeconds=FrameEllapsedTime.wall/1000000000.0;
    newLarva.capture_times.push_back(FrameEllapsedSeconds);

    // Initialize the speed to 0
    newLarva.midpoint_speed_x.push_back(0);
    newLarva.midpoint_speed_y.push_back(0);

    newLarva.centroids.push_back(centroid);
    newLarva.centroids_full.push_back(centroidf);

    if (detected_clusters.find(ID)==detected_clusters.end())
    {

      ++newLarva.lifetimeWithStats;
      newLarva.lastBlobWithStats=0;

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
      newLarva.centroid_distance_x.push_back(0);
      newLarva.centroid_distance_y.push_back(0);
      newLarva.centroid_distance_x_sum=0;
      newLarva.centroid_distance_y_sum=0;

      newLarva.roundness.push_back((perimeter*perimeter)/(2*CV_PI*blob.area));

      // In this block we compute the inner spine for the larva
      std::vector<cv::Point2f> cntPoints;
      blobToPointVector(blob,cntPoints);
      larvaDistanceMap Distances(cntPoints);
      computeSpine(blob,Distances);
      newLarva.lrvDistances.push_back(Distances);
      newLarva.length.push_back(Distances.MaxDist);
      newLarva.length_mean = Distances.MaxDist;
      newLarva.length_sum= Distances.MaxDist;
      newLarva.length_max = Distances.MaxDist;
      newLarva.length_min = Distances.MaxDist;
      PointPair MAXPair=Distances.MaxDistPoints;

      cv::Point2f Head,Tail;
      std::vector<cv::Point2f> startPoints;
      startPoints.push_back(cv::Point2f(
            MAXPair.first.x-blob.minx+ROI_PADDING,
            MAXPair.first.y-blob.miny+ROI_PADDING)
          );
      startPoints.push_back(cv::Point2f(
            MAXPair.second.x-blob.minx+ROI_PADDING,
            MAXPair.second.y-blob.miny+ROI_PADDING)
          );
      newLarva.angular_speed.push_back(0);

      findHeadTail(startPoints,newLarva,Head,Tail);
      newLarva.heads.push_back(Head);
      newLarva.tails.push_back(Tail);

      cv::Point2f MP;
      MP.x=Distances.MidPoint.x-newLarva.blobs.back().minx;
      MP.y=Distances.MidPoint.y-newLarva.blobs.back().miny;
      cv::Point2f AxS(MP.x,Tail.y);
      newLarva.headBodyAngle.push_back(angleD(Head,MP,Tail));
      newLarva.orientationAngle.push_back(cvb::cvAngle(&blob));

      newLarva.width.push_back(Distances.WidthDist);
      newLarva.width_mean = Distances.WidthDist;
      newLarva.width_sum= Distances.WidthDist;
      newLarva.width_max = Distances.WidthDist;
      newLarva.width_min = Distances.WidthDist;


      if(DEBUG_INFO!=0)
      {
        std::cout << CURRENT_FRAME <<
          " , " << newLarva.larva_ID <<
          " , " << newLarva.headBodyAngle.back() <<
          " , " << newLarva.orientationAngle.back() <<
          " , " << newLarva.midpoint_speed_x.back() <<
          " , " << newLarva.midpoint_speed_y.back() <<
          " , " << newLarva.centroids.size()-1 <<
          " , " << newLarva.centroids.back().x + blob.minx <<
          " , " << newLarva.centroids.back().y + blob.miny <<
          " , " << newLarva.inCluster.size()-1 <<
          " , " << newLarva.inCluster.back() <<
          " , " << newLarva.isCluster <<
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
      for ( cl=detected_clusters[ID].begin()+1 ; cl!=detected_clusters[ID].end() ; ++cl)
      {
        larvaObject &clusterLarva=detected_larvae[*cl];
        clusterLarva.inCluster.push_back(ID);
        clusterLarva.blobs.push_back(blob);
        if(DEBUG_INFO!=0)
        {
          std::cout << CURRENT_FRAME <<
            " , " << clusterLarva.larva_ID <<
            " , " << clusterLarva.headBodyAngle.back() <<
            " , " << clusterLarva.orientationAngle.back() <<
            " , " << clusterLarva.midpoint_speed_x.back() <<
            " , " << clusterLarva.midpoint_speed_y.back() <<
            " , " << clusterLarva.centroids.size()-1 <<
            " , " << clusterLarva.centroids.back().x + blob.minx <<
            " , " << clusterLarva.centroids.back().y + blob.miny <<
            " , " << clusterLarva.inCluster.size()-1 <<
            " , " << clusterLarva.inCluster.back() <<
            " , " << clusterLarva.isCluster <<
            std::endl;
        }
      }

    }
    //detected_larvae[ID]=newLarva;
    //NEW[ID]=newLarva;
    tbb::concurrent_hash_map<unsigned int,larvaObject>::accessor a;
    NEW.insert(a,ID);
    a->second=newLarva;
  }
  // UPDATED LARVA OBJECT!
  else
  {
    //Reference to current larva
    larvaObject &cur_larva=(*curLarva).second;
    //Pointer for the previous blob
    cvb::CvBlob &preBlob = cur_larva.blobs.back();

    // Set the ID of the larvaObject to the ID found TODO:Probably unnecessary
    cur_larva.larva_ID=ID;

    // Add the current blob to the blobs history of the larva
    cur_larva.blobs.push_back(blob);

    // Create the skeleton of the larva and add it to the skeletons history
    cv::Point2f centroid=cv::Point2f(
        (blob.centroid.x-blob.minx+ROI_PADDING),
        (blob.centroid.y-blob.miny+ROI_PADDING)
        );

    cv::Point2f centroidf=cv::Point2f(
        (blob.centroid.x-blob.minx+ROI_PADDING),
        (blob.centroid.y-blob.miny+ROI_PADDING)
        );

    cur_larva.centroids_full.push_back(centroidf);
    cur_larva.centroids.push_back(centroid);

    double FrameEllapsedSeconds=FrameEllapsedTime.wall/1000000000.0;
    cur_larva.capture_times.push_back(CurrentTime.wall/1000000000.0);

    // Look if larva is a blob
    if ( (!cur_larva.isCluster) &&
        (detected_clusters.find(ID)==detected_clusters.end()))
    {

      // If not then:
      //  Update area values for larva.

      ++cur_larva.lifetimeWithStats;
      cur_larva.lastBlobWithStats=cur_larva.blobs.size()-1;
      cur_larva.isCluster=false;

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
      //larvaSkel newLarvaSkel(larvaROI,centroid);
      //cur_larva.lrvskels.push_back(newLarvaSkel);

      cur_larva.centroid_distance_x.push_back(fabs(blob.centroid.x - preBlob.centroid.x));
      cur_larva.centroid_distance_y.push_back(fabs(blob.centroid.y - preBlob.centroid.y));
      cur_larva.centroid_distance_x_sum+=fabs(blob.centroid.x - preBlob.centroid.x);
      cur_larva.centroid_distance_x_sum+=fabs(blob.centroid.y - preBlob.centroid.y);


      // Point coordinates for head/tail
      cv::Point2f Head,Tail;

      // Map to keep the distances of each point to all the others
      // Pair of points to keep the points with the Maximum distance (i.e. head and tail :) )
      std::vector<cv::Point2f> cntPoints;
      blobToPointVector(blob,cntPoints);
      larvaDistanceMap Distances(cntPoints);
      // Compute all the inner distances for the larva
      //computeInnerDistances(blob,Distances,newLarvaSkel.MidPoint);
      computeSpine(blob,Distances);

      if(cur_larva.inCluster.back()==0)
      {
        cur_larva.midpoint_speed_x.push_back(
            (Distances.MidPoint.x - cur_larva.lrvDistances.back().MidPoint.x)
            /FrameEllapsedSeconds);
        cur_larva.midpoint_speed_y.push_back(
            (Distances.MidPoint.y - cur_larva.lrvDistances.back().MidPoint.y)
            /FrameEllapsedSeconds);
      }
      else
      {
        cur_larva.midpoint_speed_x.push_back(
            (Distances.MidPoint.x - cur_larva.lrvDistances.back().MidPoint.x)
            /(CurrentTime.wall - cur_larva.capture_times.back())
            );
        cur_larva.midpoint_speed_y.push_back(
            (Distances.MidPoint.y - cur_larva.lrvDistances.back().MidPoint.y)
            /(CurrentTime.wall - cur_larva.capture_times.back())
            );
      }

      cur_larva.lrvDistances.push_back(Distances);
      cur_larva.length.push_back(Distances.MaxDist);
      cur_larva.length_mean=(cur_larva.length_mean+Distances.MaxDist)/2;
      cur_larva.length_sum=cur_larva.length_sum+Distances.MaxDist;
      if (cur_larva.length_max < Distances.MaxDist)
      {
        cur_larva.length_max=Distances.MaxDist;
      }
      if (cur_larva.length_min > Distances.MaxDist)
      {
        cur_larva.length_min=Distances.MaxDist;
      }
      PointPair MAXPair=Distances.MaxDistPoints;
      
      cur_larva.width.push_back(Distances.WidthDist);
      cur_larva.width_mean=(cur_larva.width_mean+Distances.WidthDist)/2;
      cur_larva.width_sum=cur_larva.width_sum+Distances.WidthDist;
      if (cur_larva.width_max < Distances.WidthDist)
      {
        cur_larva.width_max=Distances.WidthDist;
      }
      if (cur_larva.width_min > Distances.WidthDist)
      {
        cur_larva.width_min=Distances.WidthDist;
      }

      double greyVal=getGreyValue(larvaROI,blob,grey_frame);
      cur_larva.grey_value.push_back(greyVal);
      cur_larva.grey_value_mean=(cur_larva.grey_value_mean+greyVal)/2;
      cur_larva.grey_value_sum=cur_larva.grey_value_sum+greyVal;
      if (cur_larva.grey_value_max < greyVal)
      {
        cur_larva.grey_value_max=greyVal;
      }
      if (cur_larva.grey_value_min > greyVal)
      {
        cur_larva.grey_value_min=greyVal;
      }

      double perimeter=getPerimeter(blob);
      cur_larva.perimeter.push_back(perimeter);
      cur_larva.perimeter_mean=(cur_larva.perimeter_mean+perimeter)/2;
      cur_larva.perimeter_sum=cur_larva.perimeter_sum+perimeter;
      if (cur_larva.perimeter_max < perimeter)
      {
        cur_larva.perimeter_max=perimeter;
      }
      if (cur_larva.perimeter_min > perimeter)
      {
        cur_larva.perimeter_min=perimeter;
      }

      cur_larva.roundness.push_back((perimeter*perimeter)/(2*CV_PI*blob.area));

      // Construct a vector of points including both
      // head and tail to decide which is which
      std::vector<cv::Point2f> startPoints;

      // Points must be corrected to match the ROI including the padding set
      startPoints.push_back(cv::Point2f(
            MAXPair.first.x-blob.minx+ROI_PADDING,
            MAXPair.first.y-blob.miny+ROI_PADDING)
          );
      startPoints.push_back(cv::Point2f(
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

      cv::Point2f MP;
      MP.x=Distances.MidPoint.x-cur_larva.blobs.back().minx;
      MP.y=Distances.MidPoint.y-cur_larva.blobs.back().miny;
      cur_larva.headBodyAngle.push_back(angleD(Head,MP,Tail));
      cv::Point2f AxS(MP.x,Tail.y);
      cur_larva.headBodyAngle.push_back(angleD(Head,MP,Tail));
      cur_larva.orientationAngle.push_back(cvb::cvAngle(&blob));

      double curAngle=cvb::cvAngle(&blob);
      double preAngle=cvb::cvAngle(&preBlob);

      cur_larva.angular_speed.push_back(cv::fast_abs(curAngle-preAngle)/FrameEllapsedSeconds);

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
          " , " << cur_larva.headBodyAngle.back() <<
          " , " << cur_larva.orientationAngle.back() <<
          " , " << cur_larva.midpoint_speed_x.back() <<
          " , " << cur_larva.midpoint_speed_y.back() <<
          " , " << cur_larva.centroids.size()-1 <<
          " , " << cur_larva.centroids.back().x + blob.minx <<
          " , " << cur_larva.centroids.back().y + blob.miny <<
          " , " << cur_larva.inCluster.size()-1 <<
          " , " << cur_larva.inCluster.back() <<
          " , " << cur_larva.isCluster <<
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
      if(detected_clusters[ID].size()<2)
      {
        std::cerr << "We are working with a cluster but the size of registered nodes is too small" << std::endl;
        return;
      }
      
      for ( cl=detected_clusters[ID].begin()+1 ; cl!=detected_clusters[ID].end() ; ++cl)
      {
        larvaObject &clusterLarva=detected_larvae[*cl];
        clusterLarva.inCluster.push_back(ID);
        clusterLarva.blobs.push_back(blob);
        if(DEBUG_INFO!=0)
        {
          std::cout << CURRENT_FRAME <<
            " , " << clusterLarva.larva_ID <<
            " , " << clusterLarva.headBodyAngle.back() <<
            " , " << clusterLarva.orientationAngle.back() <<
            " , " << clusterLarva.midpoint_speed_x.back() <<
            " , " << clusterLarva.midpoint_speed_y.back() <<
            " , " << clusterLarva.centroids.size()-1 <<
            " , " << clusterLarva.centroids.back().x + blob.minx <<
            " , " << clusterLarva.centroids.back().y + blob.miny <<
            " , " << clusterLarva.inCluster.size()-1 <<
            " , " << clusterLarva.inCluster.back() <<
            " , " << clusterLarva.isCluster <<
            std::endl;
        }
      }
    }
  }
}

class larvaeUpdateBody : public cv::ParallelLoopBody
{
private:
  cvb::CvBlobs &In;
  cvb::CvBlobs &Prev;
  tbb::concurrent_hash_map<unsigned int, larvaObject> &NEW;
  cvb::CvBlobs::iterator it;

public:
  larvaeUpdateBody(cvb::CvBlobs &IIn, 
      cvb::CvBlobs &IPrev,
      tbb::concurrent_hash_map<unsigned int, larvaObject> &n
      ): 
    In(IIn),
    Prev(IPrev),
    NEW(n)
  {}
    void operator ()(const cv::Range& range) const
    {
        for (int i = range.start; i < range.end; ++i)
        {
          cvb::CvBlobs::iterator it=In.begin();
          std::advance(it,i);
            updateOneLarva(In,Prev,it,NEW);
        }
    }
};


void updateLarvae(cvb::CvBlobs &In, cvb::CvBlobs &Prev)
{
  cvb::CvBlobs::iterator it=In.begin();
  tbb::concurrent_hash_map<unsigned int, larvaObject> NEW;
  larvaeUpdateBody body(In,Prev,NEW);
  //cv::parallel_for_(cv::Range(0, In.size()), body);
  
  for(;it!=In.end();++it)
  {
    updateOneLarva(In,Prev,it,NEW);
  }

  tbb::concurrent_hash_map<unsigned int, larvaObject>::iterator nit=
    NEW.begin();
  while (nit!=NEW.end())
  {
    detected_larvae[nit->first]=nit->second;
    ++nit;
    }
}

double is_larva(cvb::CvBlob *blob)
{
  std::stringstream DEBUG;
  std::vector<cv::Point2f> cntPointsF;
  std::vector<cv::Point> cntPoints;
  std::vector<cv::Point> SimplePoints;
  std::vector<cv::Point> hull;
  std::vector<int> hullPoints;
  double defectSUM=0;
  blobToPointVector(*blob,cntPointsF);

  std::vector<cv::Point2f>::iterator f=cntPointsF.begin();
  while(f!=cntPointsF.end())
  {
    cv::Point i;
    i.x=(int) f->x;
    i.y=(int) f->y;
    cntPoints.push_back(i);
    ++f;
  }

  cv::approxPolyDP(cntPoints,SimplePoints,0.9,true);
  cv::convexHull(SimplePoints,hullPoints);
  cv::convexHull(SimplePoints,hull);
  
  std::vector<cv::Vec4i> defects;
  convexityDefects(SimplePoints,hullPoints, defects);

  for (unsigned int i=0;i<defects.size();i++)
  {
    defectSUM+=defects[i][3];
  }

  double ret=defectSUM * (0.5*defects.size() * blob->area/80 );
  if ( ret > 1300 )
  {
    DEBUG << "ISLARVA: Blob [" << blob->label << "] is likely a blob (ret: " << ret << ")"; 
    verbosePrint(DEBUG);
  }
  else if (ret < 900 )
  {
    DEBUG << "ISLARVA: Blob [" << blob->label << "] is likely a larva (ret: " << ret << ")"; 
    verbosePrint(DEBUG);
  }
  else
  {
    DEBUG << "ISLARVA: Blob [" << blob->label << "] is a bit vague... (ret: " << ret << ")"; 
    verbosePrint(DEBUG);
  }
  return ret; 
}

inline double SQUARE(double n)
{
  return n*n;
}

//Currently only for the 4 sized vector
inline void vmin(std::vector<double>& vals, std::vector<int>& seq)
{
  for (unsigned int i=0 ; i<vals.size() ; ++i )
    {
      int sum=0;
      for (unsigned int j=0; j<vals.size() ; ++j)
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
    std::string printVector(std::vector<number> vec,int position=0)
    {
      if (vec.size()==0)
        return "";
      std::stringstream sstm;
      //bool const is_number= std::is_arithmetic<number>::value;
      //static_assert( is_number, "Provided type is not an arithmetic type");
      sstm << "[" ;
      typename std::vector<number>::const_iterator i=vec.begin()+position;
      sstm << *i ;
      ++i;
      for( ; i != vec.end(); ++i)
        sstm << ","<< *i;
      sstm << "]";
      return sstm.str();
    }
}

double calculate_MH_sum(std::map<unsigned int, unsigned int> &AM,
  std::map<unsigned int, cv::Mat> &candidateCovarMat,
  std::map<unsigned int, cv::Mat> &candidateMeanMat,
  std::map<unsigned int, cv::Mat> &newMeanMat,
  cvb::CvBlobs &NEW
                        )
{
  /*
  std::cout << "Duration threshold: " << COLLISION_DURATION_THRESHOLD << std::endl;
  */
  double SUM=0.0;
  std::map <unsigned int, unsigned int>::iterator mIT=AM.begin();
  while(mIT!=AM.end())
  {
    double speed_x=
      detected_larvae[mIT->second].blobs.back().centroid.x-
      NEW[mIT->first]->centroid.x/
      (CurrentTime.wall-detected_larvae[mIT->second].capture_times.back());
    double speed_y=
      detected_larvae[mIT->second].blobs.back().centroid.y-
      NEW[mIT->first]->centroid.y/
      (CurrentTime.wall-detected_larvae[mIT->second].capture_times.back());
    cv::Mat nm;
    newMeanMat[mIT->first].copyTo(nm);
    nm.push_back(speed_x);
    nm.push_back(speed_y);
    SUM+=cv::Mahalanobis(nm,
        candidateMeanMat[mIT->second],
        candidateCovarMat[mIT->second]);
    /*cv::Mat res;
    res=newMeanMatInst-candidateMeanMat[mIT->second];
    cv::Scalar r=cv::sum(res);
    SUM+=fabs(r[0]);*/

    /*
    std::cout << "CALCULATE_MH_SUM: " << newMeanMat[mIT->first] << std::endl;
    std::cout << "CALCULATE_MH_SUM: " << candidateMeanMat[mIT->second] << std::endl;
    std::cout << "SUM: " << SUM << std::endl;
    */
    ++mIT;
  }
  return SUM;
}


void assign_combinations(std::vector<unsigned int> &NEW,
    std::map <unsigned int, unsigned int> &oldAssignments,
    std::map <unsigned int, unsigned int> &AMcur,
    std::map <unsigned int, unsigned int> &AMmin,
    std::map<unsigned int, cv::Mat> &candidateCovarMat,
    std::map<unsigned int, cv::Mat> &candidateMeanMat,
    std::map<unsigned int, cv::Mat> &newMeanMat,
    double &minSUM,
    cvb::CvBlobs &NEWBlobs
    )
{
  std::stringstream DEBUG;
  if(NEW.size()==1)
  {
    std::map <unsigned int, unsigned int>::iterator c=oldAssignments.begin();
    while(c!=oldAssignments.end())
    {
      if(c->second==0)
      {
        AMcur[NEW[0]]=c->first;
        c->second=NEW[0];
        
        double curSUM = calculate_MH_sum(AMcur, 
            candidateCovarMat,
            candidateMeanMat,
            newMeanMat,
            NEWBlobs
            );

        if (curSUM<minSUM)
        {
          DEBUG << "Minimum SUM found: " << curSUM << " with assignment " <<
            printUIMap(AMcur);
          verbosePrint(DEBUG);
          minSUM=curSUM;
          AMmin=AMcur;
        }

        c->second=0;
      }
      ++c;
    }
  }

  if(NEW.size()>1)
  {
    std::map <unsigned int, unsigned int>::iterator c=oldAssignments.begin();
    while(c!=oldAssignments.end())
    {
      if(c->second==0)
      {
        AMcur[NEW[0]]=c->first;
        c->second=NEW[0];
        std::vector<unsigned int> sub(NEW.begin()+1,NEW.end());

        assign_combinations(
            sub,
            oldAssignments,
            AMcur,
            AMmin,
            candidateCovarMat,
            candidateMeanMat,
            newMeanMat,
            minSUM,
            NEWBlobs
            );
        
        c->second=0;
      }
      ++c;

    }
  }

}

int diverge_match_short(
  std::vector<unsigned int> &candidateLarvae,
  std::vector<unsigned int> &newLarvae,
  std::map<unsigned int, unsigned int> &newAssignments,
  cvb::CvBlobs &NEW)
{
  std::vector<unsigned int>::iterator NewIT=newLarvae.begin();
  int matched=0;
  while(NewIT!=newLarvae.end())
  {
    std::vector<unsigned int>::iterator PreIT=candidateLarvae.begin();
    while(PreIT!=candidateLarvae.end())
    {
      unsigned int lastIdx=detected_larvae[*PreIT].lastBlobWithStats;
      cvb::CvBlob &blobP=detected_larvae[*PreIT].blobs[lastIdx];

      if(centresMatch(&blobP,NEW[*NewIT],0.20) && 
          blobSizeIsRelevant(&blobP,NEW[*NewIT]))
      {
        newAssignments[*NewIT]=*PreIT;
        matched++;
      }
      ++PreIT;
    }
    ++NewIT;
  }
  return matched;
}

void diverge_match_new(
  std::vector<unsigned int> &candidateLarvae,
  std::vector<unsigned int> &newLarvae,
  std::map<unsigned int, unsigned int> &newAssignments,
  double duration,
  cvb::CvBlobs &NEW)
{

  if(duration<0.51)
  {
    int matched=
      diverge_match_short(candidateLarvae,newLarvae,newAssignments,NEW);
    if(matched>0)
    {
      std::cerr << printUIMap(newAssignments) << std::endl;
      return;
    }
  }

  std::stringstream DEBUG;
  std::map<unsigned int, cv::Mat> candidateCovarMat;
  std::map<unsigned int, cv::Mat> candidateMeanMat;

  std::map<unsigned int, cv::Mat> newMeanMat;

  std::vector<unsigned int>::iterator cIT=candidateLarvae.begin();

  verbosePrint("Setting up candidate larvae data");
  while(cIT!=candidateLarvae.end())
  {
    /*
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
*/
    double size_avg=avgNVec(detected_larvae[*cIT].area);
    double grey_value_avg=avgNVec(detected_larvae[*cIT].grey_value);
    double length_avg=avgNVec(detected_larvae[*cIT].length);
    double perimeter_avg=avgNVec(detected_larvae[*cIT].perimeter);
    double width_avg=avgNVec(detected_larvae[*cIT].width);
    double speed_x=avgNVec(detected_larvae[*cIT].midpoint_speed_x);
    double speed_y=avgNVec(detected_larvae[*cIT].midpoint_speed_y);


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

    cv::hconcat(InputArray,
        cv::Mat(detected_larvae[*cIT].midpoint_speed_x),
        InputArray);

    cv::hconcat(InputArray,
        cv::Mat(detected_larvae[*cIT].midpoint_speed_y),
        InputArray);

    std::vector<double> meanVec;
    meanVec.push_back(size_avg);
    meanVec.push_back(grey_value_avg);
    meanVec.push_back(length_avg);
    meanVec.push_back(perimeter_avg);
    meanVec.push_back(width_avg);
    meanVec.push_back(speed_x);
    meanVec.push_back(speed_y);
    
    cv::Mat meanMat(meanVec);
    cv::Mat meanTMat;
    cv::transpose(meanMat,meanTMat);

    cv::Mat covarMat;

    cv::calcCovarMatrix(InputArray, covarMat,meanTMat,CV_COVAR_ROWS|CV_COVAR_NORMAL|CV_COVAR_USE_AVG);
    cv::invert(covarMat,covarMat,cv::DECOMP_SVD);
    covarMat.copyTo(candidateCovarMat[*cIT]);
    meanMat.copyTo(candidateMeanMat[*cIT]);

    ++cIT;
  }

  verbosePrint("Setting up new larvae data");
  cIT=newLarvae.begin();
  while(cIT!=newLarvae.end())
  {
    //Setup of each new larva
    cv::Mat larvaROI;
    std::vector<cv::Point2f> newLarvaPoints;
    blobToPointVector(*NEW[*cIT],newLarvaPoints);
    createLarvaContour(larvaROI,(*NEW[*cIT]));
    larvaDistanceMap dstLarva(newLarvaPoints);
    cv::Point2f centroid;
    centroid.x=NEW[*cIT]->centroid.x;
    centroid.y=NEW[*cIT]->centroid.y;
    computeSpine(*NEW[*cIT],dstLarva);

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
    meanMat.copyTo(newMeanMat[*cIT]);
    ++cIT;
  }

  double minSUM=999999999.0;
  std::map <unsigned int, unsigned int> AMmin;
  std::map <unsigned int, unsigned int> AMcur;
  std::map <unsigned int, unsigned int> oldAssignments;

  std::vector<unsigned int> newLarvaeVec; //vector of objects that we are
                                          // more certain are larvae
  for(unsigned int i=0 ; i<candidateLarvae.size() ; i++)
  {
    oldAssignments[candidateLarvae[i]]=0;
  }

  cv::Mat samples;

  for(unsigned int i=0; i<newLarvae.size();++i)
  {
    if(is_larva(NEW[newLarvae[i]])<IS_LARVA_THRESHOLD)
    {
      newLarvaeVec.push_back(newLarvae[i]);
    }
  }

  
  assign_combinations(
      newLarvaeVec,
      oldAssignments,
      AMcur,
      AMmin,
      candidateCovarMat,
      candidateMeanMat,
      newMeanMat,
      minSUM,
      NEW
      );

  //if(newLarvaeVec.size()>1 || minSUM/newLarvaeVec.size()<LARVA_MAHALANOBIS_THRESHOLD)
  if( (AMmin.size()==newLarvaeVec.size() && AMmin.size() > 1 &&
        minSUM/AMmin.size()<3*LARVA_MAHALANOBIS_THRESHOLD
       ) 
      || minSUM/newLarvaeVec.size()<LARVA_MAHALANOBIS_THRESHOLD)
  {
    newAssignments=AMmin;
    std::cerr << printUIMap(AMmin) << std::endl;
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
  cv::Point2f mdpLarva1, mdpLarva2;
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
  std::vector<cv::Point2f> newLarva1Points;
  std::vector<cv::Point2f> newLarva2Points;
  blobToPointVector(*newLarva1,newLarva1Points);
  blobToPointVector(*newLarva2,newLarva2Points);
  larvaDistanceMap dstLarva1(newLarva1Points), dstLarva2(newLarva2Points);

  cv::Point2f centroid1;
  centroid1.x=newLarva1->centroid.x;
  centroid1.y=newLarva1->centroid.y;
  cv::Point2f centroid2;
  centroid2.x=newLarva2->centroid.x;
  centroid2.y=newLarva2->centroid.y;

  //larvaSkel newLarvaSkel1(larvaROI1,centroid1);
  //larvaSkel newLarvaSkel2(larvaROI2,centroid2);

  //computeInnerDistances(*newLarva1,dstLarva1,newLarvaSkel1.MidPoint);
  //computeInnerDistances(*newLarva2,dstLarva2,newLarvaSkel2.MidPoint);
  computeSpine(*newLarva1,dstLarva1);
  computeSpine(*newLarva2,dstLarva2);


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

  cv::calcCovarMatrix(InputArrayA, covarMatA,meanTMatA,CV_COVAR_ROWS|CV_COVAR_NORMAL|CV_COVAR_USE_AVG);
  cv::calcCovarMatrix(InputArrayB, covarMatB,meanTMatB,CV_COVAR_ROWS|CV_COVAR_NORMAL|CV_COVAR_USE_AVG);

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
		            std::vector<unsigned int> &nearbyLarvae,bool pre=true,
                double PADRatio=2)
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
        cBlob->centroid.y > (Blob->miny - MaxDist/2) &&
        ( pre==false && assignedNew[It->first].size()<=0 || 
          pre==true  && assignedPrevious[It->first].size()<=0)
        )
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
  //unsigned int MINID=0;
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

void findInDistance(cvb::CvBlobs &In, cvb::CvBlob *BLOB, double dist, std::map <unsigned int, double> neighbours)
{
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

void assign_one(unsigned int preID,unsigned int postID)
{
  std::stringstream DEBUG;
  assignedPrevious[preID].push_back(postID);
  assignedNew[postID].push_back(preID);
  DEBUG << "Assigning: " << postID << " -> " << preID;
  verbosePrint(DEBUG);
  assignedPreMap[preID]=1;
}

void assign_one(unsigned int preID,
                std::vector<unsigned int> postID)
{
  std::stringstream DEBUG;
  std::vector<unsigned int>::iterator postIT=postID.begin();
  assignedPreMap[preID]=postID.size();
  while(postIT!=postID.end())
  {
    unsigned int NewID=++LARVAE_COUNT;
    assignedNew[*postIT].push_back(NewID);
    DEBUG << "Assigning: " << *postIT << " -> " << NewID;
    verbosePrint(DEBUG);
    assignedPrevious[preID].push_back(NewID);
    ++postIT;
  }
  DEBUG << "Cluster " << preID << " diverged into new larvae: " << printVector(assignedPrevious[preID]);
  verbosePrint(DEBUG);
  assignedPreMap[preID]=postID.size();
}

void assign_one(std::vector<unsigned int> preID,
                unsigned int postID,unsigned int newID)
{
  std::stringstream DEBUG;
  std::vector<unsigned int>::iterator preIT=preID.begin();
  assignedNew[postID].push_back(newID);
  DEBUG << "New cluster " << postID << " from larvae: " << printVector(preID);
  verbosePrint(DEBUG);
  while(preIT!=preID.end())
  {
    assignedNew[postID].push_back(*preIT);
    assignedPrevious[*preIT].push_back(postID);
    assignedPreMap[*preIT]=1;
    ++preIT;
  }
}

void assign_diverging(cvb::CvBlobs &New,
                      unsigned int CLUSTER_ID,
                      std::vector<unsigned int> &IDs
                      )
{
  std::stringstream DEBUG;
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
    DEBUG << "Cluster " << CLUSTER_ID << " is new. Assigning new IDs for diverged larvae";
    verbosePrint(DEBUG);
    assign_one(CLUSTER_ID,IDs);
    detected_larvae[CLUSTER_ID].isCluster=true;
    //TODO: IDs is the labels as returned by CvBlob we need to fix this!!
  }
  else
  {
    std::vector<unsigned int> candidateLarvae(dcIT->second.begin()+1,
        dcIT->second.end());
    // New assignments are returned by the diverge_match_new function
    // [ NEW_ID -> OLD_ID ]
    std::map<unsigned int, unsigned int> newAssignments;
    /*
    std::cout << "Callgin diverge match: CF: " << CURRENT_FRAME <<
      " BLOB SF: " << detected_larvae[dcIT->first].start_frame << std::endl;
      */
    diverge_match_new(candidateLarvae,
        IDs,
        newAssignments,
        (CURRENT_FRAME - detected_larvae[dcIT->first].start_frame)/24.0,
        New);

    DEBUG << "New assignments return by diverge matching: " << std::endl << printUIMap(newAssignments);
    verbosePrint(DEBUG);
    // Perform the assignments

    //newCluster to register the assignments and what is left in the cluster
    std::vector<unsigned int> 
      newCluster(dcIT->second.begin()+1,dcIT->second.end()); 

    //Initiate a vector the the IDs to asign to keep track of what is assigned.
    std::vector<unsigned int> newIDs=IDs;

    for (std::map<unsigned int, unsigned int>::iterator naIT=
        newAssignments.begin();
        naIT!=newAssignments.end();
        ++naIT)
    {
      if(naIT->second!=0)
      {
        assign_one(naIT->second,naIT->first);
        assignedPreMap[CLUSTER_ID]=assignedPreMap[CLUSTER_ID]+1;
        assignedPrevious[CLUSTER_ID].push_back(naIT->second);

        for(std::vector<unsigned int>::iterator erIT=newCluster.begin();
            erIT!=newCluster.end();++erIT)
        {
          if(*erIT==naIT->second)
          {
            newCluster.erase(erIT);
            break;
          }
        }
        for(std::vector<unsigned int>::iterator erIT=newIDs.begin();
            erIT!=newIDs.end();++erIT)
        {
          if(*erIT==naIT->first)
          {
            newIDs.erase(erIT);
            break;
          }
        }
      }
    }

    // We have performed the assignents and what is left is contained in
    // the newIDs and the newCluster vectors
    if(newIDs.size()==newCluster.size() && newCluster.size()==0)
    {
      //perfect assignment. return from function
      detected_clusters.erase(dcIT);
      return;
    }
    if(newIDs.size()==newCluster.size() && newCluster.size()==1)
    {
      // Object diverged and left only one unassigned and 
      // only one is left in the cluster
      // We force the assignment
      assign_one(newCluster[0],newIDs[0]);
      assignedPreMap[CLUSTER_ID]=assignedPreMap[CLUSTER_ID]+1;
      assignedPrevious[CLUSTER_ID].push_back(newCluster[0]);
      detected_clusters.erase(dcIT);
      return;
    }
    if(newIDs.size()==1 && newCluster.size()>1)
    {
      //only one object remains to be assigned but our cluster was bigger
      unsigned int CLUSTER_ID=++LARVAE_COUNT;

      detected_clusters[CLUSTER_ID].push_back(newCluster.size());
      detected_clusters[CLUSTER_ID].insert(
          detected_clusters[CLUSTER_ID].end(),
          newCluster.begin(),
          newCluster.end());

      assign_one(newCluster,newIDs[0],CLUSTER_ID);
      detected_clusters.erase(dcIT);
    }
    //TODO: Figure out the case for newIDs.size()>1
    //  a) give new numbers to those that look like larvae
    //  b) create new clusters with dummy? larvae (if necessary i.e. if )
    //     plus those that were there
    //
    // IMPORTANT!! LIMIT THE SIZE
    if(newIDs.size()==2 && newCluster.size()==2)
    {
      if(is_larva(New[newIDs[0]])>1000 && is_larva(New[newIDs[1]])>1000)
      {
        //unsigned int CLUSTER_ID=++LARVAE_COUNT;
        //detected_clusters[CLUSTER_ID]=newCluster;
        //assign_one(newCluster,newIDs[0],CLUSTER_ID);
        detected_clusters.erase(dcIT);
      }
      if(is_larva(New[newIDs[0]])<1000 && is_larva(New[newIDs[1]])>1000)
      {
        // Unhandled
      }
    }
    // Unhandled should produce an message so that we debug it
    std::cerr << "Unhandled case were more than 1 IDs were unassigned"
      << std::endl;
  }
}

void assign_clustering(
                      unsigned int POST_ID,
                      std::vector<unsigned int> &IDs
                      )
{
  unsigned int CLUSTER_ID=++LARVAE_COUNT;
  newClusters.push_back(CLUSTER_ID);
  std::vector<unsigned int> contents;
  contents.push_back(IDs.size());
  detected_clusters[CLUSTER_ID]=contents;
  std::map<unsigned int,std::vector<unsigned int> >::iterator dcIT;
  assign_one(IDs,POST_ID,CLUSTER_ID);
  std::vector<unsigned int>::iterator IT=IDs.begin();
  while (IT!=IDs.end())
  {
    if((dcIT=detected_clusters.find(*IT))!=detected_clusters.end())
    {
      detected_clusters[CLUSTER_ID].insert(
          detected_clusters[CLUSTER_ID].end(),
          dcIT->second.begin()+1,
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
    }
    ++IT;
  }
}

bool centresMatch(
    cvb::CvBlobs &In, 
    cvb::CvBlob *blob,
    std::vector<unsigned int> &larvae, 
    double factor=LARVA_CENTRE_COMPARISON_FACTOR-1)
{
  std::stringstream DEBUG;
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
  DEBUG << "CentresMatch [" << blob->label << ", " << printVector(larvae) << "]: Length: " << objectLength << " Difx: " << fabs(blob->centroid.x - xcomb) << " Dify: " << fabs(blob->centroid.y - ycomb) << " Threshold: " << factor*objectLength;
  verbosePrint(DEBUG);
  if (fabs(blob->centroid.x - xcomb)< factor*objectLength &&
      fabs(blob->centroid.y - ycomb)< factor*objectLength )
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
  verbosePrint("Trying to detect diverging clusters");

  std::stringstream DEBUG;
  DEBUG<< "Size of newLarvaeNearby: "  << newLarvaeNearby.size();
  verbosePrint(DEBUG);

  if(newLarvaeNearby.size()<=1)
    return -1; // No diverging is possible
  // Check the complete set first
  std::vector<unsigned int>::iterator pIT=preLarvaeNearby.begin();
  while(pIT!=preLarvaeNearby.end())
  {
    DEBUG << "Checking if nodes " << printVector(newLarvaeNearby) << " diverged from: " << *pIT;
    verbosePrint(DEBUG);
    if(centresMatch(New,Pre[*pIT],newLarvaeNearby))
    {
      // Centres of all candidates match with new blob
      // cluster contains all. We can return
      DEBUG << "Node " << *pIT << " matches " << printVector(newLarvaeNearby);
      verbosePrint(DEBUG);
      assign_diverging(New,*pIT,newLarvaeNearby);
      break;
    }
    if(newLarvaeNearby.size()>2)
    {
      DEBUG << "Checking powersets of " << printVector(newLarvaeNearby) << " that diverged from: " << *pIT;
      verbosePrint(DEBUG);
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
    }
    ++pIT;
  }
  return 0;
}

int detect_clustering(std::vector<unsigned int> &preLarvaeNearby,
                       std::vector<unsigned int> &newLarvaeNearby,
                       cvb::CvBlobs &Pre,
                       cvb::CvBlobs &New)
{
  verbosePrint("Trying to detect clusters");

  std::stringstream DEBUG;
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
      DEBUG << "Centres of new larva: " << *nIT << " match centres of " 
       << printVector(preLarvaeNearby) << ". Assigning clustering";
      verbosePrint(DEBUG);

      assign_clustering(*nIT,preLarvaeNearby);
      //continue;
      break;
    }
    std::vector<std::vector<unsigned int> > pSETS;
    powersets(preLarvaeNearby,pSETS);
    std::vector<std::vector<unsigned int> >::iterator pSIT=pSETS.begin();
    while(pSIT!=pSETS.end())
    {
      verbosePrint("Trying to detect subsets for clustering");
      if(centresMatch(Pre,New[*nIT],*pSIT))
      {
        // Centres of all candidates match with new blob
        // cluster contains subset pSIT. We can return
        DEBUG << "Centres of new larva: " << *nIT << " match centres of subset " 
          << printVector(*pSIT) << ". Assigning clustering";
        verbosePrint(DEBUG);
        assign_clustering(*nIT,*pSIT);
        break;
      }
      ++pSIT;
    }
    ++nIT;
  }
  return 0;
}

bool checkLarvae()
{
  std::map<unsigned int,larvaObject>::iterator lit=detected_larvae.begin();
  while(lit!=detected_larvae.end())
  {
    if (lit->first != lit->second.larva_ID)
    {
      std::cerr << "INCONSISTENCY!!! [ID] doesn't match larva_ID" << std::endl;
      return false;
    }
    if (lit->first != lit->second.blobs.back().label)
    {
      std::cerr << "INCONSISTENCY!!! [ID: " << lit->first << "] doesn't match blob label: " << lit->second.blobs.back().label << " !!!" << std::endl;
      return false;
    }
    if (lit->second.larva_ID != lit->second.blobs.back().label)
    {
      std::cerr << "INCONSISTENCY!!! larva_ID doesn't match blob label" << std::endl;
      return false;
    }
    ++lit;
  }
  return true;
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
  newInFrame.clear();
  newClusters.clear();
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
      getNearbyLarvae(Prev,preBlob,preLarvaeNearby,true);
      getNearbyLarvae(In,preBlob,postLarvaeNearby,false);

      DEBUG << "Around larva " << prevIt->first << " we found " << preLarvaeNearby.size() << " larvae from the previous frame and " << postLarvaeNearby.size() << " from the new frame.";
      verbosePrint(DEBUG);

      DEBUG << "Around larva " << prevIt->first << " we found pre: " << printVector(preLarvaeNearby) << " and post: " << printVector(postLarvaeNearby);
      verbosePrint(DEBUG);
      // Great case collection now:

      // CASE 1 -> 1
      if(preLarvaeNearby.size()==1 && postLarvaeNearby.size()==1)
      {
        if(blobSizeIsRelevant(In[postLarvaeNearby[0]],preBlob))
        {
          if(centresMatch(In,Prev[preLarvaeNearby[0]],postLarvaeNearby,0.3))
          {
            DEBUG << "Larva " << prevIt->first << " from previous frame matches " << postLarvaeNearby[0] << " from new frame. Assigning.";
            verbosePrint(DEBUG);
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
          bool erased=false;
          std::vector<unsigned int>::iterator postNearIt=postLarvaeNearby.begin();
          while(postNearIt!=postLarvaeNearby.end())
          {
            if (blobSizeIsRelevant(Prev[*preNearIt],In[*postNearIt]) &&
                centresMatch(Prev[*preNearIt],In[*postNearIt],0.2))
            {
              //1-1 and one extra in sight
              //New Larva matches preID larva therefore assign and ignore the
              //other.
              DEBUG << "Larva " << *preNearIt << " from previous frame matches " << *postNearIt << " from new frame. Assigning";
              verbosePrint(DEBUG);
              assign_one(*preNearIt,*postNearIt);
              // erase the obvious assignments
              try{
                preLarvaeNearby.erase(preNearIt);
                postLarvaeNearby.erase(postNearIt);
              }
              catch(...)
              {
                std::cerr << "Problem erasing: " << *preNearIt << " or ";
                std::cerr << *postNearIt << std::endl;
              }
              erased=true;
            }
            else{
              DEBUG << "Larva " << *preNearIt << " from previous frame doesn't match " << *postNearIt << " from new frame.";
              verbosePrint(DEBUG);
              erased=false;
              ++postNearIt;
            }
          }
          if(!erased)
          {
            ++preNearIt;
            break;
          }
        }
        // The rest are either appearing/disappearing/clustering/diverging
        // BUT if the main one was indeed assigned then we are in luck 
        // and we move on
        DEBUG << "After initial assignment assignedPreMap: " << printUIMap(assignedPreMap) << " for larva with Id: " << preID;
        verbosePrint(DEBUG);
        if(assignedPreMap[preID]>0)
        {
          continue;
        }

        detect_clustering(preLarvaeNearby,postLarvaeNearby,Prev, In);

        DEBUG << "After clustering check assignedPreMap: " << printUIMap(assignedPreMap) << " for larva with Id: " << preID;
        verbosePrint(DEBUG);
        if(assignedPreMap[preID]>0)
        {
          continue;
        }

        detect_diverging(preLarvaeNearby,postLarvaeNearby,Prev, In);
        DEBUG << "After diverging check assignedPreMap: " << printUIMap(assignedPreMap) << " for larva with Id: " << preID;
        verbosePrint(DEBUG);
        if(assignedPreMap[preID]>0)
        {
          continue;
        }
        else
          verbosePrint("FOUND TESTCASE FOR MIXED DIVERGING/CONVERGING");
      }
    }
    ++prevIt;
    verbosePrint("assignedPrevious: ");
    std::map<unsigned int, std::vector<unsigned int> >::iterator prI=assignedPrevious.begin();
    while(prI!=assignedPrevious.end())
    {
      if(prI->second.size()>0)
      {
        DEBUG << prI->first << " -> " << printVector(prI->second);
        verbosePrint(DEBUG);
      }
      else
      {
        DEBUG << prI->first << " -> []" ;
        verbosePrint(DEBUG);
      }
      prI++;
    }

    verbosePrint("assignedNew: ");
    prI=assignedNew.begin();
    while(prI!=assignedNew.end())
    {
      if(prI->second.size()>0)
      {
        DEBUG << prI->first << " -> " << printVector(prI->second);
        verbosePrint(DEBUG);
      }
      else
      {
        DEBUG << prI->first << " -> []" ;
        verbosePrint(DEBUG);
      }
      prI++;
    }
  }


  cvb::CvBlobs::iterator bIT=In.begin();
  while(bIT!=In.end())
  {
    if (assignedNew[bIT->first].size()>0)
    {
      bIT->second->label=assignedNew[bIT->first][0];
      out[assignedNew[bIT->first][0]]=bIT->second;
    }
    else
    {
      unsigned int NEWID=++LARVAE_COUNT;
      verbosePrint("Extra Assignment:");
      DEBUG << bIT->first << " -> " << NEWID;
      newInFrame.push_back(NEWID);
      verbosePrint(DEBUG);
      bIT->second->label=NEWID;
      out[NEWID]=bIT->second;
    }
    bIT++;
  }
  updateLarvae(out,Prev);
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
      //double XVAL,YVAL,cur=0;
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
      )

      ("batch-testing,B",
       po::value<bool> (&LRVTRACK_SHOW_HEAD_TAIL),
       "Show no output except from the matchings. Used for testing internal parameters."
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

  return 0;
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

  return 0;
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
            continue;

          ++larvaeToConsider;
          lifespanSUM+=(CURRENT_FRAME-cl.start_frame)/VIDEO_FPS;
          // avg speed
          double xvel,yvel;
          if(cl.midpoint_speed_x.size()>0)
            xvel=cl.midpoint_speed_x.back();
          else
            xvel=0.0;
          if(cl.midpoint_speed_y.size()>0)
            yvel=cl.midpoint_speed_y.back();
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
                   cl.lrvDistances.back().p20,
                   cl.tails.back());
          a2=angle(cl.heads.back(),
                   cl.lrvDistances.back().p80,
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
      
      bool dpc=false;
      if(newInFrame.size()>0)
        {
          if (!dpc)
          {
            summary << " %% ";
            dpc=true;
          }
          for (unsigned int i=0; i<newInFrame.size(); ++i)
            {
              summary << "0 " << newInFrame[i] << " ";
            }
        }
      if(newClusters.size()>0)
      {
        if (!dpc)
        {
          summary << " %% ";
          dpc=true;
        }
        for (unsigned int i=0; i<newClusters.size(); ++i)
        {
          std::vector<unsigned int> cluster_nodes=
            detected_clusters[newClusters[i]];
          for(unsigned int j=1; j<cluster_nodes.size() ; ++j)
          {
            summary << cluster_nodes[j] << " " << newClusters[i] << " ";
          }
        }
      }
      if(assignedPreMap.size()>0)
      {
        std::map<unsigned int,unsigned int>::iterator It=assignedPreMap.begin();
        while(It!=assignedPreMap.end())
        {
          if(It->second==0)
          {
            if (!dpc)
            {
              summary << " %% ";
              dpc=true;
            }
            summary << It->first << " 0 " ;
          }
          else
          {
            std::vector<unsigned int> &A = assignedPrevious[It->first];
            if (A.size()>1)
            {
              if (!dpc)
              {
                summary << " %% ";
                dpc=true;
              }
              for(unsigned int i=0; i<A.size();++i)
              {
                summary << It->first << " " << A[i] << " ";
              }
            }
          }
          ++It;
        }
      }
      summary << std::endl;
    }
}

void printBlobFile(larvaObject &lrv)
{
  
  if(lrv.isCluster)
    return;
  std::ostringstream BLOBFILENAME;
  std::ofstream blobFile;

  //double elapsed=(tC.tv_sec - tS.tv_sec) + ((tC.tv_usec - tS.tv_usec)/1000000.0);
  cpu_times elapsed(tS.elapsed());
  BLOBFILENAME <<
               LRVTRACK_NAME <<
               "_" <<
               std::fixed <<
               std::setfill('0') <<
               std::setw(5) <<
               lrv.larva_ID <<
               ".blob";

  fs::path experimentFolderPath(LRVTRACK_RESULTS_FOLDER);
  experimentFolderPath = experimentFolderPath / LRVTRACK_DATE;

  fs::path blobFilePath(
    experimentFolderPath / BLOBFILENAME.str()
  );

  blobFile.open(blobFilePath.c_str(),std::ofstream::out | std::ofstream::app);

  blobFile << CURRENT_FRAME-START_FRAME+1 << " " ;
  blobFile << std::left
               << std::fixed
               << std::setfill('0')
               << std::setprecision(3)
               << (double) elapsed.wall/1000000000.0
               << "  ";

  blobFile<< std::left
               << std::fixed
               << std::setfill('0')
               << std::setprecision(3)
               << lrv.blobs.back().centroid.x
               << " ";

  blobFile << std::left
               << std::fixed
               << std::setfill('0')
               << std::setprecision(3)
               << lrv.blobs.back().centroid.y
               << "  ";

  blobFile << lrv.blobs.back().area << "  ";

  cv::Point2f major,minor;
  double long_axis,short_axis;

  principalAxes(lrv.blobs.back(),major,minor,long_axis,short_axis);

  blobFile << std::left
               << std::fixed
               << std::setfill('0')
               << std::setprecision(3)
               << major.x
               << " ";

  blobFile << std::left
               << std::fixed
               << std::setfill('0')
               << std::setprecision(3)
               << major.y
               << "  ";

  blobFile << std::left
               << std::fixed
               << std::setfill('0')
               << std::setprecision(3)
               << sqrt(minor.y*minor.y+minor.x*minor.x)
               << "  ";

  blobFile << std::left
               << std::fixed
               << std::setfill('0')
               << std::setprecision(1)
               << long_axis
               << " ";

  blobFile << std::left
               << std::fixed
               << std::setfill('0')
               << std::setprecision(1)
               << short_axis
               << " ";

  /*
  blobFile << "% ";
   
  larvaDistanceMap &dm=lrv.lrvDistances.back();
  cv::Point2f c;
  c.x=lrv.blobs.back().centroid.x;
  c.y=lrv.blobs.back().centroid.y;
  if(dm.spineSegments != dm.spinePoints.size())
  {
    std::cerr << "PROBLEM COMPUTING SPINE" << std::endl;
  }
  for (unsigned int i=0; i<dm.spineSegments ;i++)
  {
    blobFile << (int) (dm.spinePoints[i].x-c.x) << " " ;
    blobFile << (int) (dm.spinePoints[i].y-c.y) << " " ;
  }
*/
  blobFile << "%% ";

  std::string cntstr;
  cv::Point first;
  unsigned int cntSize;
  createLarvaContourPacked(first,cntSize,cntstr,lrv.blobs.back());

  blobFile << first.x << " ";
  blobFile << first.y << " ";
  blobFile << cntSize << " ";
  blobFile << cntstr;
  blobFile << std::endl;
}

void colorInvBW(cv::Mat &A)
{
  for(int i=0;i<A.rows;i++)
    for(int j=0;j<A.cols;j++)
      for(int k=0;k < A.channels() ;k++)  //loop to read for each channel
        A.data[i*A.step+j*A.channels()+k]=255-A.data[i*A.step+j*A.channels()+k];    //inverting the image

}

int main(int argc, char* argv[])
{
  bool SHOWTAGS=true;
  bool TRACK=false;
  bool STEP=true;

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
          //std::cout << VIDEO_FPS << std::endl;
        }
      else
        {
          VIDEO_FPS=capture.get(CV_CAP_PROP_FPS);
          //std::cout << VIDEO_FPS << std::endl;
        }
    }

  cv::Mat grey_bgFrame;
  //unsigned int fg_threshold=10;
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
  int votes=180;
  while (circles.size()==0 && votes >= 100)
  {
    cv::HoughCircles(grey_bgFrame, circles, CV_HOUGH_GRADIENT,
        1,   // accumulator resolution (size of the image / 2)
        850,  // minimum distance between two circles
        50, // Canny high threshold
        votes, // minimum number of votes
        bgFrame.rows/3, bgFrame.rows/1.95); // min and max radiusV

    votes-=10;

  }
  std::vector<cv::Vec3f> cups;
  if (LRVTRACK_ODOUR_CUPS>0 || LRVTRACK_ODOUR_CUPS==-1)
    {
      cv::Mat fgROI;
      fgROI=cv::Mat::zeros(grey_bgFrame.rows , grey_bgFrame.cols,grey_bgFrame.depth());
      if(circles.size()>0)
        {
          cv::circle(fgROI, cv::Point2f(circles[0][0],circles[0][1]),int(circles[0][2]/1.1),cv::Scalar(255),-1);
        }

      cv::Mat thr;

      thr=grey_bgFrame&fgROI;

      cv::morphologyEx(thr,thr,cv::MORPH_OPEN,cv::Mat(),cv::Point2f(-1,-1),5);
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
                     cv::Point2f((*itc)[0], (*itc)[1]), // circle centre
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

  /*cv::Ptr<cv::superres::SuperResolution> srobj =cv::superres::createSuperResolution_BTVL1_GPU();
  cv::Ptr<cv::superres::DenseOpticalFlowExt> opflow = cv::superres::createOptFlow_Farneback_GPU();
  srobj->set("scale", 2);
  srobj->set("opticalFlow", opflow);

  cv::Ptr<cv::superres::FrameSource> frameSource = 
    cv::superres::createFrameSource_Video_GPU(LRVTRACK_FILE_INPUT);
  srobj->setInput(frameSource);
  cv::Mat frame_super;
*/

  cvb::CvBlobs preBlobs;
  while (!stop)
    {
      // read next frame if any
      if (!capture.read(frame))
        break;
      //srobj->nextFrame(frame);
      std::stringstream F;
      F<< CURRENT_FRAME;
      cv::putText(frame,
          F.str(),
          cv::Point2f(20,40),
          cv::FONT_HERSHEY_PLAIN,
          0.8,
          cv::Scalar(255,255,255),
          1,
          CV_AA);

      //cv::Mat image;
      //cv::GaussianBlur(frame, image, cv::Size(0, 0), 3);
      //cv::addWeighted(frame, 1.7, image, -0.7, 0, image);
      
      if(TRACK)
        {
          if(START_FRAME==0)
            {
              std::cout << "FRAME,ID,ANGLE,ORIENTATION,SPEEDX,SPEEDY,CENTROIDX,CENTROIDY,INBLOB,ISBLOB" << std::endl;
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
              cv::circle(fgROI, cv::Point2f(circles[0][0],circles[0][1]),int(circles[0][2]/0.99),cv::Scalar(255),-1);
              cv::circle(frame, cv::Point2f(circles[0][0],circles[0][1]),int(circles[0][2]/0.99),cv::Scalar(0,255,0),1);
            }

          if(cups.size()>0)
            {
              for(unsigned int i=0; i<cups.size(); ++i)
                {
                  cv::circle(fgROI, cv::Point2f(cups[i][0],cups[i][1]),int(cups[i][2]),cv::Scalar(0),-1);
                  cv::circle(frame, cv::Point2f(cups[i][0],cups[i][1]),int(cups[i][2]),cv::Scalar(0,0,255),1);
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

          //unsigned int result=
          cvLabel(&ipl_thresholded, labelImg, blobs);
          cvb::cvFilterByArea(blobs, 36, 900);
          
          //----------------------------------------------------------
        /*
          cv::Mat debugImg;
          frame.copyTo(debugImg);
          cvb::CvBlobs::iterator dbg=blobs.begin();
          dbg=blobs.begin();
          while (dbg!=blobs.end())
          {
            std::stringstream sstm;
            cvb::CvBlob *blob=(*dbg).second;
            sstm << (*dbg).first;
            cv::putText(debugImg,
                sstm.str(),
                cv::Point2f2d(blob->centroid.x+12,blob->centroid.y+12),
                cv::FONT_HERSHEY_PLAIN,
                0.8,
                cv::Scalar(255,255,255),
                1,
                CV_AA);
            ++dbg;
          }
          cv::imshow("DEBUG",debugImg);
          cv::waitKey(1);
          */
          //------------------------------------------------------------


          cvb::CvBlobs tracked_blobs;
          cvb::CvBlobs blob1;
          cvb::CvBlobs::iterator it=blobs.begin();

          if(preBlobs.size()>0)
            {
              FrameEllapsedTime = tP.elapsed();
              CurrentTime= tS.elapsed();
              //larvae_track(blobs,preBlobs,tracked_blobs);
              newLarvaeTrack(blobs,preBlobs,tracked_blobs);
              tP.start();

              printSummary(preBlobs,blobs,false);

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
                  cvb::CvBlob *blob=it->second;
                  printBlobFile(detected_larvae[it->first]);
                  try
                    {
                      std::vector<unsigned int> &clusterMBs=detected_clusters.at(it->first);
                      sstm << it->first << printVector(clusterMBs,1);
                      cv::putText(frame,
                                  sstm.str(),
                                  cv::Point2f(blob->centroid.x+12,blob->centroid.y+12),
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
                                  cv::Point2f(blob->centroid.x+12,blob->centroid.y+12),
                                  cv::FONT_HERSHEY_PLAIN,
                                  0.8,
                                  cv::Scalar(255,255,255),
                                  1,
                                  CV_AA);
                    }

                  // Here we define an extra PADDING
                  int PAD=2;
                  cv::Mat larvaROI;
                  try{
                    //(0 <= roi.x && 0 <= roi.width && roi.x + roi.width <= m.cols && 0 <= roi.y && 0 <= roi.height && roi.y + roi.height <= m.rows)
                  larvaROI=cv::Mat(frame,
                                   cv::Rect(blob->minx-ROI_PADDING-PAD,
                                            blob->miny-ROI_PADDING-PAD,
                                            blob->maxx-blob->minx+2*(ROI_PADDING+PAD),
                                            blob->maxy-blob->miny+2*(ROI_PADDING+PAD))
                                  );
                  }
                  catch(...)
                  {
                    std::cerr << "larvaROI failed. continuing" << std::endl;
                    break;
                  }
                  cv::circle(frame,
                             cv::Point2f(blob->centroid.x,blob->centroid.y),
                             1,
                             cv::Scalar(255,0,0),
                             -1);

                  if(detected_larvae[it->first].isCluster==false)
                    {
                      /*
                      cv::circle(larvaROI,
                                 cv::Point2f(
                                   detected_larvae[it->first].lrvskels.back().Point20.x+PAD,
                                   detected_larvae[it->first].lrvskels.back().Point20.y+PAD),
                                 1,
                                 cv::Scalar(255,255,0),
                                 -1);
                      */
                      cv::circle(larvaROI,
                                 cv::Point2f(
                                   detected_larvae[it->first].heads.back().x+PAD,
                                   detected_larvae[it->first].heads.back().y+PAD),
                                 1,
                                 cv::Scalar(0,255,0),
                                 -1);

                      cv::circle(larvaROI,
                                 cv::Point2f(
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

      cv::imshow("Extracted Frame",frame);
      if (LRVTRACK_SAVE_PROCESSED_VIDEO!="")
        {
          vidPOut << frame;
        }

      int k;

      if (STEP==true)
        k=46;
      else
          k=cv::waitKey(1);
      if (k>=0)
        {
          if (k==32)
            {
              if(!TRACK)
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
