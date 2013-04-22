#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/features2d/features2d.hpp>
#include <opencv2/nonfree/nonfree.hpp>
#include <opencv2/nonfree/features2d.hpp>
#include <opencv2/video/background_segm.hpp>
#include <opencv2/video/tracking.hpp>
#include <math.h>
#include "cvblob.h"
#include <iostream>
#include <stdlib.h>
#include <vector>
#include <string>

#define MIN_DISCARD_DISTANCE 30
#define ROI_PADDING 7

static double VIDEO_FPS=24.0;
static unsigned int CURRENT_FRAME=0;
unsigned int LARVAE_COUNT;
std::map<unsigned int, std::vector<unsigned int> > detected_clusters;
cv::Mat frame;
cv::Mat thresholded_frame;
cv::Mat grey_frame;
cv::Mat bgFrame;
IplImage *labelImg;

//minor statistics functions
double interquantile_range(std::vector<double> sizes)
{
  double median,q1,q3;
  int middle,fq;
  std::vector<double> copy=sizes;
  sort(copy.begin(),copy.end());
  middle=copy.size()/2;
  if (copy.size()%2==0)
    {
      median=(copy[middle]+copy[middle-1])/2;
      fq=middle/2;
      if (middle%2==0)
        {
          q1=(copy[fq]+copy[fq-1])/2;
          q3=(copy[middle+fq]+copy[middle+fq-1])/2;
        }
      else
        {
          q1=copy[fq];
          q3=copy[middle+fq];
        }
    }
  else
    {
      median=copy[middle];
      fq=middle/2;
      if (middle%2==0)
        {
          q1=(copy[fq]+copy[fq-1])/2;
          q3=(copy[middle+fq]+copy[middle+fq-1])/2;
        }
      else
        {
          q1=copy[fq];
          q3=copy[middle+fq+1];
        }
    }
  return q3-q1;
}

void convertVectorToPoint2f(std::vector<cv::Point>& input, std::vector<cv::Point2f>& output)
{
  output.clear();

  for (unsigned int i = 0; i < input.size(); i++)
    {
      output.push_back(cv::Point2f((float)input.at(i).x, (float)input.at(i).y));
    }
}

void convertVectorToPoint(std::vector<cv::Point2f>& input, std::vector<cv::Point>& output)
{
  output.clear();

  for (unsigned int i = 0; i < input.size(); i++)
    {
      output.push_back(cv::Point((int)input.at(i).x, (int)input.at(i).y));
    }
}

/* Skeletonization based on the algorithm of Zhang and Suen 
 *
 * Implementation copied from here by ZachTM:
 * http://answers.opencv.org/question/3207/what-is-a-good-thinning-algorithm-for-getting-the/
 * It contains three functions:
 *  ThinSubiteration1, ThinSubiteration2 and larvae_skel (original: normalizeLetter)
 * */

class larvaSkel{
  private:
      cv::Mat skelImg;
      std::vector<cv::Point> skelPoints;

    void ThinSubiteration1(cv::Mat & pSrc, cv::Mat & pDst) {
      int rows = pSrc.rows;
      int cols = pSrc.cols;
      pSrc.copyTo(pDst);
      for(int i = 0; i < rows; i++) {
        for(int j = 0; j < cols; j++) {
          if(pSrc.at<float>(i, j) == 1.0f) {
            /// get 8 neighbors
            /// calculate C(p)
            int neighbor0 = (int) pSrc.at<float>( i-1, j-1);
            int neighbor1 = (int) pSrc.at<float>( i-1, j);
            int neighbor2 = (int) pSrc.at<float>( i-1, j+1);
            int neighbor3 = (int) pSrc.at<float>( i, j+1);
            int neighbor4 = (int) pSrc.at<float>( i+1, j+1);
            int neighbor5 = (int) pSrc.at<float>( i+1, j);
            int neighbor6 = (int) pSrc.at<float>( i+1, j-1);
            int neighbor7 = (int) pSrc.at<float>( i, j-1);
            int C = int(~neighbor1 & ( neighbor2 | neighbor3)) +
              int(~neighbor3 & ( neighbor4 | neighbor5)) +
              int(~neighbor5 & ( neighbor6 | neighbor7)) +
              int(~neighbor7 & ( neighbor0 | neighbor1));
            if(C == 1) {
              /// calculate N
              int N1 = int(neighbor0 | neighbor1) +
                int(neighbor2 | neighbor3) +
                int(neighbor4 | neighbor5) +
                int(neighbor6 | neighbor7);
              int N2 = int(neighbor1 | neighbor2) +
                int(neighbor3 | neighbor4) +
                int(neighbor5 | neighbor6) +
                int(neighbor7 | neighbor0);
              int N = cv::min(N1,N2);
              if ((N == 2) || (N == 3)) {
                /// calculate criteria 3
                int c3 = ( neighbor1 | neighbor2 | ~neighbor4) & neighbor3;
                if(c3 == 0) {
                  pDst.at<float>( i, j) = 0.0f;
                }
              }
            }
          }
        }
      }
    }


    void ThinSubiteration2(cv::Mat & pSrc, cv::Mat & pDst) {
      int rows = pSrc.rows;
      int cols = pSrc.cols;
      pSrc.copyTo( pDst);
      for(int i = 0; i < rows; i++) {
        for(int j = 0; j < cols; j++) {
          if (pSrc.at<float>( i, j) == 1.0f) {
            /// get 8 neighbors
            /// calculate C(p)
            int neighbor0 = (int) pSrc.at<float>( i-1, j-1);
            int neighbor1 = (int) pSrc.at<float>( i-1, j);
            int neighbor2 = (int) pSrc.at<float>( i-1, j+1);
            int neighbor3 = (int) pSrc.at<float>( i, j+1);
            int neighbor4 = (int) pSrc.at<float>( i+1, j+1);
            int neighbor5 = (int) pSrc.at<float>( i+1, j);
            int neighbor6 = (int) pSrc.at<float>( i+1, j-1);
            int neighbor7 = (int) pSrc.at<float>( i, j-1);
            int C = int(~neighbor1 & ( neighbor2 | neighbor3)) +
              int(~neighbor3 & ( neighbor4 | neighbor5)) +
              int(~neighbor5 & ( neighbor6 | neighbor7)) +
              int(~neighbor7 & ( neighbor0 | neighbor1));
            if(C == 1) {
              /// calculate N
              int N1 = int(neighbor0 | neighbor1) +
                int(neighbor2 | neighbor3) +
                int(neighbor4 | neighbor5) +
                int(neighbor6 | neighbor7);
              int N2 = int(neighbor1 | neighbor2) +
                int(neighbor3 | neighbor4) +
                int(neighbor5 | neighbor6) +
                int(neighbor7 | neighbor0);
              int N = cv::min(N1,N2);
              if((N == 2) || (N == 3)) {
                int E = (neighbor5 | neighbor6 | ~neighbor0) & neighbor7;
                if(E == 0) {
                  pDst.at<float>(i, j) = 0.0f;
                }
              }
            }
          }
        }
      }
    }

  public:
    std::vector <cv::Point> startPoints;
  

    larvaSkel(cv::Mat &inputarray) {
      bool bDone = false;
      int rows = inputarray.rows;
      int cols = inputarray.cols;
      cv::Mat img_thr;
      inputarray.copyTo(img_thr);
      img_thr.convertTo(img_thr,CV_32FC1);
      img_thr.copyTo(skelImg);
      
      skelImg.convertTo(skelImg,CV_32FC1);

      /// pad source
      cv::Mat p_enlarged_src = cv::Mat(rows + 2, cols + 2, CV_32FC1);
      for(int i = 0; i < (rows+2); i++) {
        p_enlarged_src.at<float>(i, 0) = 0.0f;
        p_enlarged_src.at<float>( i, cols+1) = 0.0f;
      }
      for(int j = 0; j < (cols+2); j++) {
        p_enlarged_src.at<float>(0, j) = 0.0f;
        p_enlarged_src.at<float>(rows+1, j) = 0.0f;
      }
      for(int i = 0; i < rows; i++) {
        for(int j = 0; j < cols; j++) {
          if (img_thr.at<float>(i, j) >= 20.0f) {
            p_enlarged_src.at<float>( i+1, j+1) = 1.0f;
          }
          else
            p_enlarged_src.at<float>( i+1, j+1) = 0.0f;
        }
      }

      /// start to thin
      cv::Mat p_thinMat1 = cv::Mat::zeros(rows + 2, cols + 2, CV_32FC1);
      cv::Mat p_thinMat2 = cv::Mat::zeros(rows + 2, cols + 2, CV_32FC1);
      cv::Mat p_cmp = cv::Mat::zeros(rows + 2, cols + 2, CV_8UC1);

      while (bDone != true) {
        /// sub-iteration 1
        ThinSubiteration1(p_enlarged_src, p_thinMat1);
        /// sub-iteration 2
        ThinSubiteration2(p_thinMat1, p_thinMat2);
        /// compare
        compare(p_enlarged_src, p_thinMat2, p_cmp, CV_CMP_EQ);
        /// check
        int num_non_zero = countNonZero(p_cmp);
        if(num_non_zero == (rows + 2) * (cols + 2)) {
          bDone = true;
        }
        /// copy
        p_thinMat2.copyTo(p_enlarged_src);
      }
      //Point2 *skelPointCurve=(Point2 *) malloc(rows*cols*sizeof(Point2Struct));
      //cv::GaussianBlur(p_enlarged_src,p_enlarged_src,cv::Size(2,2),0,0);
      // copy result
      int cPoints=0;
      for(int i = 0; i < rows; i++) {
        for(int j = 0; j < cols; j++) {
          skelImg.at<float>( i, j) = p_enlarged_src.at<float>( i+1, j+1);
          if (skelImg.at<float>( i, j) > 0)
          {
            skelPoints.push_back(cv::Point(j,i));
        //    skelPointCurve[cPoints].x=i;
        //    skelPointCurve[cPoints++].y=cols-j+1;
          }
        }
      }

      
      //skelImg=cv::Mat::zeros(skelImg.size(),skelImg.depth());
      //FitCurve(skelPointCurve,skelPoints.size(),0.1,skelImg);
      //free(skelPointCurve);
      //cv::imshow("Curve", skelImg);
      //cv::moveWindow("Threshold",0,100);
      //cv::moveWindow("Curve",0,100);
      //cv::moveWindow("Skeleton",0,200);


      getStartingPoints();
    }
    

    void drawSkeleton(cv::Mat &img,cv::Scalar col=cv::Scalar(0,0,0))
    {
      for (int i=0;i<skelPoints.size();i++)
      {
          cv::circle(img,
              skelPoints[i], // circle centre
              0,       // circle radius
              col, // color
              -1);              // thickness
      }
    }

    void getStartingPoints()
    {
      for (int i=0;i<skelPoints.size();i++)
      {
        int nonZeroNBs=0;
        int x=skelPoints[i].x;
        int y=skelPoints[i].y;
        if ( x-1>=0 && y+1<skelImg.rows && skelImg.at<float>(y+1,x-1)!=0)
        {
          nonZeroNBs++;
        }
        if (y+1<skelImg.rows && skelImg.at<float>(y+1,x)!=0)
        {
          nonZeroNBs++;
        }
        if (y+1<skelImg.rows && x+1<skelImg.cols && skelImg.at<float>(y+1,x+1)!=0)
        {
          nonZeroNBs++;
        }
        if ( x-1>=0 && skelImg.at<float>(y,x-1)!=0)
        {
          nonZeroNBs++;
        }
        if ( x+1<skelImg.cols && skelImg.at<float>(y,x+1)!=0)
        {
          nonZeroNBs++;
        }
        if ( x-1>=0 && y-1>=0 && skelImg.at<float>(y-1,x-1)!=0)
        {
          nonZeroNBs++;
        }
        if (y-1>=0 && skelImg.at<float>(y-1,x)!=0)
        {
          nonZeroNBs++;
        }
        if ( x+1<skelImg.cols && y-1>=0 && skelImg.at<float>(y-1,x+1)!=0)
        {
          nonZeroNBs++;
        }
        if (nonZeroNBs==1)
        {
          // We have a candidate starting point
          startPoints.push_back(cvPoint(skelPoints[i].x,skelPoints[i].y)); 
        }
      }
    }
};

/* Skeletonization Part ends here */

double getSurroundingSize(cv::Point point, cvb::CvBlob* blob)
{
  cv::Mat bigLarvaImg,larvaImg;
  int PADDR=4;
  grey_frame.copyTo(larvaImg);
  IplImage* ipl_frame = new IplImage(larvaImg);
  cv::Mat ROI=larvaImg(cv::Rect(blob->minx-ROI_PADDING,blob->miny-ROI_PADDING,blob->maxx-blob->minx+2*ROI_PADDING,blob->maxy-blob->miny+2*ROI_PADDING));
  cv::circle(ROI, cv::Point(point.x,point.y),1,cv::Scalar(0),-1);
  cv::Mat cROI(ROI.size(),ROI.depth());
  cROI=cv::Scalar(0);
  //cv::normalize(ROI, ROI, 0, 255, CV_MINMAX );
  cv::circle(cROI, cv::Point(point.x,point.y),7,cv::Scalar(255),-1);
  cv::Mat area=ROI&cROI;
  cv::imshow("Skeleton",area);
  cv::imshow("Blob",ROI);
  cv::waitKey(1);
  double nz=cv::norm(area,cv::NORM_L1);
  //double nz=cv::countNonZero(area);
  return nz;
  //return 1;
}

/*
 * Class containing information about each larva and its history on the plate.
 */
class larvaObject
{
  public:
  unsigned int start_frame;
  unsigned int cur_frame;
  unsigned int larva_ID;
  std::vector<cvb::CvBlob*> blobs; //Blob for each frame for a given larva
  double area_mean;
  double area_max;
  double area_min;
  double grey_value_mean;
  double grey_value_max;
  double grey_value_min;
  std::vector<double> centroid_speed_x;
  std::vector<double> centroid_speed_y;
  std::vector<bool> inBlob;
  std::vector<larvaSkel> lrvskels;
  std::vector<cv::Point> heads;
  std::vector<cv::Point> tails;
  //std::vector<larvaSpline*> splines;
  //std::vector<cv::Point> head;
  //std::vector<std::vector<cv::Point> > skeleton;
  larvaObject(): 
    start_frame(0),
    cur_frame(0),
    area_mean(0),
    area_max(0),
    area_min(0),
    grey_value_mean(0),
    grey_value_max(0),
    grey_value_min(0)
  {}
};

void flattenedClusters(std::vector<unsigned int> &inClusters)
{
  std::map<unsigned int, std::vector<unsigned int> >::iterator cluster=detected_clusters.begin();
  while (cluster!=detected_clusters.end())
  {
    std::vector<unsigned int>::iterator elem=(*cluster).second.begin();
    while (elem!=(*cluster).second.end())
    {
      inClusters.push_back(*elem);
      elem++;
    }
    cluster++;
  }
  std::cout << "LARVAE IN CLUSTERS: ";
  for (int i=0; i<inClusters.size() ; i++)
  {
    std::cout << inClusters[i] << " ";
  }
  std::cout << std::endl;
}

bool findInVector(std::vector<unsigned int> &flattenedCluster, unsigned int ID)
{
  for (std::vector<unsigned int>::iterator i=flattenedCluster.begin();
       i!=flattenedCluster.end();i++)
  {
    if (*i==ID)
    {
      std::cout << "ID: " << ID << " is in a cluster. Iterator: " << *i << std::endl;
      return true;
    }
  }
  std::cout << "ID: " << ID << " is not in a cluster." << std::endl;
  return false;
}

std::map<unsigned int,larvaObject> detected_larvae;

double diff(cv::Point &a, cv::Point &b)
{
  return (fabs(a.x-b.x)+fabs(a.y-b.y));
}

void findHeadTail(std::vector<cv::Point> &startPoints, larvaObject *lrv ,cv::Point &Head, cv::Point &Tail)
{
  int max=0,min=65535;
  int i=0;
  if(lrv->start_frame==CURRENT_FRAME)
  {
    for (i=0; i<startPoints.size(); i++)
    {
      unsigned int tmpsz=getSurroundingSize(startPoints[i],lrv->blobs.back());
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
    for (i=0; i<startPoints.size(); i++)
    {
      double hdiff=diff(lrv->heads.back(),startPoints[i]);
      double tdiff=diff(lrv->tails.back(),startPoints[i]);
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

void createLarvaROI(cv::Mat &frame, cv::Mat &ROI, cvb::CvBlob *blob)
{
    ROI=cv::Mat(frame, cv::Rect(blob->minx-ROI_PADDING,blob->miny-ROI_PADDING,blob->maxx-blob->minx+2*ROI_PADDING,blob->maxy-blob->miny+2*ROI_PADDING));
      
}

void createLarvaContour(cv::Mat &lrvROI, cvb::CvBlob *blob,int type=CV_8UC1)
{
    int sizes[1];
    cvb::CvContourPolygon cntPoly=*cvb::cvConvertChainCodesToPolygon(&blob->contour);
    lrvROI=cv::Mat::zeros(blob->maxy-blob->miny+(2*ROI_PADDING),blob->maxx-blob->minx+(2*ROI_PADDING),type);
    cv::Point *ContourPoints[1];
    ContourPoints[0]=(cv::Point*) malloc(blob->contour.chainCode.size()*sizeof(cv::Point));
    sizes[0]=cntPoly.size();
    for (int i=0;i<=cntPoly.size();i++)
    {
      ContourPoints[0][i].x=cntPoly[i].x-blob->minx+ROI_PADDING;
      ContourPoints[0][i].y=cntPoly[i].y-blob->miny+ROI_PADDING;
    }
    if(lrvROI.channels()==3)
    {
      cv::fillPoly(lrvROI,(const cv::Point**) ContourPoints,sizes,1,cv::Scalar(255,255,255));
    }
    else
    {
      cv::fillPoly(lrvROI,(const cv::Point**) ContourPoints,sizes,1,cv::Scalar(255));
    }
}

void updateLarvae(cvb::CvBlobs &In, cvb::CvBlobs &Prev)
{
  cvb::CvBlobs::iterator it=In.begin();
  std::vector<unsigned int> larvaeInClusters;
  flattenedClusters(larvaeInClusters);
  while (it!=In.end())
  {
    unsigned int ID=(*it).first;
    cvb::CvBlob* blob=(*it).second;
    //cv::Mat larvaROI(grey_frame, cv::Rect(blob->minx-ROI_PADDING,blob->miny-ROI_PADDING,blob->maxx-blob->minx+2*ROI_PADDING,blob->maxy-blob->miny+2*ROI_PADDING));
    cv::Mat larvaROI;
    createLarvaContour(larvaROI,blob);
    //cv::imshow("Skeleton",larvaROI);
    std::map<unsigned int,larvaObject>::iterator curLarva;
    if ((curLarva=detected_larvae.find(ID))==detected_larvae.end())
    {
      larvaObject newLarva;
      newLarva.start_frame=CURRENT_FRAME;
      newLarva.larva_ID=ID;
      newLarva.blobs.push_back(blob);
      newLarva.area_mean=blob->area;
      newLarva.area_max=newLarva.area_min=blob->area;
      newLarva.area_min=newLarva.area_min=blob->area;
      newLarva.grey_value_mean = 0; //TODO
      newLarva.grey_value_max = 0; //TODO
      newLarva.grey_value_min = 0; //TODO
      newLarva.centroid_speed_x.push_back(0);
      newLarva.centroid_speed_y.push_back(0);
      newLarva.inBlob.push_back(false);
      larvaSkel newLarvaSkel(larvaROI);
      newLarva.lrvskels.push_back(newLarvaSkel);
      cv::Point Head,Tail;
      findHeadTail(newLarva.lrvskels.back().startPoints,&newLarva,Head,Tail);
      newLarva.heads.push_back(Head);
      newLarva.tails.push_back(Tail);
      //cv::imshow("Skeleton", newLarvaSkel.skelImgColor);
      //cv::waitKey(1);
      //larvaSpline *newLarvaSpline=new larvaSpline(larvaROI);
      //newLarva.splines.push_back(newLarvaSpline);
      /*
      std::cout << "UL: Spline for ID: " << ID  << std::endl;
      std::cout << "UL: Spline: " << newLarvaSpline->bsplines[0] << std::endl;
      std::cout << "UL: Size: " << newLarvaSpline->bspline[0]->size << std::endl;
      std::cout << "UL: Stride: " << newLarvaSpline->bspline[0]->stride << std::endl;
      std::cout << "UL: Data: " << newLarvaSpline->bspline->data << std::endl;
      std::cout << "UL: Block: " << newLarvaSpline->bspline->block << std::endl;
      std::cout << "UL: Owner: " << newLarvaSpline->bspline->owner << std::endl;
      std::cout <<  std::endl;
      */
      detected_larvae[ID]=newLarva;
      //std::cout << "New Larva: " << ID << " area: " << blob->area << std::endl;
    }
    else
    {
      larvaObject *cur_larva=&(*curLarva).second;
      cvb::CvBlob *preBlob = cur_larva->blobs.back();
      cur_larva->larva_ID=ID;
      cur_larva->blobs.push_back(blob);
      larvaSkel newLarvaSkel(larvaROI);
      cur_larva->lrvskels.push_back(newLarvaSkel);
      if (!findInVector(larvaeInClusters,ID))
      {
        cur_larva->area_mean=((cur_larva->area_mean+blob->area)/2);
      
        if (cur_larva->area_max < blob->area)
        {
          cur_larva->area_max=blob->area;
        }
        if (cur_larva->area_min > blob->area)
        {
          cur_larva->area_min=blob->area;
        }
	cur_larva->centroid_speed_x.push_back((blob->centroid.x - preBlob->centroid.x)/VIDEO_FPS);
	cur_larva->centroid_speed_y.push_back((blob->centroid.y - preBlob->centroid.y)/VIDEO_FPS);
        //larvaSpline *newLarvaSpline=new larvaSpline(larvaROI);
        //cur_larva->splines.push_back(newLarvaSpline);
        //std::cout << "Updated Larva: " << ID << " area: " << blob->area << " area_mean: " << cur_larva->area_mean << " inblob: NO" << std::endl ;
        cv::Point Head,Tail;
        findHeadTail(cur_larva->lrvskels.back().startPoints,cur_larva,Head,Tail);
        cur_larva->heads.push_back(Head);
        cur_larva->tails.push_back(Tail);
      }
      else
      {
        //std::cout << "Updated Larva: " << ID << " area: " << blob->area << " area_mean: " << cur_larva->area_mean << " inblob: YES" << std::endl;
        //std::cout << "Updated Larva: " << detected_clusters[ID][1]<< " area: " << blob->area << " area_mean: " << cur_larva->area_mean << " inblob: YES" << std::endl;
        detected_larvae[detected_clusters[ID][1]].blobs.push_back(blob);
        cur_larva->heads.push_back(cvPoint(blob->centroid.x,blob->centroid.y));
        cur_larva->tails.push_back(cvPoint(blob->centroid.x,blob->centroid.y));
      }
      cur_larva->grey_value_mean = 0; //TODO
      cur_larva->grey_value_max = 0; //TODO
      cur_larva->grey_value_min = 0; //TODO
      cur_larva->inBlob.push_back(false);
    }
    it++;
  }
}

// Matching of diverging larvae
void diverge_match(unsigned int &Larva1, unsigned int &Larva2, cvb::CvBlob  *curLarva1, cvb::CvBlob *curLarva2)
{
  double SizeDiff=fabs(detected_larvae[Larva1].area_mean - curLarva1->area) - 
           fabs(detected_larvae[Larva1].area_mean-curLarva2->area);
  
  // detected_larvae[Larva1].centroid_speed_x
  //SpeedDiff=
  if (fabs(detected_larvae[Larva1].area_mean-curLarva1->area) <
      fabs(detected_larvae[Larva1].area_mean-curLarva2->area))
  {
    Larva1=curLarva1->label;
    Larva2=curLarva2->label;
  }
  else
  {
    Larva1=curLarva2->label;
    Larva2=curLarva1->label;
  }
}

void larvae_track(cvb::CvBlobs &In,cvb::CvBlobs &Prev,cvb::CvBlobs &out)
{

  // Map to keep track of the larvae that were tracked succesfully so far
  std::map<unsigned int,std::vector<unsigned int> > used_map;
  
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
      // Stores the minimal distance found (actually the square of the distance).
      double min=200;

      // Here we check if the *prevIt was already assigned a match.
      // This can occur when diverging where the two diverging worms are
      // both assigned at the same time.
      if(out.find((*prevIt).first)!=out.end())
      {
        prevIt++;
        continue;
      }

      // Here we start looking for the blob in "In" (i.e. the latest larvae set)
      // that is closest to the preBlob. Typically distances of matching larvae should be
      // too small to calculate even with double precision.
      cvb::CvBlobs::iterator it=In.begin();
      while (it!=In.end())
        {
          cvb::CvBlob *blob;
          blob=(*it).second;
          // MIN_DISCARD_DISTANCE is used to quickly exclude larvae that are
          // too far to be considered
          if (((XVAL=fabs(blob->centroid.x - preBlob->centroid.x)) < MIN_DISCARD_DISTANCE ) &&
              ((YVAL=fabs(blob->centroid.y - preBlob->centroid.y)) < MIN_DISCARD_DISTANCE ))
            {
              // we only use the square of the distance (sqrt is expensive)
              cur=XVAL+YVAL;
              if (cur < min)
                {
                  min=cur;
                  minLabel=(*it).first;
                }
            }
          it++;
        }

      //std::cout << "minimum distance PreID:" << (*prevIt).first << " NewID: " << minLabel << " = " << min << std::endl;
      // We went through all the current batch of larvae and found minimum distance min
      // for the larva with label minLabel
      if (min<=1.5)
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
      else if (min<150)
        {
          // Rather large jump. Indicates:
          // 1) larvae converged
          // 2) larvae diverged
          // 3) TODO: Larva dissapeared and matched wrongly with another...

          // Check for converging
          // We check if more than one larvae match the new centroid and
          //  whether the middle of their centroids matches the centroid we have.
          
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
            double newx = ((*prevIt).second->centroid.x + Prev[used_map[minLabel][0]]->centroid.x);
            double newy = ((*prevIt).second->centroid.y + Prev[used_map[minLabel][0]]->centroid.y);
            double XDIF,YDIF;
            // TODO: Use correct centroid calculation based on area not just middle
            if (((XDIF=fabs(2*blob->centroid.x - newx)) < 10 ) &&
                ((YDIF=fabs(2*blob->centroid.y - newy)) < 10))
            {
              // We're in luck! The center's match :) We have a cluster.
              used_map[minLabel].push_back((*prevIt).first);
	      // Add both larvae to the detected clusters
              detected_clusters[blob->label].push_back(used_map[minLabel][0]);
              detected_clusters[blob->label].push_back(used_map[minLabel][1]);
	      // Create a second "instance" of the blob to represent both larvae
	      //out[(*prevIt).first]=In[minLabel];
            }
            else
            {
              // May be that one larva dissapeared and matched wrongly with another nearby
              // In this case we should only keep the one that is closest!
              std::cout << "PROBLEM!!! : Check how to resolve these cases where centers after do not match" << std::endl;
            }

          }
          // Check for diverging:
          // If the detected_clusters map contains a cluster for the previous blob with the one
          // we are examining we might be looking into larvae that have separated.
          // TODO: nested clusters would also be nice...
          else if ((cluster=detected_clusters.find((*prevIt).first))!=detected_clusters.end())
          {
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
                if (((XVAL=fabs(blob->centroid.x - preBlob->centroid.x)) < 90) &&
                    ((YVAL=fabs(blob->centroid.y - preBlob->centroid.y)) < 90))
                {
                  cur=XVAL+YVAL;
                  if (cur < min)
                  {
                    min=cur;
                    secondDivergent=(*it).first;
                  }
                }
              }
              it++;
            }
            if(secondDivergent==0)
            {
              std::cout << "PROBLEM!!! : Did not find second divergent." << std::endl;
            }
            double newx = (In[minLabel]->centroid.x + In[secondDivergent]->centroid.x);
            double newy = (In[minLabel]->centroid.y + In[secondDivergent]->centroid.y);
            // TODO: Use correct centroid calculation based on area not just middle
            if (((XVAL=fabs(2*(*prevIt).second->centroid.x - newx)) < 20 ) &&
                ((YVAL=fabs(2*(*prevIt).second->centroid.y - newy)) < 20))
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
              detected_clusters.erase(cluster);
            }
            else
            {
              std::cout << "PROBLEM!!! : Check how to resolve these cases where centers do not match" << std::endl;
            }
          }
          else
          {
            // Here we are in two cases:
            // 1) First larva spotted belonging to a cluster
            // 2) Divergence of larvae that were clustered from the start <TODO>
            // 3) A small jump. Perhaps framerate related... <TODO>
              used_map[minLabel].push_back((*prevIt).first);
              out[(*prevIt).first]=In[minLabel];
              out[(*prevIt).first]->label=(*prevIt).first;
          }
          //out[(*prevIt).first]=In[minLabel];
          //out[(*prevIt).first]->label=(*prevIt).first;
          //used_map[minLabel]++;
          //detected_clusters[minLabel].push_back((*prevIt).first);
        }
      prevIt++;
    }

  cvb::CvBlobs::iterator it=In.begin();
  while (it!=In.end())
    {
      int ID=(*it).first;
      if (used_map.find(ID)==used_map.end())
        {
            out[++LARVAE_COUNT]=(*it).second;
            out[LARVAE_COUNT]->label=LARVAE_COUNT;
        }
      it++;
    }
  updateLarvae(out,Prev);

}


int main(int argv, char* argc[])
{
  bool SHOWTAGS=true;
  bool TRACK=true;
  bool STEP=false;
  //cv::VideoCapture capture(CV_CAP_DC1394);
  cv::VideoCapture capture("/Users/epaisios/Desktop/LarvaeCapture201302211115.mp4");
  //cv::VideoCapture capture("/Users/alasondro/Desktop/LarvaeCapture201302211115.mp4");
  //cv::VideoCapture capture("/Users/alasondro/Desktop/LarvaeCapture201302211054.mp4");
  //cv::VideoCapture capture("/Users/epaisios/Desktop/LarvaeCapture201302201531.mp4");

  if (capture.get(CV_CAP_PROP_FPS)==0)
  {
    VIDEO_FPS=24.1;
  }
  else
  {
    if ( capture.get(CV_CAP_PROP_FPS) > 26 )
    {
      VIDEO_FPS=capture.get(CV_CAP_PROP_FPS)/2;
    }
    else
    {
      VIDEO_FPS=capture.get(CV_CAP_PROP_FPS);
    }
  }
  cv::Mat preframe;
  cv::Mat grey_bgFrame;
  unsigned int fg_threshold=10;
  unsigned int thresholdlow=35;
  unsigned int thresholdhigh=255;
  unsigned int cannylow=70;
  unsigned int cannyhigh=255;
  unsigned int contourMin=25;
  unsigned int contourMax=250;
  LARVAE_COUNT=0;

  if (!capture.isOpened())
    return 1;
  // Get the frame rate
  //double rate=capture.get(CV_CAP_PROP_FPS);
  //std::cout << rate << std::endl;
  //std::cout << capture.get(CV_CAP_PROP_FRAME_WIDTH) << std::endl;
  //std::cout << capture.get(CV_CAP_PROP_FRAME_HEIGHT) << std::endl;
  //std::cout << capture.get(CV_CAP_PROP_MODE) << std::endl;
  //std::cout << capture.get(CV_CAP_PROP_FORMAT) << std::endl;

  bool stop(false);
  //memset(inactive, 0, sizeof(int)*50);
  cv::namedWindow("Extracted Frame");
  //cv::namedWindow("Skeleton");
  //cv::namedWindow("Blob");
  //cv::namedWindow("Threshold");
  //cv::moveWindow("Blob",0,100);
  //cv::moveWindow("Threshold",0,200);
  //cvCreateTrackbar("FOO", "Extracted Frame" , (int *) &LARVAE_COUNT, 30 , NULL);

  //cv::namedWindow("Label Frame");

  // Delay between each frame in ms
  // corresponds to video frame rate
  //int delay= 1;
  // for all frames in video

  capture.read(bgFrame);
  cvtColor(bgFrame,grey_bgFrame,CV_BGR2GRAY);
  std::vector<cv::Vec3f> circles;
  cv::HoughCircles(grey_bgFrame, circles, CV_HOUGH_GRADIENT,
                   1,   // accumulator resolution (size of the image / 2)
                   850,  // minimum distance between two circles
                   50, // Canny high threshold
                   200, // minimum number of votes
                   bgFrame.rows/3, bgFrame.rows/1.9); // min and max radiusV

  std::vector<cv::Vec3f>::
  const_iterator itc= circles.begin();
  while (itc!=circles.end())
    {
      cv::circle(bgFrame,
                 cv::Point((*itc)[0], (*itc)[1]), // circle centre
                 (*itc)[2]/2,       // circle radius
                 cv::Scalar(), // color
                 -1);              // thickness
      ++itc;
    }
  cvtColor(bgFrame,grey_bgFrame,CV_BGR2GRAY);
  grey_bgFrame.convertTo(bgFrame,CV_32F);
  cvb::CvBlobs preBlobs;
  while (!stop)
  {
    // read next frame if any
    if (!capture.read(frame))
      break;

    if(TRACK)
    {

      CURRENT_FRAME=capture.get(CV_CAP_PROP_POS_FRAMES);
      cv::Mat fg_frame;
      cv::Mat image;
      cv::Mat processed_frame;
      cv::Mat working_frame;
      cv::Mat fgROI;
      cv::Mat fg_image;
      cv::Mat mask;

      //cv::GaussianBlur(frame, image, cv::Size(0, 0), 2);
      //cv::addWeighted(frame, 1.5, image, -0.5, 0, image);
      //frame=image;
      cvtColor(frame,grey_frame,CV_BGR2GRAY);
      cv::addWeighted(grey_frame, 1.0, grey_bgFrame, -3.0, 0.0, fg_frame);
      fgROI=cv::Mat::zeros(fg_frame.rows , fg_frame.cols,fg_frame.depth());
      cv::circle(fgROI, cv::Point(circles[0][0],circles[0][1]),int(circles[0][2]/1.0),cv::Scalar(255),-1);
      cv::circle(frame, cv::Point(circles[0][0],circles[0][1]),int(circles[0][2]/1.0),cv::Scalar(0,255,0),1);
      fg_image=fg_frame&fgROI;
      cv::Mat fg_image_norm;
      //cv::normalize(fg_image,fg_image_norm,0,255,cv::NORM_MINMAX);
      //cv::imshow("Skeleton",fg_image_norm);
      //

      cv::inRange(fg_image,cv::Scalar(0),cv::Scalar(30),mask);
      cv::bitwise_not(mask,mask);
      fg_image.copyTo(image,mask);
      cv::equalizeHist(image,fg_image);
      cv::normalize(fg_image,fg_image_norm,0,255,cv::NORM_MINMAX);
      //cv::imshow("Skeleton",fg_image_norm);

      cv::threshold(fg_image_norm,
          thresholded_frame,
          thresholdlow,
          thresholdhigh,
          cv::THRESH_BINARY);
      //cv::dilate(thresholded_frame,working_frame,cv::Mat(),cv::Point(-1,-1),1);
      //thresholded_frame=working_frame;
      //frame.copyTo(preframe);
      cvb::CvBlobs blobs;
      //IplImage ipl_preframe = preframe;
      IplImage ipl_thresholded = thresholded_frame;
      labelImg=cvCreateImage(cvGetSize(&ipl_thresholded), IPL_DEPTH_LABEL, 1);
      //cv::imshow("Skeleton",cv::Mat(labelImg));
      unsigned int result=cvLabel(&ipl_thresholded, labelImg, blobs);
      //std::cout << "###########  Total Blobs: " << blobs.size() << std::endl;
      cvb::cvFilterByArea(blobs, 31, 250);
      cvb::CvBlobs tracked_blobs;
      cvb::CvBlobs blob1;
      cvb::CvBlobs::iterator it=blobs.begin();

      if(preBlobs.size()>0)
      {
        larvae_track(blobs,preBlobs,tracked_blobs);
        preBlobs.clear();
        preBlobs=tracked_blobs;
        if(tracked_blobs.find(7)!=tracked_blobs.end())
          blob1[7]=tracked_blobs[7];
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
        LARVAE_COUNT=preBlobs.size();
        //blob1[8]=blobs[8];
      }
      //updateLarvae(preBlobs);
      //if(tracked_blobs.find(7)!=tracked_blobs.end())

      if (SHOWTAGS!=0)
      {
        it=tracked_blobs.begin();
        while (it!=tracked_blobs.end())
        {
          std::stringstream sstm;
          sstm << (*it).first;
          cvb::CvBlob *blob=(*it).second;
          cv::putText(frame,
              sstm.str(),
              cv::Point2d(blob->centroid.x+10,blob->centroid.y+10),
              cv::FONT_HERSHEY_PLAIN,
              0.8,
              cv::Scalar(255,255,255),
              1,
              CV_AA);
          cv::Mat larvaROI(frame, cv::Rect(blob->minx-ROI_PADDING,blob->miny-ROI_PADDING,blob->maxx-blob->minx+2*ROI_PADDING,blob->maxy-blob->miny+2*ROI_PADDING));
          
          detected_larvae[it->first].lrvskels.back().drawSkeleton(larvaROI,cv::Scalar(0,0,255));
          cv::circle(frame,
              cv::Point2d(blob->centroid.x,blob->centroid.y), // circle centre
              1,       // circle radius
              cv::Scalar(255,0,0), // color
              -1);              // thickness
          /*
          detected_larvae[it->first].lrvskels.back().drawSkeleton(larvaROI,cv::Scalar(0,0,255));
          cv::circle(frame,
              cv::Point2d(blob->contour.startingPoint.x,blob->contour.startingPoint.y), // circle centre
              1,       // circle radius
              cv::Scalar(0,255,0), // color
              -1);              // thickness
              */
          cv::circle(larvaROI,
              cv::Point2d(detected_larvae[it->first].heads.back().x,detected_larvae[it->first].heads.back().y), // circle centre
              1,       // circle radius
              cv::Scalar(0,255,0), // color
              -1);              // thickness
        /*
          cv::circle(frame,
              cv::Point2d(detected_larvae[it->first].tails.back().x+blob->minx,detected_larvae[it->first].tails.back().y+blob->miny), // circle centre
              1,       // circle radius
              cv::Scalar(0,0,255), // color
              -1);              // thickness
          
          */
          it++;
        }
      }

      //cvUpdateTracks(blobs, tracks, 2, 900);
      cvReleaseImage(&labelImg);

      //if(preframe.empty())
      //  frame.copyTo(preframe);
      //cv::imshow("Previous Frame",preframe);

      //frame.copyTo(preframe);
      // introduce a delay
      // or press key to stop
    }
    cv::imshow("Extracted Frame",frame);
    int k;

    if (STEP==true)
      k=46;
    else
      k=cv::waitKey(1);
    if (k>=0)
    {
      if (k==32)
      {
        while (cv::waitKey(1)!=32){
          //No OP
        }
      }
      if (k==46)
      {
        STEP=true;
        while ((k=cv::waitKey(1))>=127 || k<0){
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
}
