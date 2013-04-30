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
#include <queue>

#define MIN_DISCARD_DISTANCE 30
#define ROI_PADDING 2
#define MAX_HORIZ_RESOLUTION 22000 //Pixels

/*#define LARVAE_HASH(l,s,d)\
 ( Wlength*l+Wsize*s+Wdir*d)
*/
#define LARVAE_HASH(l,s)\
 ( Wlength*l+Wsize*s)

#ifndef NO_SSE
#define ltsqrt SSESqrt_Recip_Times_X
#else
#define ltsqrt sqrt
#endif

static double VIDEO_FPS=24.0;
static unsigned int CURRENT_FRAME=0;
unsigned int LARVAE_COUNT;
std::map<unsigned int, std::vector<unsigned int> > detected_clusters;
cv::Mat frame;
cv::Mat thresholded_frame;
cv::Mat grey_frame;
cv::Mat bgFrame;
IplImage *labelImg;

double Wlength=1.0;
double Wsize=0.2;


typedef std::pair<cv::Point_<int>,cv::Point_<int> > PointPair;
typedef std::unordered_map<PointPair, double> DistanceMap;

// SSE Optimized reciprocal sqrt and mutliplies it by pIn to get the sqrt
// Optimized sqrt from http://stackoverflow.com/questions/1528727/why-is-sse-scalar-sqrtx-slower-than-rsqrtx-x
inline void SSESqrt_Recip_Times_X( float * pOut, float * pIn )
{
  __m128 in = _mm_load_ss( pIn );
  _mm_store_ss( pOut, _mm_mul_ss( in, _mm_rsqrt_ss( in ) ) );
  // compiles to movss, movaps, rsqrtss, mulss, movss
}

namespace std{
template <typename T >
struct hash<cv::Point_<T> > {
  public:
    size_t operator()(cv::Point_<T> a) const throw() {
      size_t h = (47+ a.x )*47+a.y;
      return h;
    }
};


template <typename T >
struct hash<pair<cv::Point_<T>, cv::Point_<T> > > {
  public:
    size_t operator()(pair<cv::Point_<T>,cv::Point_<T> > a) const throw() {
      size_t h = (997 + MAX(hash<cv::Point_<T> >()(a.first),
                           hash<cv::Point_<T> >()(a.second)))*997+
                       MIN(hash<cv::Point_<T> >()(a.first),
                           hash<cv::Point_<T> >()(a.second));
      return h;
    }
};
}

void blobToPointVector(cvb::CvBlob &p,std::vector<cv::Point> &points)
{
      cvb::CvContourPolygon *cntPoly=
                             cvb::cvConvertChainCodesToPolygon(&p.contour);
      cvb::CvContourPolygon::iterator a=cntPoly->begin();
      while(a!=cntPoly->end())
      {
        points.push_back(*a++);
      }
      delete cntPoly;
}

class LarvaDistanceMap{
  private:
    std::vector<double> distances; 
    cv::Mat px;
    cv::Mat py;
  public:
    double MaxDist;
    PointPair MaxDistPoints;
    std::vector<cv::Point> &points;
    class my2ndPoint 
    {
      private:
        LarvaDistanceMap &parent;
        int p1;
      public:
      my2ndPoint(LarvaDistanceMap& p, int fstPoint):parent(p),p1(fstPoint){}

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

    LarvaDistanceMap(std::vector<cv::Point> &ps):points(ps){
      distances.reserve(points.size()*points.size());
      /*
      for (int i=0;i<=ps.size();i++)
      {
        px=ps[i].x;
        py=ps[i].y;
      }*/
    }
    void getPxPy(cv::Mat &x,cv::Mat &y){
      x=px;
      y=py;
    }
    void getDistances(cv::Point p1)
    {
    }

};

void createLarvaROI(cv::Mat &frame, cv::Mat &ROI, cvb::CvBlob &blob)
{
    ROI=cv::Mat(frame, 
                cv::Rect(blob.minx-ROI_PADDING,
                         blob.miny-ROI_PADDING,
                         blob.maxx-blob.minx+2*ROI_PADDING,
                         blob.maxy-blob.miny+2*ROI_PADDING
                         )
              );
    //cv::normalize(ROI,ROI,0,255,cv::NORM_MINMAX);  
}

void createLarvaContour(cv::Mat &lrvROI, 
                        cvb::CvBlob &blob,
                        int type=CV_8UC1)
{
    int sizes[1];
    cvb::CvContourPolygon *cntPoly=
                             cvb::cvConvertChainCodesToPolygon(&blob.contour);

    lrvROI=cv::Mat::zeros(blob.maxy-blob.miny+(2*ROI_PADDING),
                          blob.maxx-blob.minx+(2*ROI_PADDING),
                          type);

    cv::Point *ContourPoints[1];
    ContourPoints[0]=(cv::Point*) malloc(
                                blob.contour.chainCode.size()*sizeof(cv::Point)
                                  );
    sizes[0]=cntPoly->size();
    for (int i=0;i<cntPoly->size();i++)
    {
      ContourPoints[0][i].x=(*cntPoly)[i].x-blob.minx+ROI_PADDING;
      ContourPoints[0][i].y=(*cntPoly)[i].y-blob.miny+ROI_PADDING;
    }
    cv::fillPoly(lrvROI,
        (const cv::Point**) ContourPoints,
        sizes,
        1,
        cv::Scalar(255)
        );
    free(ContourPoints[0]);
    delete(cntPoly);

}

void createLarvaContourPoints(cv::Mat &lrvROI, 
                              cvb::CvBlob &blob,
                              int type=CV_8UC1)
{
    int sizes[1];
    cvb::CvContourPolygon *cntPoly=
                             cvb::cvConvertChainCodesToPolygon(&blob.contour);
    lrvROI=cv::Mat::zeros(blob.maxy-blob.miny+(2*ROI_PADDING),
                          blob.maxx-blob.minx+(2*ROI_PADDING),
                          type);
    cv::Point *ContourPoints[1];
    ContourPoints[0]=(cv::Point*) malloc(
                               blob.contour.chainCode.size()*sizeof(cv::Point)
        );

    sizes[0]=cntPoly->size();
    for (int i=0;i<cntPoly->size();i++)
    {
      ContourPoints[0][i].x=(*cntPoly)[i].x-blob.minx+ROI_PADDING;
      ContourPoints[0][i].y=(*cntPoly)[i].y-blob.miny+ROI_PADDING;
      cv::circle(lrvROI,
          ContourPoints[0][i], // circle centre
          0,       // circle radius
          cv::Scalar(255), // color
          -1);              // thickness
    }
    free(ContourPoints[0]);
    delete(cntPoly);
}


inline void lBFS(int p1, 
          std::vector<cv::Point> &Points ,
          LarvaDistanceMap &Distances,
          cv::Point MidPoint
          )
{
  std::queue<int> Q;
  Q.push(p1);
  double dst_p1_MP = (Points[p1].x-MidPoint.x)*(Points[p1].x-MidPoint.x)+
                     (Points[p1].y-MidPoint.y)*(Points[p1].y-MidPoint.y);
  double MAX=Distances.MaxDist;
  while (!Q.empty())
  {
    int cur=Q.front();
    Q.pop();
    for(int i=0; i<Points.size() ; i++)
    {
      PointPair cur_pi = PointPair(Points[cur],Points[i]);
      PointPair p1_pi= PointPair(Points[p1],Points[i]);
      PointPair p1_cur= PointPair(Points[p1],Points[cur]);
      double dst_pi_MP=(Points[i].x-MidPoint.x)*(Points[i].x-MidPoint.x)+
                       (Points[i].y-MidPoint.y)*(Points[i].y-MidPoint.y);
      if (Distances[cur][i]>0)
      {
        if (Distances[p1][i]<0)
          Q.push(i);
        double newDst = Distances[cur][i]+Distances[p1][cur] + 
                        2*(Distances[cur][i]*Distances[p1][cur]);
        float multres=dst_pi_MP * dst_p1_MP;
        float sqrtres[1];
        ltsqrt(sqrtres,&multres);
        double MPDst = dst_pi_MP+dst_p1_MP + 2*sqrtres[0];
        if (newDst < MPDst)
          newDst=MPDst;
        if (Distances[p1][i]>newDst)
        {
          Distances[p1][i]=newDst;
          Distances[i][p1]=newDst;
        }
        if (Distances[p1][i]>MAX)
        {
          MAX=Distances[p1][i];
          Distances.MaxDistPoints=p1_pi;
          Distances.MaxDist=MAX;
        }
      }
    }
  }
}

void computeInnerDistances(cvb::CvBlob &blob,
                           LarvaDistanceMap &Distances,
                           cv::Point &MidPoint)
{
  std::vector<cv::Point> Points;
  cv::Mat contour;
  createLarvaContour(contour,blob);
  cv::Mat workingContour;
  contour.copyTo(workingContour);
  std::vector<cv::Point> &points=Distances.points;
  int origNZ=countNonZero(contour);
  double MAX=0;
  for (int i=0;i<points.size();i++)
  {
    cv::Point p1=cv::Point(points[i].x,points[i].y);
    Distances[i][i]=0;
    for (int j=i+1;j<points.size();j++)
    {
      cv::Point p2(points[j].x,points[j].y);
      cv::line(workingContour,
               p1,
               p2,
               cv::Scalar(255));
      if (countNonZero(workingContour)>origNZ)
      {
        Distances[i][j]=-1;
        Distances[j][i]=-1;
      }
      else
      {
        PointPair p1p2 = PointPair(p1,p2);
        double xdiff=points[i].x-points[j].x;
        double ydiff=points[i].y-points[j].y;
        double sqDst=xdiff*xdiff+ydiff*ydiff;
        Distances[i][j]=sqDst;
        Distances[j][i]=sqDst;

        if (MAX<Distances[i][j])
        {
          MAX=Distances[i][j];
          Distances.MaxDistPoints=p1p2;
          Distances.MaxDist=MAX;
        }
      }
      contour.copyTo(workingContour);
    }
  }

  for (int i=0;i<points.size();i++)
  {
    cv::Point p1(points[i].x,points[i].y);
    lBFS(i,points,Distances,MidPoint);
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
    cv::Point MidPoint;
    bool emptySkeleton;
    larvaSkel(bool emptySkel):emptySkeleton(true) {}
    larvaSkel(cv::Mat &inputarray, cv::Point &centroid):emptySkeleton(false) {
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

      int mini;
      double min=65535;
      for (int i=0; i<skelPoints.size(); i++)
      {
        double xdiff=fabs(skelPoints[i].x-centroid.x);
        double ydiff=fabs(skelPoints[i].y-centroid.y);
        double dst=xdiff+ydiff;
        if (dst<min)
        {
          mini=i;
          min=dst;
        }
      }
      MidPoint=skelPoints[mini];
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
};

/* Skeletonization Part ends here */

double getSurroundingSize(cv::Point &point, cvb::CvBlob &blob)
{
  cv::Mat larvaImg,lrvROI;
  grey_frame.copyTo(larvaImg);

  cv::Mat ROI=larvaImg(cv::Rect(blob.minx-ROI_PADDING,
                                blob.miny-ROI_PADDING,
                                blob.maxx-blob.minx+2*ROI_PADDING,
                                blob.maxy-blob.miny+2*ROI_PADDING)
                      );
  createLarvaContour(lrvROI, blob);
  //cv::dilate(lrvROI,lrvROI,cv::Mat(),cv::Point(-1,-1),1);
  ROI=ROI&lrvROI;
  cv::Mat cROI(ROI.size(),ROI.depth());
  cROI=cv::Scalar(0);
  //cv::normalize(ROI, ROI, 0, 255, CV_MINMAX );
  //cv::circle(ROI, cv::Point(point.x,point.y),0,cv::Scalar(255),-1);
  cv::circle(cROI, cv::Point(point.x,point.y),4,cv::Scalar(255),-1);
  cv::Mat area=ROI&cROI;
  /*
  cv::imshow("Contour",lrvROI);
  cv::imshow("Blob",ROI);
  cv::imshow("Skeleton",area);
  cv::waitKey(1);
  */
  double nz=cv::norm(area,cv::NORM_L1);
  //double nz=cv::countNonZero(area);
  return nz;
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
  std::vector<cvb::CvBlob> blobs; //Blob for each frame for a given larva
  double area_mean;
  double area_max;
  double area_min;
  double grey_value_mean;
  double grey_value_max;
  double grey_value_min;
  double length_mean;
  double length_max;
  double length_min;
  std::vector<cv::Point> centroids;
  std::vector<double> centroid_speed_x;
  std::vector<double> centroid_speed_y;
  std::vector<bool> inBlob;
  std::vector<larvaSkel> lrvskels;
  std::vector<cv::Point> heads;
  std::vector<cv::Point> tails;
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
  std::map<unsigned int, std::vector<unsigned int> >::iterator cluster;
  cluster=detected_clusters.begin();
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
      std::cout << "ID: " << ID << " is in a cluster. Iterator: " 
        << *i << std::endl;
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
      lrv.inBlob.back()==true || 
      force_SurroundingValSearch ||
      (lrv.heads.back().x == 0 && lrv.heads.back().y==0 && 
       lrv.tails.back().x == 0 && lrv.tails.back().y==0)
      )
  {
    for (i=0; i<startPoints.size(); i++)
    {
      cv::Point p=startPoints[i];
      cvb::CvBlob blob=lrv.blobs.back();
      double tmpsz=getSurroundingSize(p,blob);
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
  flattenedClusters(larvaeInClusters);
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

      // Set the frame of it's existence
      newLarva.start_frame=CURRENT_FRAME;

      // Give the larva the necessary ID
      newLarva.larva_ID=ID;

      // Add the blob of the larva to its blob history
      newLarva.blobs.push_back(blob);

      // Initialize the area_mean of the larva
      newLarva.area_mean=blob.area;

      // Initialize the max and min areas as well
      newLarva.area_max=newLarva.area_min=blob.area;
      newLarva.area_min=newLarva.area_min=blob.area;


      newLarva.grey_value_mean = 0; //TODO
      newLarva.grey_value_max = 0; //TODO
      newLarva.grey_value_min = 0; //TODO
      
      // Initialize the speed to 0
      newLarva.centroid_speed_x.push_back(0);
      newLarva.centroid_speed_y.push_back(0);

      // State that the larva is not in a blob
      newLarva.inBlob.push_back(false);
      
      // Create a skeleton of the larva 
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
      LarvaDistanceMap Distances(cntPoints);
      computeInnerDistances(blob,Distances,newLarvaSkel.MidPoint);
      newLarva.length_mean = Distances.MaxDist;
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
      findHeadTail(startPoints,newLarva,Head,Tail);
      newLarva.heads.push_back(Head);
      newLarva.tails.push_back(Tail);

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
      larvaSkel newLarvaSkel(larvaROI,centroid);
      cur_larva.centroids.push_back(centroid);
      cur_larva.lrvskels.push_back(newLarvaSkel);

      
      // Look if larva was found in a blob
      if (!findInVector(larvaeInClusters,ID))
      {
        // If not then:
        //  Update area values for larva.
        cur_larva.area_mean=((cur_larva.area_mean+blob.area)/2);
      
        if (cur_larva.area_max < blob.area)
        {
          cur_larva.area_max=blob.area;
        }
        if (cur_larva.area_min > blob.area)
        {
          cur_larva.area_min=blob.area;
        }
        // Update centroid_speeds (in pixel per second per axis)
	cur_larva.centroid_speed_x.push_back(
                            (blob.centroid.x - preBlob.centroid.x)/VIDEO_FPS);
	cur_larva.centroid_speed_y.push_back(
                            (blob.centroid.y - preBlob.centroid.y)/VIDEO_FPS
            );
        
        // Try to find Head and Tail
        
        // Point coordinates for head/tail
        cv::Point Head,Tail;
        
        // Map to keep the distances of each point to all the others
        // Pair of points to keep the points with the Maximum distance (i.e. head and tail :) )
        std::vector<cv::Point> cntPoints;
        blobToPointVector(blob,cntPoints);
        LarvaDistanceMap Distances(cntPoints);

        // Compute all the inner distances for the larva
        computeInnerDistances(blob,Distances,newLarvaSkel.MidPoint);

        cur_larva.length_mean=((cur_larva.length_mean+Distances.MaxDist)/2);
      
        if (cur_larva.area_max < Distances.MaxDist)
        {
          cur_larva.area_max=Distances.MaxDist;
        }
        if (cur_larva.area_min > Distances.MaxDist)
        {
          cur_larva.area_min=Distances.MaxDist;
        }
        PointPair MAXPair=Distances.MaxDistPoints;

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
        cur_larva.inBlob.push_back(false);
      }
      else
      {
        // Larva was found in a blob
        
        //Keep the larva blob for the compound shape
        detected_larvae[detected_clusters[ID][1]].blobs.push_back(blob);
        // Set head/tail to a null value
        cur_larva.heads.push_back(cvPoint(0,0));
        cur_larva.tails.push_back(cvPoint(0,0));
        cur_larva.heads.push_back(cvPoint(0,0));
        cur_larva.tails.push_back(cvPoint(0,0));
        cur_larva.centroid_speed_x.push_back(0.0);
        cur_larva.centroid_speed_y.push_back(0.0);
        larvaSkel empty(true); // Create an empty skeleton
        cur_larva.lrvskels.push_back(empty);
        //
        // Update the current larva blob history to true
        cur_larva.inBlob.push_back(true);
      }
      cur_larva.grey_value_mean = 0; //TODO
      cur_larva.grey_value_max = 0; //TODO
      cur_larva.grey_value_min = 0; //TODO
    }
    it++;
  }
}


// Matching of diverging larvae
// TODO: Generalize for more than two
void diverge_match(unsigned int &candidate_larva_a, 
                   unsigned int &candidate_larva_b,
                   cvb::CvBlob  *newLarva1, 
                   cvb::CvBlob *newLarva2)
{
  
  cv::Point mdpLarva1, mdpLarva2;
  std::vector<cv::Point> newLarva1Points;
  std::vector<cv::Point> newLarva2Points;
  blobToPointVector(*newLarva1,newLarva1Points);
  blobToPointVector(*newLarva2,newLarva2Points);
  LarvaDistanceMap dstLarva1(newLarva1Points), dstLarva2(newLarva2Points);

  int i=detected_larvae[candidate_larva_a].inBlob.size()-1;
  cv::Point centroid_a;
  cv::Point centroid_b;
  cv::Point centroid_blob = detected_larvae[candidate_larva_a].centroids.back();
  int frameCount=0;
  for (;i>=0;i--)
  {
    ++frameCount;
    int ibreak=0;
    //look for last frame where we were out of the blob
    if(detected_larvae[candidate_larva_a].inBlob[i]==false)
    {
      centroid_a=detected_larvae[candidate_larva_a].centroids[i];
      ibreak++;
    }
    if(detected_larvae[candidate_larva_b].inBlob[i]==false)
    {
      centroid_b=detected_larvae[candidate_larva_b].centroids[i];
      ibreak++;
    }
    if (ibreak==2)
      ibreak;
  }

//  double dist_a=(centroid_blob.x-centroid_a.x) + (centroid_blob.y - centroid_a.y);
 // double dist_b=(centroid_blob.x-centroid_b.x) + (centroid_blob.y - centroid_b.y);

//  double dist_1=(centroid_1.x-centroid_blob.x) + (centroid_1.y - centroid_blob.y);
//  double dist_2=(centroid_2.x-centroid_blob.x) + (centroid_2.y - centroid_blob.y);

  computeInnerDistances(*newLarva1,dstLarva1,mdpLarva1);
  computeInnerDistances(*newLarva2,dstLarva2,mdpLarva2);

  double size_a=detected_larvae[candidate_larva_a].area_mean;
  double size_b=detected_larvae[candidate_larva_b].area_mean;
  double size_1=newLarva1->area;
  double size_2=newLarva2->area;

  double length_a=detected_larvae[candidate_larva_a].length_max;
  double length_b=detected_larvae[candidate_larva_b].length_max;
  double length_1=dstLarva1.MaxDist;
  double length_2=dstLarva2.MaxDist;
  
  double lrvhash_a = LARVAE_HASH(length_a,size_a);
  double lrvhash_b = LARVAE_HASH(length_b,size_b);
  double lrvhash_1 = LARVAE_HASH(length_1,size_1);
  double lrvhash_2 = LARVAE_HASH(length_2,size_2);

  if (fabs(lrvhash_a-lrvhash_1)+fabs(lrvhash_b-lrvhash_2) <
      fabs(lrvhash_a-lrvhash_2)+fabs(lrvhash_b-lrvhash_1))
  {
    candidate_larva_a=newLarva1->label;
    candidate_larva_b=newLarva2->label;
  }
  else
  {
    candidate_larva_a=newLarva2->label;
    candidate_larva_b=newLarva1->label;
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
      cvb::CvBlob preBlob=*((*prevIt).second);
     
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
          cvb::CvBlob blob;
          blob=*((*it).second);
          // MIN_DISCARD_DISTANCE is used to quickly exclude larvae that are
          // too far to be considered
          if (((XVAL=fabs(blob.centroid.x - preBlob.centroid.x)) < MIN_DISCARD_DISTANCE ) &&
              ((YVAL=fabs(blob.centroid.y - preBlob.centroid.y)) < MIN_DISCARD_DISTANCE ))
            {
              // we only use the sum of the differences (sqrt is expensive) which in this case should be sufficient
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
            cvb::CvBlob blob;
            blob=*(In[minLabel]);
            // One or more larvae have matched to the new blob with ID minLabel
            // We will now test whether the centroids of the other matches give a
            //   good aproximation of the new centroid.
            //   TODO: Take care of case of more than 2 merging at one point.
            double newx = ((*prevIt).second->centroid.x + Prev[used_map[minLabel][0]]->centroid.x);
            double newy = ((*prevIt).second->centroid.y + Prev[used_map[minLabel][0]]->centroid.y);
            double XDIF,YDIF;
            // TODO: Use correct centroid calculation based on area not just middle
            if (((XDIF=fabs(2*blob.centroid.x - newx)) < 10 ) &&
                ((YDIF=fabs(2*blob.centroid.y - newy)) < 10))
            {
              // We're in luck! The center's match :) We have a cluster.
              used_map[minLabel].push_back((*prevIt).first);
	      // Add both larvae to the detected clusters
              detected_clusters[blob.label].push_back(used_map[minLabel][0]);
              detected_clusters[blob.label].push_back(used_map[minLabel][1]);
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
          // If the detected_clusters map contains a cluster for the 
          // previous blob with the one we are examining we might be 
          // looking into larvae that have separated.
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
                cvb::CvBlob blob;
                blob=*((*it).second);
                if (((XVAL=fabs(blob.centroid.x - preBlob.centroid.x)) < 90) &&
                    ((YVAL=fabs(blob.centroid.y - preBlob.centroid.y)) < 90))
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
            else
            {
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
  //cv::VideoCapture capture("/Users/epaisios/Desktop/LarvaeCapture201302211115.mp4");
  //cv::VideoCapture capture("/Users/epaisios/Downloads/lc2-processed.mp4");
  cv::VideoCapture capture("/Users/epaisios/Downloads/lc3-processed.mp4");
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
  cv::Mat grey_bgFrame;
  unsigned int fg_threshold=10;
  unsigned int thresholdlow=30;
  unsigned int thresholdhigh=255;
  LARVAE_COUNT=0;

  if (!capture.isOpened())
    return 1;

  bool stop(false);
  cv::namedWindow("Extracted Frame");

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

      // Captured frame to BW (necessary for various filters and speed)
      cvtColor(frame,grey_frame,CV_BGR2GRAY);
      
      //Background subtraction is done by eliminating areas outside of the
      //detected petri-dish
      // Here we subtract from the BW captured frame the background frame.
      cv::addWeighted(grey_frame, 1.0, grey_bgFrame, -3.0, 0.0, fg_frame);

      fgROI=cv::Mat::zeros(fg_frame.rows , fg_frame.cols,fg_frame.depth());
      
      cv::circle(fgROI, cv::Point(circles[0][0],circles[0][1]),int(circles[0][2]/1.0),cv::Scalar(255),-1);
      cv::circle(frame, cv::Point(circles[0][0],circles[0][1]),int(circles[0][2]/1.0),cv::Scalar(0,255,0),1);
      
      fg_image=fg_frame&fgROI;
      
      cv::Mat fg_image_norm;
      
      cv::normalize(fg_image,fg_image_norm,0,255,cv::NORM_MINMAX);
      
      cv::threshold(fg_image_norm,
          thresholded_frame,
          thresholdlow,
          thresholdhigh,
          cv::THRESH_BINARY);

      cvb::CvBlobs blobs;
      IplImage ipl_thresholded = thresholded_frame;
      labelImg=cvCreateImage(cvGetSize(&ipl_thresholded), IPL_DEPTH_LABEL, 1);

      unsigned int result=cvLabel(&ipl_thresholded, labelImg, blobs);
      cvb::cvFilterByArea(blobs, 31, 400);

      cvb::CvBlobs tracked_blobs;
      cvb::CvBlobs blob1;
      cvb::CvBlobs::iterator it=blobs.begin();

      if(preBlobs.size()>0)
      {
        larvae_track(blobs,preBlobs,tracked_blobs);
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
        LARVAE_COUNT=preBlobs.size();
      }

      if (SHOWTAGS!=0)
      {
        it=tracked_blobs.begin();
        while (it!=tracked_blobs.end())
        {
          std::stringstream sstm;
          sstm << (*it).first;
          cvb::CvBlob blob=*((*it).second);
          cv::putText(frame,
              sstm.str(),
              cv::Point2d(blob.centroid.x+12,blob.centroid.y+12),
              cv::FONT_HERSHEY_PLAIN,
              0.8,
              cv::Scalar(255,255,255),
              1,
              CV_AA);
        
          cv::Mat larvaROI(frame, cv::Rect(blob.minx-ROI_PADDING,blob.miny-ROI_PADDING,blob.maxx-blob.minx+2*ROI_PADDING,blob.maxy-blob.miny+2*ROI_PADDING));
          detected_larvae[it->first].lrvskels.back().drawSkeleton(larvaROI,cv::Scalar(0,0,255));
          cv::circle(frame,
              cv::Point2d(blob.centroid.x,blob.centroid.y), // circle centre
              1,       // circle radius
              cv::Scalar(255,0,0), // color
              -1);              // thickness
          /*
          detected_larvae[it->first].lrvskels.back().drawSkeleton(larvaROI,cv::Scalar(0,0,255));
          cv::circle(frame,
              cv::Point2d(blob.contour.startingPoint.x,blob.contour.startingPoint.y), // circle centre
              1,       // circle radius
              cv::Scalar(0,255,0), // color
              -1);              // thickness
              */
          cv::circle(larvaROI,
              cv::Point2d(detected_larvae[it->first].lrvskels.back().MidPoint.x,detected_larvae[it->first].lrvskels.back().MidPoint.y), // circle centre
              1,       // circle radius
              cv::Scalar(255,255,0), // color
              -1);              // thickness
          cv::circle(larvaROI,
              cv::Point2d(detected_larvae[it->first].heads.back().x,detected_larvae[it->first].heads.back().y), // circle centre
              1,       // circle radius
              cv::Scalar(0,255,0), // color
              -1);              // thickness
          cv::circle(larvaROI,
              cv::Point2d(detected_larvae[it->first].tails.back().x,detected_larvae[it->first].tails.back().y), // circle centre
              1,       // circle radius
              cv::Scalar(0,0,255), // color
              -1);              // thickness
          
          it++;
        }
      }

      cvReleaseImage(&labelImg);

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
