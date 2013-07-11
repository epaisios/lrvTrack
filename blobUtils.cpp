#include "blobUtils.hpp"
#include "lrvTrackBase.hpp"

double diff(cv::Point &a, cv::Point &b)
{
  return (fabs((double) a.x-b.x)+fabs((double)a.y-b.y));
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
                        int type)
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
  sizes[0]=static_cast<int> (cntPoly->size());
  for (unsigned int i=0; i<cntPoly->size(); ++i)
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
                              int type)
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

  sizes[0]=static_cast<int> (cntPoly->size());
  for (unsigned int i=0; i<cntPoly->size(); ++i)
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

double angle( cv::Point &pt1, cv::Point &pt0, cv::Point &pt2 )
{
  double dx1 = pt1.x - pt0.x;
  double dy1 = pt1.y - pt0.y;
  double dx2 = pt2.x - pt0.x;
  double dy2 = pt2.y - pt0.y;
  return acos((dx1*dx2 + dy1*dy2)/sqrt((dx1*dx1 + dy1*dy1)*(dx2*dx2 + dy2*dy2) + 1e-10));
}

double plotAngle(cvb::CvBlob *blob,cv::Mat &ROIimg,int PAD)
{
  double angle = cvb::cvAngle(blob);

  double x1,y1,x2,y2;
  double cx,cy;
  double lengthLine = MAX(blob->maxx-blob->minx, blob->maxy-blob->miny)/2.;

  cx=blob->centroid.x-blob->minx+PAD;
  cy=blob->centroid.y-blob->miny+PAD;

  x1=cx-lengthLine*cos(angle);
  y1=cy-lengthLine*sin(angle);
  x2=cx+lengthLine*cos(angle);
  y2=cy+lengthLine*sin(angle);
  cv::line(ROIimg,
           cv::Point(int(x1),int(y1)),
           cv::Point(int(x2),int(y2)),
           cv::Scalar(0,255,0));
  return angle;
}

double getGreyValue(cv::Mat &larvaROI, cvb::CvBlob &blob,cv::Mat &grey_frame)
{
  cv::Mat ROI;
  //TODO: Fix when the Padding exceeds the image size!!!
  cv::Mat ROIcopy=grey_frame(cv::Rect(blob.minx-ROI_PADDING,
                                      blob.miny-ROI_PADDING,
                                      blob.maxx-blob.minx+2*ROI_PADDING,
                                      blob.maxy-blob.miny+2*ROI_PADDING)
                            );
  ROIcopy.copyTo(ROI);
  ROI=ROI&larvaROI;
  lrvTrackNormalize(ROI, ROI, 0, 255, CV_MINMAX );
  double nz=cv::norm(ROI,cv::NORM_L1);
  return nz/blob.area;
}

double getPerimeter(cvb::CvBlob &blob)
{
  std::vector<cv::Point> cntPoints;
  blobToPointVector(blob,cntPoints);
  return arcLength(cntPoints, true);
}

double getSurroundingSize(cv::Point &point, cvb::CvBlob &blob, cv::Mat &grey_frame)
{
  cv::Mat larvaImg,lrvROI;
  grey_frame.copyTo(larvaImg);

  cv::Mat ROI=larvaImg(cv::Rect(blob.minx-ROI_PADDING,
                                blob.miny-ROI_PADDING,
                                blob.maxx-blob.minx+2*ROI_PADDING,
                                blob.maxy-blob.miny+2*ROI_PADDING)
                      );
  createLarvaContour(lrvROI, blob);
  ROI=ROI&lrvROI;
  cv::Mat cROI(ROI.size(),ROI.depth());
  cROI=cv::Scalar(0);
  cv::circle(cROI, cv::Point(point.x,point.y),3,cv::Scalar(255),-1);
  cv::Mat area=ROI&cROI;
  double nz=cv::norm(area,cv::NORM_L1);
  //double nz=cv::countNonZero(area);
  return nz;
}
