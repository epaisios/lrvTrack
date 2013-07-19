#include "blobUtils.hpp"
#include "lrvTrackBase.hpp"

double diff(cv::Point2f &a, cv::Point2f &b)
{
  return (fabs((double) a.x-b.x)+fabs((double)a.y-b.y));
}

void blobToPointVector(cvb::CvBlob &p,std::vector<cv::Point2f> &points)
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
                        int type,
                        int PADDING,
                        bool FILL)
{
  int sizes[1];
  cvb::CvContourPolygon *cntPoly=
    cvb::cvConvertChainCodesToPolygon(&blob.contour);

  lrvROI=cv::Mat::zeros(blob.maxy-blob.miny+(2*ROI_PADDING)+(2*PADDING),
                        blob.maxx-blob.minx+(2*ROI_PADDING)+(2*PADDING),
                        type);

  cv::Point *ContourPoints[1];
  ContourPoints[0]=(cv::Point*) malloc(
                     blob.contour.chainCode.size()*sizeof(cv::Point)
                   );
  sizes[0]=static_cast<int> (cntPoly->size());

  cv::Scalar color;
  if(type==CV_8UC1)
    color=cv::Scalar(255);
  else
    color=cv::Scalar(255,255,255);
  
  for (unsigned int i=0; i<cntPoly->size(); ++i)
  {
    ContourPoints[0][i].x=(*cntPoly)[i].x-blob.minx+ROI_PADDING+PADDING;
    ContourPoints[0][i].y=(*cntPoly)[i].y-blob.miny+ROI_PADDING+PADDING;
    if(!FILL && i>0)
    {
      cv::line(lrvROI,ContourPoints[0][i-1],ContourPoints[0][i],color);
    }
  }


  if(FILL)
  {
    cv::fillPoly(lrvROI,
        (const cv::Point**) ContourPoints,
        sizes,
        1,
        color 
        );
  }
  free(ContourPoints[0]);
  delete(cntPoly);

}

void createLarvaContourPoints(cv::Mat &lrvROI,
                              cvb::CvBlob &blob,
                              int type,
                              int PADDING)
{
  int sizes[1];
  cvb::CvContourPolygon *cntPoly=
    cvb::cvConvertChainCodesToPolygon(&blob.contour);
  lrvROI=cv::Mat::zeros(blob.maxy-blob.miny+(2*ROI_PADDING)+(2*PADDING),
                        blob.maxx-blob.minx+(2*ROI_PADDING)+(2*PADDING),
                        type);
  cv::Point2f *ContourPoints[1];
  ContourPoints[0]=(cv::Point2f*) malloc(
                     blob.contour.chainCode.size()*sizeof(cv::Point2f)
                   );
  cv::Scalar color;
  if(type==CV_8UC1)
    color=cv::Scalar(255);
  else
    color=cv::Scalar(255,255,255);

  sizes[0]=static_cast<int> (cntPoly->size());
  for (unsigned int i=0; i<cntPoly->size(); ++i)
    {
      ContourPoints[0][i].x=(*cntPoly)[i].x-blob.minx+ROI_PADDING+PADDING;
      ContourPoints[0][i].y=(*cntPoly)[i].y-blob.miny+ROI_PADDING+PADDING;
      cv::circle(lrvROI,
                 ContourPoints[0][i], // circle centre
                 0,       // circle radius
                 color, // color
                 -1);              // thickness
    }
  free(ContourPoints[0]);
  delete(cntPoly);
}

double angle( cv::Point2f &pt1, cv::Point2f &pt0, cv::Point2f &pt2 )
{
  double dx1 = pt1.x - pt0.x;
  double dy1 = pt1.y - pt0.y;
  double dx2 = pt2.x - pt0.x;
  double dy2 = pt2.y - pt0.y;
  return acos((dx1*dx2 + dy1*dy2)/sqrt((dx1*dx1 + dy1*dy1)*(dx2*dx2 + dy2*dy2) + 1e-10));
}

double angleD( cv::Point2f &pt1, cv::Point2f &pt0, cv::Point2f &pt2 )
{
    cv::Point2f ab = pt0-pt1;
    cv::Point2f cb = pt0-pt2;

    // dot product  
    float dot = (ab.x * cb.x + ab.y * cb.y);

    // length square of both vectors
    float abSqr = ab.x * ab.x + ab.y * ab.y;
    float cbSqr = cb.x * cb.x + cb.y * cb.y;

    // square of cosine of the needed angle    
    float cosSqr = dot * dot / abSqr / cbSqr;

    // this is a known trigonometric equality:
    // cos(alpha * 2) = [ cos(alpha) ]^2 * 2 - 1
    float cos2 = 2 * cosSqr - 1;

    // Here's the only invocation of the heavy function.
    // It's a good idea to check explicitly if cos2 is within [-1 .. 1] range

    const float pi = 3.141592f;

    float alpha2 =
        (cos2 <= -1) ? pi :
        (cos2 >= 1) ? 0 :
        acosf(cos2);

    float rslt = alpha2 / 2;

    float rs = rslt * 180. / pi;


    // Now revolve the ambiguities.
    // 1. If dot product of two vectors is negative - the angle is definitely
    // above 90 degrees. Still we have no information regarding the sign of the angle.

    // NOTE: This ambiguity is the consequence of our method: calculating the cosine
    // of the double angle. This allows us to get rid of calling sqrt.

    if (dot < 0)
        rs = 180 - rs;

    // 2. Determine the sign. For this we'll use the Determinant of two vectors.

    float det = (ab.x * cb.y - ab.y * cb.y);
    if (det < 0)
        rs = -rs;

    return (int) floor(rs + 0.5);
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
           cv::Point2f(int(x1),int(y1)),
           cv::Point2f(int(x2),int(y2)),
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
  return nz;
}

double getPerimeter(cvb::CvBlob &blob)
{
  std::vector<cv::Point2f> cntPoints;
  blobToPointVector(blob,cntPoints);
  return arcLength(cntPoints, true);
}

double getSurroundingSize(cv::Point2f &point, cvb::CvBlob &blob, cv::Mat &grey_frame)
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
  cv::circle(cROI, cv::Point2f(point.x,point.y),3,cv::Scalar(255),-1);
  cv::Mat area=ROI&cROI;
  double nz=cv::norm(area,cv::NORM_L1);
  //double nz=cv::countNonZero(area);
  return nz;
}
