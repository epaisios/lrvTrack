#include "alglib/stdafx.h"
#include "blobUtils.hpp"
#include "gnuplot_i.hpp"
#include "lrvTrackBase.hpp"
#include "boost/dynamic_bitset.hpp"
#include "alglib/interpolation.h"

//cv::Point2f px3Smooth(cv::Mat &f, cv::Point2f &a, cv::Point2f &b, cv::Point2f &c)
cv::Point2f px3Smooth(cv::Mat &f, cv::Point2f &e, cv::Point2f &a, cv::Point2f &b, cv::Point2f &c, cv::Point2f &d)
{
  cv::Vec3b va=f.at<cv::Vec3b>(a);
  cv::Vec3b vb=f.at<cv::Vec3b>(b);
  cv::Vec3b vc=f.at<cv::Vec3b>(c);
  
  float ua=(va[0]+va[1]+va[2])/3;
  float ub=(vb[0]+vb[1]+vb[2])/3;
  float uc=(vc[0]+vc[1]+vc[2])/3;

  float sum=ua+ub+uc;
  float wa=ua/sum;
  float wb=ub/sum;
  float wc=uc/sum;

  return wa*a+wb*b+wc*c;
}

void derivVec(std::vector<cv::Point2f> &in,
              std::vector<float> &p,
    std::vector<cv::Point2f> &d)
{
  d.push_back(cv::Point(0,0));
  for(unsigned int i=1;i<in.size();i++)
  {
    double R=1/(p[i]-p[i-1]);
    d.push_back((in[i]-in[i-1])*R);
  }
  d[0]=(in[0]-in.back())*(1/(p[p.size()-1]-p[p.size()-2]));
}

void deriv2Vec(
    std::vector<cv::Point2f> &d1,
              std::vector<float> &p,
    std::vector<cv::Point2f> &d2)
{
  d2.push_back(cv::Point(0,0));
  for(unsigned int i=1;i<d1.size();i++)
  {
    double R=1/(p[i]-p[i-1]);
    d2.push_back((d1[i]-d1[i-1])*R);
  }
  d2[0]=(d1[0]-d1.back())*(1/(p[p.size()-1]-p[p.size()-2]));
}

void curvVec(std::vector<cv::Point2f> &in,
             std::vector<float> &p,
             std::vector<float> &c)
{

  std::vector<cv::Point2f> d1,d2;
  derivVec(in,p,d1);
  deriv2Vec(d1,p,d2);
  for(unsigned int i=0;i<in.size();i++)
  {
    c.push_back(
        (d1[i].x*d2[i].y-d1[i].y*d2[i].x)/
        pow(d1[i].x*d1[i].x+d1[i].y*d1[i].y,1.5));
  }
}

cv::Point2f px5Smooth(cv::Mat &f, cv::Point2f &a, cv::Point2f &b, cv::Point2f &c, cv::Point2f &d, cv::Point2f &e)
{
  cv::Vec3b va=f.at<cv::Vec3b>(a);
  cv::Vec3b vb=f.at<cv::Vec3b>(b);
  cv::Vec3b vc=f.at<cv::Vec3b>(c);
  cv::Vec3b vd=f.at<cv::Vec3b>(d);
  cv::Vec3b ve=f.at<cv::Vec3b>(e);
  
  float ua=(va[0]+va[1]+va[2])/3;
  float ub=(vb[0]+vb[1]+vb[2])/3;
  float uc=(vc[0]+vc[1]+vc[2])/3;
  float ud=(vd[0]+vd[1]+vd[2])/3;
  float ue=(ve[0]+ve[1]+ve[2])/3;

  float sum=ua+ub+uc+ue+ud;
  float wa=ua/sum;
  float wb=ub/sum;
  float wc=uc/sum;
  float wd=ud/sum;
  float we=ue/sum;

  return wa*a+wb*b+wc*c+wd*d+we*e;
}


void smoothPoints(std::vector<cv::Point2f> &origPoints, std::vector<cv::Point2f> &smoothened,cv::Mat &frame, unsigned int rep)
{
  std::vector<cv::Point2f> temp1=origPoints;
  for (unsigned int k=0;k<rep;k++)
  {
    cv::Point2f NP=px3Smooth(frame, temp1[temp1.size()-2], temp1.back(), temp1[0], temp1[1], temp1[2]); 
    smoothened.push_back(NP);
    for(unsigned int i=1;i<temp1.size();i++)
    {
      if(temp1[i-1]!=temp1[i])
      {
        cv::Point2f NP;
        if(i==1)
          NP=px3Smooth(frame,temp1.back(),temp1[i-1],temp1[i],temp1[i+1],temp1[i+2]);
        else if(i<temp1.size()-2)
          NP=px3Smooth(frame,temp1[i-2],temp1[i-1],temp1[i],temp1[i+1],temp1[i+2]);
        else if (i==temp1.size()-1)
          NP=px3Smooth(frame,temp1[i-2],temp1[i-1],temp1[i],temp1[0],temp1[1]);
        else if (i==temp1.size()-2)
          NP=px3Smooth(frame,temp1[i-2],temp1[i-1],temp1[i],temp1[i+1],temp1[0]);

        if(NP!=smoothened.back())
        {
          smoothened.push_back(NP);
        }
      }
    }
    if(k!=rep-1)
    {
      temp1=smoothened;
      smoothened.clear();
    }
  }
}


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

void blobToContourVector(cvb::CvBlob &blob,
                         cv::Mat &frame, 
                         int PAD,
                         std::vector<cv::Point2f> &points)
{
  std::vector<cv::Point2f> kpoints;
  blobToPointVector(blob,kpoints);

  cv::Mat ROI=cv::Mat(frame,
              cv::Rect(blob.minx-PAD,
                       blob.miny-PAD,
                       blob.maxx-blob.minx+1+2*PAD,
                       blob.maxy-blob.miny+1+2*PAD
                      )
             );

  cv::Point2f PADpoint(blob.minx-PAD,blob.miny-PAD);
  for(unsigned int i=1;i<kpoints.size();i++)
  {
    cv::LineIterator it(ROI,kpoints[i-1]-PADpoint,kpoints[i]-PADpoint,8);
    for(int j = 1; j < it.count; j++, ++it)
    {
      cv::Point2f newPoint;
      newPoint.x=it.pos().x;
      newPoint.y=it.pos().y;
      points.push_back(newPoint+PADpoint);
    }
  }

  cv::LineIterator it(ROI,kpoints.back()-PADpoint,kpoints[0]-PADpoint,8);
  for(int i = 0; i < it.count; i++, ++it)
  {
    cv::Point2f newPoint;
    newPoint.x=it.pos().x;
    newPoint.y=it.pos().y;
    points.push_back(newPoint+PADpoint);
  }
}

void lengthAreaPerimeter(double a,double b,double &length, double &width)
{
double x1=-0.0459441*sqrt(a*a+25.1327*b)+0.0649747*sqrt(-(43.5312*a*b)/sqrt(a*a+25.1327*b)+2*a*a-(1.73205*a*a*a)/sqrt(a*a+25.1327*b)-62.8319*b)+0.0795775*a;
double y1=(0.0530516*(0.028715*a*a*sqrt(a*a+25.1327*b)-(0.866025*a*a*b)/sqrt(a*a+25.1327*b)+0.00574301*pow((a*a+25.1327*b),1.5)-1.01036*b*sqrt(a*a+25.1327*b)-(0.0344581*a*a*a*a)/sqrt(a*a+25.1327*b)+0.0324874*a*a*sqrt(-(43.5312*a*b)/sqrt(a*a+25.1327*b)+2*a*a-(1.73205*a*a*a)/sqrt(a*a+25.1327*b)-62.8319*b)-0.0281349*a*sqrt(a*a+25.1327*b)*sqrt(-(43.5312*a*b)/sqrt(a*a+25.1327*b)+2*a*a-(1.73205*a*a*a)/sqrt(a*a+25.1327*b)-62.8319*b)-0.0162437*pow((-(43.5312*a*b)/sqrt(a*a+25.1327*b)+2*a*a-(1.73205*a*a*a)/sqrt(a*a+25.1327*b)-62.8319*b),1.5)-2.24537*b*sqrt(-(43.5312*a*b)/sqrt(a*a+25.1327*b)+2*a*a-(1.73205*a*a*a)/sqrt(a*a+25.1327*b)-62.8319*b)+1.5*a*b))/b;

double x2=-0.0459441*sqrt(a*a+25.1327*b)-0.0649747*sqrt(-(43.5312*a*b)/sqrt(a*a+25.1327*b)+2*a*a-(1.73205*a*a*a)/sqrt(a*a+25.1327*b)-62.8319*b)+0.0795775*a;
double y2=(0.0530516*(0.028715*a*a*sqrt(a*a+25.1327*b)-(0.866025*a*a*b)/sqrt(a*a+25.1327*b)+0.00574301*pow((a*a+25.1327*b),1.5)-1.01036*b*sqrt(a*a+25.1327*b)-(0.0344581*a*a*a*a)/sqrt(a*a+25.1327*b)-0.0324874*a*a*sqrt(-(43.5312*a*b)/sqrt(a*a+25.1327*b)+2*a*a-(1.73205*a*a*a)/sqrt(a*a+25.1327*b)-62.8319*b)+0.0281349*a*sqrt(a*a+25.1327*b)*sqrt(-(43.5312*a*b)/sqrt(a*a+25.1327*b)+2*a*a-(1.73205*a*a*a)/sqrt(a*a+25.1327*b)-62.8319*b)+0.0162437*pow((-(43.5312*a*b)/sqrt(a*a+25.1327*b)+2*a*a-(1.73205*a*a*a)/sqrt(a*a+25.1327*b)-62.8319*b),1.5)+2.24537*b*sqrt(-(43.5312*a*b)/sqrt(a*a+25.1327*b)+2*a*a-(1.73205*a*a*a)/sqrt(a*a+25.1327*b)-62.8319*b)+1.5*a*b))/b;

double x3=0.0459441*sqrt(a*a+25.1327*b)+0.0649747*sqrt((43.5312*a*b)/sqrt(a*a+25.1327*b)+2*a*a+(1.73205*a*a*a)/sqrt(a*a+25.1327*b)-62.8319*b)+0.0795775*a;
double y3=(0.0530516*(-0.028715*a*a*sqrt(a*a+25.1327*b)+(0.866025*a*a*b)/sqrt(a*a+25.1327*b)-0.00574301*pow((a*a+25.1327*b),1.5)+1.01036*b*sqrt(a*a+25.1327*b)+(0.0344581*a*a*a*a)/sqrt(a*a+25.1327*b)+0.0324874*a*a*sqrt((43.5312*a*b)/sqrt(a*a+25.1327*b)+2*a*a+(1.73205*a*a*a)/sqrt(a*a+25.1327*b)-62.8319*b)+0.0281349*a*sqrt(a*a+25.1327*b)*sqrt((43.5312*a*b)/sqrt(a*a+25.1327*b)+2*a*a+(1.73205*a*a*a)/sqrt(a*a+25.1327*b)-62.8319*b)-0.0162437*pow(((43.5312*a*b)/sqrt(a*a+25.1327*b)+2*a*a+(1.73205*a*a*a)/sqrt(a*a+25.1327*b)-62.8319*b),1.5)-2.24537*b*sqrt((43.5312*a*b)/sqrt(a*a+25.1327*b)+2*a*a+(1.73205*a*a*a)/sqrt(a*a+25.1327*b)-62.8319*b)+1.5*a*b))/b;

double x4=0.0459441*sqrt(a*a+25.1327*b)-0.0649747*sqrt((43.5312*a*b)/sqrt(a*a+25.1327*b)+2*a*a+(1.73205*a*a*a)/sqrt(a*a+25.1327*b)-62.8319*b)+0.0795775*a;
double y4=(0.0530516*(-0.028715*a*a*sqrt(a*a+25.1327*b)+(0.866025*a*a*b)/sqrt(a*a+25.1327*b)-0.00574301*pow((a*a+25.1327*b),1.5)+1.01036*b*sqrt(a*a+25.1327*b)+(0.0344581*a*a*a*a)/sqrt(a*a+25.1327*b)-0.0324874*a*a*sqrt((43.5312*a*b)/sqrt(a*a+25.1327*b)+2*a*a+(1.73205*a*a*a)/sqrt(a*a+25.1327*b)-62.8319*b)-0.0281349*a*sqrt(a*a+25.1327*b)*sqrt((43.5312*a*b)/sqrt(a*a+25.1327*b)+2*a*a+(1.73205*a*a*a)/sqrt(a*a+25.1327*b)-62.8319*b)+0.0162437*pow(((43.5312*a*b)/sqrt(a*a+25.1327*b)+2*a*a+(1.73205*a*a*a)/sqrt(a*a+25.1327*b)-62.8319*b),1.5)+2.24537*b*sqrt((43.5312*a*b)/sqrt(a*a+25.1327*b)+2*a*a+(1.73205*a*a*a)/sqrt(a*a+25.1327*b)-62.8319*b)+1.5*a*b))/b;

if(x1>0 && y1>0)
{
  length = std::max(x1,y1);
  width = std::min(x1,y1);
}
if(x2>0 && y2>0)
{
  length = std::max(x2,y2);
  width = std::min(x2,y2);
}
if(x3>0 && y3>0)
{
  length = std::max(x3,y3);
  width = std::min(x3,y3);
}
if(x4>0 && y4>0)
{
  length = 2 * std::max(x4,y4);
  width = 2 * std::min(x4,y4);
}

}

void createLarvaROI(cv::Mat &frame, cv::Mat &ROI, cvb::CvBlob &blob)
{
  ROI=cv::Mat(frame,
              cv::Rect(blob.minx-ROI_PADDING,
                       blob.miny-ROI_PADDING,
                       blob.maxx-blob.minx+1+2*ROI_PADDING,
                       blob.maxy-blob.miny+1+2*ROI_PADDING
                      )
             );
  //cv::normalize(ROI,ROI,0,255,cv::NORM_MINMAX);
}

/* BROKEN!! IF FIXED MAY BE FASTER
void createLarvaContourCV(cv::Mat &lrvROI,
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

  std::vector<cv::Point2f> cntVec;

  for (unsigned int i=0; i<cntPoly->size(); ++i)
  {
    ContourPoints[0][i].x=(*cntPoly)[i].x-blob.minx+ROI_PADDING+PADDING;
    ContourPoints[0][i].y=(*cntPoly)[i].y-blob.miny+ROI_PADDING+PADDING;
    cntVec.push_back(ContourPoints[0][i]);
  }

  sizes[0]=static_cast<int> (cntPoly->size());
  std::vector<std::vector<cv::Point2f> > cvec;
  cvec.push_back(cntVec);

  cv::Scalar color;
  if(type==CV_8UC1)
    color=cv::Scalar(255);
  else
    color=cv::Scalar(255,255,255);

  drawContours(lrvROI,cvec,0,color,1,4,cv::noArray(),INT_MAX,cv::Point(PADDING,PADDING));
}
*/

void createLarvaContour(cv::Mat &lrvROI,
                        cvb::CvBlob &blob,
                        int type,
                        int PADDING,
                        bool FILL,
                        cv::Scalar color,
                        int connectivity)
{
  int sizes[1];
  cvb::CvContourPolygon *cntPoly=
    cvb::cvConvertChainCodesToPolygon(&blob.contour);

  if (lrvROI.empty())
    lrvROI=cv::Mat::zeros(blob.maxy-blob.miny+1+(2*ROI_PADDING)+(2*PADDING),
                        blob.maxx-blob.minx+1+(2*ROI_PADDING)+(2*PADDING),
                        type);

  cv::Point *ContourPoints[1];
  ContourPoints[0]=(cv::Point*) malloc(
                     blob.contour.chainCode.size()*sizeof(cv::Point)
                   );
  sizes[0]=static_cast<int> (cntPoly->size());

  for (unsigned int i=0; i<cntPoly->size(); ++i)
  {
    ContourPoints[0][i].x=(*cntPoly)[i].x-blob.minx+ROI_PADDING+PADDING;
    ContourPoints[0][i].y=(*cntPoly)[i].y-blob.miny+ROI_PADDING+PADDING;
    if(!FILL && i>0)
    {
      cv::line(lrvROI,ContourPoints[0][i-1],ContourPoints[0][i],color,1,connectivity);
    }
  }
  if(!FILL)
    cv::line(lrvROI,ContourPoints[0][cntPoly->size()-1],ContourPoints[0][0],color,1,connectivity);


  if(FILL)
  {
    cv::fillPoly(lrvROI,
        (const cv::Point**) ContourPoints,
        sizes,
        1,
        color,
        connectivity
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

void createLarvaContourPacked(cv::Point &first, 
                              unsigned int &size,
                              std::string &STR,
                              cvb::CvBlob &blob)
{
  std::vector<cv::Point2f> contourPoints;
  blobToPointVector(blob, contourPoints);
  first.x=(int) contourPoints.back().x;
  first.y=(int) contourPoints.back().y;
  std::stringstream outline;
  unsigned int acc=0;
  unsigned int cntSize=0;
  cv::Mat lrvROI;
  createLarvaContour(lrvROI,blob,CV_8UC1,0,false);
  int cnt3=0;

  std::vector<cv::Point2f>::reverse_iterator P=contourPoints.rbegin()+1;
  std::vector<cv::Point2f>::reverse_iterator pP=contourPoints.rbegin();
  for(;pP!=contourPoints.rend();)
  {
    cv::Point a;
    a.x= (int) P->x - blob.minx;
    a.y= (int) P->y - blob.miny;
    cv::Point pre;
    pre.x= (int) pP->x - blob.minx;
    pre.y= (int) pP->y - blob.miny;

    cv::LineIterator it(lrvROI, pre, a, 4);
    ++it;
    std::vector<uchar> buf(it.count);
    cv::Point cur;
    for(int i = 1; i < it.count; i++, ++it)
    {
      acc <<= 2;
      //contourSize++;
      cur.x=(int) it.pos().x;
      cur.y=(int) it.pos().y;
      cv::Point d=cur-pre;
      if(d.x==1 && d.y==0)
      {
        acc |= 0x1;
      }
      else if(d.y==1 && d.x==0)
      {
        acc |= 0x3;
      }
      else if(d.x==-1 && d.y==0)
      {
        acc |= 0x0;
      }
      else if(d.y==-1 && d.x==0)
      {
        acc |= 0x2;
      }
      else
      {
        std::cerr << "Error creating 4-connected bitset!!!" << std::endl;
      }
      cnt3++;
      cntSize++;
      lrvROI.at<uchar>(cur)=100;
      lrvROI.at<uchar>(a)=50;
      if (cnt3==3)
      {
        cnt3=0;
        outline << (uchar) (acc+48);
        acc=0;
      }
      pre.x=cur.x;
      pre.y=cur.y;
    }
    ++pP;
    ++P;
    if(P==contourPoints.rend())
      P=contourPoints.rbegin();
  }

  if (cnt3==1)
  {
    acc <<= 2;
    acc <<= 2;
    acc |= 0x1;
    cnt3=0;
    outline << (char) (acc+48);
  }
  if (cnt3==2)
  {
    acc <<= 2;
    cnt3=0;
    outline << (char) (acc+48);
  }
  STR=outline.str();
  size=cntSize;
}


double angle( cv::Point2f &pt1, cv::Point2f &pt0, cv::Point2f &pt2 )
{
  double dx1 = pt1.x - pt0.x;
  double dy1 = pt1.y - pt0.y;
  double dx2 = pt2.x - pt0.x;
  double dy2 = pt2.y - pt0.y;
  return acos((dx1*dx2 + dy1*dy2)/sqrt((dx1*dx1 + dy1*dy1)*(dx2*dx2 + dy2*dy2) + 1e-10));
}
/*
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
*/

double angleD( cv::Point2f &pt1, cv::Point2f &pt0, cv::Point2f &pt2 )
{
  //Absolute Tail/Head Angle
  cv::Point2f newPT1=pt1-pt0;
  double MULVAL=180.0/CV_PI;
  double arc1=atan2(newPT1.y,newPT1.x);
  double arc2=atan2(pt2.y,pt2.x);
  cv::Point2f newPT2=pt2-pt0;
  cv::Point2f rt;
  rt.x=newPT2.x * cos(-arc1) - newPT2.y*sin(-arc1);
  rt.y=newPT2.x * sin(-arc1) + newPT2.y*cos(-arc1);
  arc2=MULVAL*atan2(rt.y,rt.x);
  return fabs(arc2);
}

double angleC( cv::Point2f &pt1, cv::Point2f &pt0, cv::Point2f &pt2 )
{
  //Absolute Tail/Head Angle
  double MULVAL=180.0/CV_PI;
  cv::Point2f newPT1=pt1-pt0;
  cv::Point2f newPT2=pt2-pt0;
  double d1=newPT1.x*newPT1.x+newPT1.y*newPT1.y;
  double d2=newPT2.x*newPT2.x+newPT2.y*newPT2.y;
  double arc1=MULVAL*acos(newPT1.dot(newPT2)/sqrt(d1*d2));
  return arc1;
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
                                      blob.maxx-blob.minx+1+2*ROI_PADDING,
                                      blob.maxy-blob.miny+1+2*ROI_PADDING)
                            );
  ROIcopy.copyTo(ROI);
  ROI=ROI&larvaROI;
  //lrvTrackNormalize(ROI, ROI, 0, 255, CV_MINMAX );
  double nz=cv::norm(ROI,cv::NORM_L1);
  return nz;
}

inline bool baseXsmaller(cv::Point2f a, cv::Point2f b)
{
  return(a.x<b.x);
}

inline bool baseYsmaller(cv::Point2f a, cv::Point2f b)
{
  return(a.y<b.y);
}

void sortPoints(std::vector<cv::Point2f> &cp,
                bool sortX)
{
  if(sortX)
  {
    std::sort(cp.begin(),cp.end(),baseXsmaller);
  }
  else
  {
    std::sort(cp.begin(),cp.end(),baseYsmaller);
  }
}

void wait_for_key ()
{
#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || defined(__TOS_WIN__)  // every keypress registered, also arrow keys
  std::cout << std::endl << "Press any key to continue..." << std::endl;

    FlushConsoleInputBuffer(GetStdHandle(STD_INPUT_HANDLE));
    _getch();
#elif defined(unix) || defined(__unix) || defined(__unix__) || defined(__APPLE__)
    std::cout << std::endl << "Press ENTER to continue..." << std::endl;

    std::cin.clear();
    std::cin.ignore(std::cin.rdbuf()->in_avail());
    std::cin.get();
#endif
    return;
}

void getBestCurvatureS(std::vector<float> &curv,
                      std::vector<unsigned int> &di,
                      std::vector<float> &dmax
                      )
{
  std::vector<float> c;
  int sf=25;
  std::vector<double> idx;
  for(unsigned int i=0;i<curv.size();i++)
  {
    idx.push_back(i);
  }
  //Gnuplot g1("lines");
  //Gnuplot g2("lines");
  //g1.plot_xy(idx,curv,"curvature");
  //wait_for_key();
  smoothVec(curv,c,sf,(float)0.0);

  //g2.plot_xy(idx,c,"curvature smoothened");
  //wait_for_key();

  bool positive=true;
  for(unsigned int i=1;i<c.size();i++)
  {
    if(c[i]-c[i-1]>=0 )
    {
      positive=true;
      continue;
    }
    else if(positive)
    {
      positive=false;
      std::vector<unsigned int>::iterator dj=di.begin();
      for(std::vector<float>::iterator j=dmax.begin();j!=dmax.end();j++)
      {
        if(c[i-1]>(*j))
        {
          dmax.insert(j,c[i-1]);
          di.insert(dj,i);
          dmax.pop_back();
          di.pop_back();
          dj++;
          break;
        }
        dj++;
      }
    }
  }
}

void getBestCurvature(std::vector<float> &c,
                      std::vector<unsigned int> &di,
                      std::vector<float> &dmax
                      )
{
  bool positive=true;
  for(unsigned int i=1;i<c.size();i++)
  {
    if(c[i]-c[i-1]>=0 )
    {
      positive=true;
      continue;
    }
    else if(positive)
    {
      positive=false;
      std::vector<unsigned int>::iterator dj=di.begin();
      for(std::vector<float>::iterator j=dmax.begin();j!=dmax.end();j++)
      {
        if(c[i-1]>(*j))
        {
          dmax.insert(j,c[i-1]);
          di.insert(dj,i);
          dmax.pop_back();
          di.pop_back();
          dj++;
          break;
        }
        dj++;
      }
    }
  }
}

void spline2(std::vector<cv::Point2f> &cp,
           std::vector<float> &d,
           std::vector<float> &w,
           int n,
           int RES,
           /*int d0,
           int dm,*/
           std::vector<cv::Point2f> &np,
           std::vector<unsigned int> &di)
{
  //alglib::real_1d_array x,y;
  alglib::real_2d_array xy;
  alglib::real_1d_array ad,aw;
  alglib::real_1d_array dummy;
  alglib::integer_1d_array dummyint;
  //alglib::ae_int_t info;
  //x.setlength(n);
  //y.setlength(n);
  xy.setlength(n,2);
  ad.setlength(n);
  aw.setlength(n);
  alglib::spline1dinterpolant sx;
  alglib::spline1dinterpolant sy;
  alglib::spline1dfitreport rep;
  //alglib::ae_int_t m=n;
  //double rho=-5.0;
  alglib::pspline2interpolant p;
  std::vector<cv::Point2f> NP;
  for(int i=0;i<n;i++)
  {
    //x[i]=cp[i].x;
    //y[i]=cp[i].y;
    xy[i][0]=cp[i].x;
    xy[i][1]=cp[i].y;
    ad[i]=d[i];
    aw[i]=w[i];
  }

  //alglib::spline1dfitpenalizedw(ad,x,aw,m+2,0.19,info,sx,rep);
  //alglib::spline1dfitpenalizedw(ad,y,aw,m+2,0.19,info,sy,rep);
  //alglib::spline1dfitpenalized(ad,x,m+2,-0.09,info,sx,rep);
  //alglib::spline1dfitpenalized(ad,y,m+2,-0.09,info,sy,rep);
  alglib::pspline2buildperiodic(xy,n,1,2,p);
  
  double t0=0;
  double t;
  double step=(1.0-t0)/RES;
  std::vector<float> dmax(6,0);
  std::vector<float> curvature;
  std::vector<float> pvals;
  for (int i=0;i<RES;i++)
  {
    t=t0+i*step;
    double x,y;
    //double dx,dy;
    //double d2x,d2y;
    alglib::pspline2calc(p,t,x,y);
    cv::Point2f N(x,y);
    NP.push_back(N);
    pvals.push_back(t);
    //alglib::pspline2diff2( p, t, x, dx, d2x, y, dy, d2y);

    //dp.push_back(cv::Point2f(d2x,d2y));
    //double val=(dx*d2y-dy*d2x)/pow(dx*dx+dy*dy,1.5);
    //curvature.push_back(val);
  }
  pvals.push_back(t+step);
  smoothVec(NP,np,13,cv::Point2f(0,0));

  std::vector<float> scurvature;
  curvVec(np,pvals,scurvature);
  getBestCurvatureS(scurvature,di,dmax);
}

void spline(std::vector<cv::Point2f> &cp,
           std::vector<float> &d,
           std::vector<float> &w,
           int n, 
           float xp1, 
           float yp1, 
           float xpn, 
           float ypn, 
           std::vector<float> &x2, 
           std::vector<float> &y2)
{
  int i,k;
  float px,sig,qnx,unx;
  float py,qny,uny;
  
  std::vector<float> yu;
  std::vector<float> xu;
  yu.reserve(n);
  xu.reserve(n);
  x2.reserve(n);
  y2.reserve(n);
  if(xp1 > 0.99e30)
  {
    x2[0]=0.0;
    xu[0]=0.0;
  }
  else
  {
    x2[0]=-0.5;
    xu[0]=(3.0*(d[1]-d[0]))*(w[0]*(cp[1].x-cp[0].x)/(d[1]-d[0])-xp1);
  }

  if(yp1 > 0.99e30)
  {
    y2[0]=0.0;
    yu[0]=0.0;
  }
  else
  {
    y2[0]=-0.5;
    yu[0]=(3.0*(d[1]-d[0]))*(w[0]*(cp[1].y-cp[0].y)/(d[1]-d[0])-yp1);
  }

  for(i=1;i<=n-2;i++)
  {
    sig=(d[i]-d[i-1])/(d[i+1]-d[i-1]);
    px=sig*x2[i-1]+2.0;
    py=sig*y2[i-1]+2.0;
    x2[i]=((sig-1.0)/px);
    y2[i]=((sig-1.0)/py);

    xu[i]=w[i]*(cp[i+1].x-cp[i].x)/(d[i+1]-d[i]) - w[i]*(cp[i].x - cp[i-1].x)/(d[i]-d[i-1]);
    yu[i]=w[i]*(cp[i+1].y-cp[i].y)/(d[i+1]-d[i]) - w[i]*(cp[i].y - cp[i-1].y)/(d[i]-d[i-1]);

    xu[i]=(6.0*xu[i]/(d[i+1]-d[i-1])-sig*xu[i-1])/px;
    yu[i]=(6.0*yu[i]/(d[i+1]-d[i-1])-sig*yu[i-1])/py;
  }
  if (xpn>0.99e30)
    qnx=unx=0.0;
  else
  {
    qnx=0.5;
    unx=(3.0/(d[n-1]-d[n-2]))*(xpn-(w[n-1]*cp[n-1].x-cp[n-2].x)/(d[n-1]-d[n-2]));
  }
  
  if (ypn>0.99e30)
    qny=uny=0.0;
  else
  {
    qny=0.5;
    uny=(3.0/(d[n-1]-d[n-2]))*(ypn-(w[n-1]*cp[n-1].y-cp[n-2].y)/(d[n-1]-d[n-2]));
  }

  x2[n-1]=(unx-qnx*xu[n-2])/(qnx*x2[n-1]+1.0);
  y2[n-1]=(uny-qny*yu[n-2])/(qny*y2[n-1]+1.0);
  for(k=n-2;k>=0;k--)
  {
    x2[k]=x2[k]*x2[k+1]+xu[k];
    y2[k]=y2[k]*y2[k+1]+yu[k];
  }
}

void fin(int n, 
         int m2,
         float h,
         float p,
         std::vector<float> &t,
         std::vector<float> &a,
         std::vector<float> &b,
         std::vector<float> &c,
         std::vector<float> &d,
         std::vector<float> &u,
         std::vector<float> &v,
         std::vector<float> &cp)
{
  for (int i=1;i<=n;i++)
  {
    a[i]=cp[i]-p*v[i];
    c[i]=u[i];
  }
  for (int i=1;i<=m2;i++)
  {
    h=t[i+1]-t[i];
    d[i]=(c[i+1]-c[i])/(3*h);
    b[i]=(a[i+1]-a[i])/h-(h*d[i]+c[i])*h;
  }
}

void splint(std::vector<cv::Point2f> &cp,
            std::vector<float> &x2a, 
            std::vector<float> &y2a, 
            std::vector<float> &d,
            int n, 
            float t, 
            int khi,
            int klo,
            cv::Point2f &np)
{
  float b,a;
  
  double h=d[khi]-d[klo];
  a=(d[khi]-t)/h;
  b=(t-d[klo])/h;
  np.x=a*cp[klo].x+b*cp[khi].x+((a*a*a-a)*x2a[klo]+(b*b*b-b)*x2a[khi])*(h*h)/6.0;
  np.y=a*cp[klo].y+b*cp[khi].y+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/6.0;
}

/*void sspline(std::vector<cv::Point2f> &cp,
           std::vector<float> &t,
           std::vector<float> &w,
           int n, 
           float s,
           newPoints np)

{
  float* xpointsx;
  float* xpointsy;
  xpointsx=(float *)malloc(n*sizeof(float));
  xpointsy=(float *)malloc(n*sizeof(float));

  for(unsigned int i=0;i<xpoints.size();i++)
  {
    xpointsx[i]=xpoints[i].x;
    xpointsy[i]=xpoints[i].y;
  }
  float **v;
  float **a;
  v=(float **)malloc(n*sizeof(float *));
  a=(float **)malloc(n*sizeof(float *));
  for(unsigned int i=0;i<xpoints.size();i++)
  {
    v[i]=(float *)malloc(7*sizeof(float));
    a[i]=(float *)malloc(4*sizeof(float));
  }

  smooth_(&x,&t[0],&w[0],&n,&s,&v,&a);
}*/

double getPerimeter(cvb::CvBlob &blob)
{
  std::vector<cv::Point2f> cntPoints;
  blobToPointVector(blob,cntPoints);
  return arcLength(cntPoints, true);
}

double getSurroundingSize(cv::Point2f &point, cvb::CvBlob &blob, cv::Mat &grey_frame)
{
  cv::Mat larvaImg,lrvROI;
  unsigned int PADDING=4;
  cv::Mat rROI,ROI;
  try{
    rROI=grey_frame(cv::Rect(blob.minx-PADDING,
          blob.miny-PADDING,
          blob.maxx-blob.minx+1+2*PADDING,
          blob.maxy-blob.miny+1+2*PADDING)
        );
    rROI.copyTo(ROI);
    createLarvaContour(lrvROI, blob,CV_8UC1,PADDING);
  }
  catch(...)
  {
    std::cerr << "getSurroundingSize: Error creating ROI" << std::endl;
    return 0;
  }
  cv::Mat element = cv::getStructuringElement(cv::MORPH_CROSS, cv::Size(3, 3));
  cv::dilate(lrvROI,lrvROI,element);
  //lrvTrackNormalize(lrvROI, lrvROI, 0, 255, CV_MINMAX );
  ROI=ROI&lrvROI;
  /*
  cv::Mat dbg;
  ROI.copyTo(dbg);
  cv::circle(dbg,
      cv::Point2d(point.x+PADDING,point.y+PADDING),
      0,
      cv::Scalar(255),
      -1);

  cv::resize(dbg,dbg,cv::Size(),8,8,cv::INTER_NEAREST);
  cv::imshow("headtail",dbg);
  cv::waitKey(1);
*/
  cv::Mat cROI(ROI.size(),ROI.depth());
  cROI=cv::Scalar(0);
  cv::circle(cROI, cv::Point2f(point.x+PADDING,point.y+PADDING),6,cv::Scalar(255),-1);
  cv::Mat area=ROI&cROI;
/*
  cv::resize(area,dbg,cv::Size(),8,8,cv::INTER_NEAREST);
  cv::imshow("headtail",dbg);
  cv::waitKey(1);
  */
  double nz=cv::norm(area,cv::NORM_L1);
  double nc=cv::countNonZero(area);
  return nz/nc;
  //return nz;
}
