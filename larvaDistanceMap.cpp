#include "larvaDistanceMap.hpp"
#include <queue>
#include "blobUtils.hpp"

void lBFS(int p1,
          std::vector<cv::Point2f> &Points ,
          larvaDistanceMap &Distances
         )
{
  std::queue<int> Q;
  Q.push(p1);
  double MAX=Distances.MaxDist;
  while (!Q.empty())
    {
      int cur=Q.front();
      Q.pop();
      for(unsigned int i=0; i<Points.size() ; ++i)
        {
          PointPair cur_pi = PointPair(Points[cur],Points[i]);
          PointPair p1_pi= PointPair(Points[p1],Points[i]);
          PointPair p1_cur= PointPair(Points[p1],Points[cur]);
          if (Distances[cur][i]>0)
            {
              if (Distances[p1][i]<0)
                Q.push(i);
              float multres=Distances[cur][i]*Distances[p1][cur];
              float sqrtres[1];
              ltsqrt(sqrtres,&multres);
              double newDst = Distances[cur][i]+Distances[p1][cur] +
                              2*sqrtres[0];
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

double p2fdist(cv::Point2f &a, cv::Point2f &b)
{
  float xdiff=a.x - b.x;
  float ydiff=a.y - b.y;
  float sqDst=xdiff*xdiff + ydiff*ydiff;
  float res;
  ltsqrt(&res,&sqDst);
  return (double) res;
}

void pointWithDist(
                   cv::Point2f &a, 
                   cv::Point2f &b, 
                   cv::Point2f &ret,
                   double Dist)
{
  double abDist=p2fdist(a,b);
  double l=Dist/abDist;
  ret.x=a.x+(l*(b.x-a.x));
  ret.y=a.y+(l*(b.y-a.y));
}

void assignPoint(cv::Point2f &a,cv::Point2f &b)
{
  a.x=b.x;
  a.y=b.y;
}

bool distanceMatch(double dist1, double dist2)
{
  double E=1e-6;
  if(dist1>dist2-E && dist1<dist2+E)
    return true;
  else
    return false;
}

void showPoints(cvb::CvBlob blob,
                std::vector<cv::Point2f > points,
                std::vector<cv::Point2f > curSpine,
                unsigned int i,
                unsigned int fwdIdx,
                unsigned int bwdIdx,
                unsigned int pfwdIdx,
                unsigned int pbwdIdx,
                cv::Point2f bwdPoint,
                cv::Point2f fwdPoint,
                unsigned int PADDING
                )
{
  cv::Mat contour;
  cv::Mat lcontour;

  createLarvaContour(contour,blob,CV_8UC3,PADDING);
  cv::circle(contour,
      cv::Point2d(points[i].x-blob.minx+PADDING,points[i].y-blob.miny+PADDING),
      0,
      cv::Scalar(0,0,255),
      -1);
  cv::circle(contour,
      cv::Point2d(points[bwdIdx].x-blob.minx+PADDING,points[bwdIdx].y-blob.miny+PADDING),
      0,
      cv::Scalar(255,0,0),
      -1);
  cv::circle(contour,
      cv::Point2d(points[pbwdIdx].x-blob.minx+PADDING,points[pbwdIdx].y-blob.miny+PADDING),
      0,
      cv::Scalar(255,0,0),
      -1);
  cv::circle(contour,
      cv::Point2d(points[fwdIdx].x-blob.minx+PADDING,points[fwdIdx].y-blob.miny+PADDING),
      0,
      cv::Scalar(0,250,0),
      -1);
  cv::circle(contour,
      cv::Point2d(points[pfwdIdx].x-blob.minx+PADDING,points[pfwdIdx].y-blob.miny+PADDING),
      0,
      cv::Scalar(0,250,0),
      -1);
  cv::circle(contour,
      cv::Point2d(bwdPoint.x-blob.minx+PADDING,bwdPoint.y-blob.miny+PADDING),
      0,
      cv::Scalar(255,0,255),
      -1);
  cv::circle(contour,
      cv::Point2d(fwdPoint.x-blob.minx+PADDING,fwdPoint.y-blob.miny+PADDING),
      0,
      cv::Scalar(255,0,255),
      -1);
  cv::circle(contour,
      cv::Point2d(curSpine.back().x-blob.minx+PADDING,curSpine.back().y-blob.miny+PADDING),
      0,
      cv::Scalar(0,0,255),
      -1);
  cv::resize(contour,lcontour,cv::Size(),8,8,cv::INTER_NEAREST);
  cv::imshow("Spine3",lcontour);
  cv::waitKey(1);

  std::cerr << "Point[i] : " << points[i] << std::endl;
  std::cerr << "Point[bwdIdx] : " << points[bwdIdx] << std::endl;
  std::cerr << "Point[pbwdIdx] : " << points[pbwdIdx] << std::endl;
  std::cerr << "Point[fwdIdx] : " << points[fwdIdx] << std::endl;
  std::cerr << "Point[pfwdIdx] : " << points[pfwdIdx] << std::endl;
  std::cerr << "bwdPoint : " << bwdPoint << std::endl;
  std::cerr << "fwdPoint : " << fwdPoint << std::endl;
  std::cerr << "curSpine.back() : " << curSpine.back() << std::endl;
  std::cerr << "i: " << i << std::endl;
  std::cerr << "bwdIdx: " << bwdIdx << std::endl;
  std::cerr << "pbwdIdx: " << pbwdIdx << std::endl;
  std::cerr << "fwdIdx: " << fwdIdx << std::endl;
  std::cerr << "pfwdIdx: " << pfwdIdx << std::endl;
  std::cerr << "Dist fwdIdx->pfwdIdx: " << p2fdist(points[fwdIdx],points[pfwdIdx]) << std::endl;

  std::cerr << "Dist bwdIdx->pbwdIdx: " << p2fdist(points[bwdIdx],points[pbwdIdx]) << std::endl;
}

void getLeftArc(unsigned int A,
                unsigned int B,
                std::vector<cv::Point2f> &contour,
                std::vector<cv::Point2f> &leftArc,
                std::vector<double> &leftArcDistanceMap
                )
{
  unsigned int I=A;
  
  while(I!=B)
  {
    if(I==A)
    {
      leftArcDistanceMap.push_back(0);
      leftArc.push_back(contour[A]);
      if(I==0)
        I=contour.size()-1;
      else
        I--;
      continue;
    }
    leftArc.push_back(contour[I]);

    if(I==0)
    {
      double newDst=leftArcDistanceMap.back()+p2fdist(contour[I],contour[I+1]);
      leftArcDistanceMap.push_back(newDst);
      I=contour.size()-1;
    }
    else
    {
      if(I==contour.size()-1)
      {
        double newDst=leftArcDistanceMap.back()+p2fdist(contour[I],contour[0]);
        leftArcDistanceMap.push_back(newDst);
      }
      else
      {
        double newDst=leftArcDistanceMap.back()+p2fdist(contour[I],contour[I+1]);
        leftArcDistanceMap.push_back(newDst);
      }
      --I;
    }
  }
  if(B!=contour.size()-1)
  {
    double newDst=leftArcDistanceMap.back()+p2fdist(contour[B],contour[B+1]);
    leftArcDistanceMap.push_back(newDst);
  }
  else
  {
    double newDst=leftArcDistanceMap.back()+p2fdist(contour[B],contour[0]);
    leftArcDistanceMap.push_back(newDst);
  }
  leftArc.push_back(contour[B]);
}

void getRightArc(unsigned int A,
                unsigned int B,
                std::vector<cv::Point2f> &contour,
                std::vector<cv::Point2f> &rightArc,
                std::vector<double> &rightArcDistanceMap
                )
{
  unsigned int I=A;
  
  while(I!=B)
  {
    if(I==A)
    {
      rightArcDistanceMap.push_back(0);
      rightArc.push_back(contour[A]);
      if(I==contour.size()-1)
        I=0;
      else
        I++;
      continue;
    }

    rightArc.push_back(contour[I]);
    if(I==contour.size()-1)
    {
      I=0;
      rightArcDistanceMap.push_back(rightArcDistanceMap.back()+
                                    p2fdist(contour[0],contour.back()));
    }
    else
    {
      if(I==0)
      {
        rightArcDistanceMap.push_back(rightArcDistanceMap.back()+
                                     p2fdist(contour[0],contour.back()));
      }
      else
      {
        rightArcDistanceMap.push_back(rightArcDistanceMap.back()+
                                      p2fdist(contour[I],contour[I-1]));
      }
      ++I;
    }
  }
  if(B!=0)
    rightArcDistanceMap.push_back(rightArcDistanceMap.back()+
                                  p2fdist(contour[B],contour[B-1]));
  else
    rightArcDistanceMap.push_back(rightArcDistanceMap.back()+
                                  p2fdist(contour[B],contour.back()));
  rightArc.push_back(contour[B]);

}

void getNextPointInArc(std::vector<cv::Point2f> &arc,
                       std::vector<double> &arcDistanceMap,
                       double distance,
                       cv::Point2f &newPoint
    )
{
  unsigned int i=0;
  while(distance>arcDistanceMap[i] && i<arcDistanceMap.size())
    i++;

  if(i>0)
  {
    double d=distance-arcDistanceMap[i-1];
    double diip1=arcDistanceMap[i]-arcDistanceMap[i-1];
    double r=d/diip1;
    newPoint=arc[i-1]+r*(arc[i]-arc[i-1]);
    r=1;
  }

}

void getCandidateSpine(unsigned int A,
                       unsigned int B, 
                       std::vector<cv::Point2f> &contour,
                       std::vector<cv::Point2f> &spine,
                       std::vector<double> &widths,
                       std::vector<double> &angles,
                       std::vector<PointPair> &pPairs,
                       double &maxAngle,
                       int &maxAngleLocation,
                       double &firstHalfWidthsSum,
                       double &secondHalfWidthsSum,
                       double &width,
                       double &length,
                       unsigned int RESOLUTION=SPINE_SEGMENTS)
{
  std::vector<cv::Point2f> leftArc;
  std::vector<double> leftArcDistanceMap;
  std::vector<cv::Point2f> rightArc;
  std::vector<double> rightArcDistanceMap;
  std::vector<cv::Point2f> wholeSpine;

  getLeftArc(A,B,contour,leftArc,leftArcDistanceMap);
  getRightArc(A,B,contour,rightArc,rightArcDistanceMap);
  
  /*if((unsigned int)abs(leftArc.size()-rightArc.size())<contour.size()/4)
  {
    length=0;
    return;
  }*/
  
  double dleft=leftArcDistanceMap.back();
  double dright=rightArcDistanceMap.back();
  /*std::cout << "P: " << cv::arcLength(contour,true) << " l: " << dleft
            << " r: " << dright << std::endl;*/
  if(dleft > 2*dright || dright > 2*dleft) 
  {
    length=0;
    return;
  }
  double stepLeft=dleft/(RESOLUTION-1);
  double stepRight=dright/(RESOLUTION-1);
  spine.clear();
  wholeSpine.clear();
  wholeSpine.push_back(contour[A]);
  //std::cout << spine.back();
  //std::cout << std::endl;
  length=0;
  width=0;
  maxAngle=0;
  firstHalfWidthsSum=0;
  secondHalfWidthsSum=0;

  for(unsigned int i=1;i<RESOLUTION-1;i++)
  {
    cv::Point2f cPointLeft;
    cv::Point2f cPointRight;

    getNextPointInArc(leftArc,
                      leftArcDistanceMap,
                      i*stepLeft,
                      cPointLeft);

    getNextPointInArc(rightArc,
                      rightArcDistanceMap,
                      i*stepRight,
                      cPointRight);

    cv::Point2f newPoint=(cPointLeft+cPointRight)*0.5;
    pPairs.push_back(std::make_pair(cPointLeft,cPointRight));
    widths.push_back(p2fdist(cPointLeft,cPointRight));
    if(widths.back()>width)
    {
      width=widths.back();
    }
    if(i>=RESOLUTION/2)
      firstHalfWidthsSum+=widths.back();

    if(i<=RESOLUTION/2)
      secondHalfWidthsSum+=widths.back();

   // std::cout << cPointLeft << " + " << cPointRight << " -> " <<  newPoint << std::endl;
    wholeSpine.push_back(newPoint);

    if(i>=2)
    {
      angles.push_back(angleC(wholeSpine[i-2],wholeSpine[i-1],wholeSpine[i]));
      if(maxAngle<angles.back())
      {
        maxAngle=angles.back();
        maxAngleLocation=i-1;
      }
    }

  }
  unsigned int i=RESOLUTION-1;
  wholeSpine.push_back(contour[B]);
  angles.push_back(angleC(wholeSpine[i-2],wholeSpine[i-1],wholeSpine.back()));
  if(wholeSpine.size()==0)
  {
    length=-1;
    return;
  }
  if(maxAngle<angles.back())
  {
    maxAngle=angles.back();
    maxAngleLocation=i-1;
  }

  //std::cout << wholeSpine.back();
  //std::cout << std::endl;
  //std::cout << std::endl;
  length=cv::arcLength(wholeSpine,false);
  std::vector<double> spineDistanceMap;

  spineDistanceMap.push_back(0);
  spineDistanceMap.push_back(p2fdist(wholeSpine[1],wholeSpine[0]));
  for (unsigned int i=2;i<wholeSpine.size()-1;i++)
    spineDistanceMap.push_back(spineDistanceMap.back()+
                               p2fdist(wholeSpine[i],wholeSpine[i-1]));

  spineDistanceMap.push_back(length);

  double stepS=length/(SPINE_SEGMENTS-1);

  spine.push_back(wholeSpine[0]);
  for(unsigned int i=1;i<SPINE_SEGMENTS;i++)
  {
    cv::Point2f cPoint;
    getNextPointInArc(wholeSpine,
                      spineDistanceMap,
                      i*stepS,
                      cPoint);
    spine.push_back(cPoint);
  }
  spine.push_back(wholeSpine.back());
}

bool checkFinish(
                unsigned int &pointsCrossed,
                std::vector<cv::Point2f> &points,
                std::vector<cv::Point2f> &curSpine,
                unsigned int i,
                unsigned int &fwdIdx,
                unsigned int &bwdIdx,
                unsigned int &pfwdIdx,
                unsigned int &pbwdIdx,
                cv::Point2f &bwdPoint,
                cv::Point2f &fwdPoint,
                double &curDist
                )
{
      if(pointsCrossed==points.size()-1)
      {
        fwdIdx++;
        if(fwdIdx==points.size())
          fwdIdx=0;
        curSpine.push_back(points[fwdIdx]);
        curDist+=p2fdist(points[fwdIdx],curSpine.back());
        return true;
      }
      if(pointsCrossed==points.size()-2)
      {
        fwdIdx++;
        if(fwdIdx==points.size())
          fwdIdx=0;

        if(bwdIdx==0)
          bwdIdx=points.size();
        bwdIdx--;

        double fangle=angleD(curSpine[curSpine.size()-2],
               curSpine[curSpine.size()-1],
               points[fwdIdx]);

        double bangle=angleD(curSpine[curSpine.size()-2],
               curSpine[curSpine.size()-1],
               points[bwdIdx]);

        double MINANGLE=140;
        if (abs(fangle)<MINANGLE&& abs(bangle)<MINANGLE)
        {
          cv::Point2f np=(points[bwdIdx]+points[fwdIdx])*0.5;
          curDist+=p2fdist(np,curSpine.back());
          curSpine.push_back(np);
          return true;
        }
        else if(abs(fangle)<MINANGLE)
        {
          double distb=p2fdist(curSpine.back(),points[bwdIdx]);
          curSpine.push_back(points[bwdIdx]);
          curDist+=distb;
          return true;
        }
        else if(abs(bangle)<MINANGLE)
        {
          double distf=p2fdist(curSpine.back(),points[fwdIdx]);
          curSpine.push_back(points[fwdIdx]);
          curDist+=distf;
          return true;
        }
        else
        {
          double distf=p2fdist(curSpine.back(),points[fwdIdx]);
          double distb=p2fdist(curSpine.back(),points[bwdIdx]);

          if(distb>distf)
          {
            curSpine.push_back(points[bwdIdx]);
            curDist+=distb;
          }
          else
          {
            curSpine.push_back(points[fwdIdx]);
            curDist+=distf;
          }
          return true;
        }
      }
      return false;
}

void getPreviousPoint(std::vector<cv::Point2f> &a,
                  unsigned int idx,
                  cv::Point2f &n)
{
  if(idx==0)
    n=a[a.size()-1];
  else
    n=a[idx-1];
}
void smoothAndExtractCentripetal(cv::Mat &frame,
                                 std::vector<cv::Point2f> &ppoints,
                                 std::vector<cv::Point2f> &xpoints,
                                 std::vector<float> &d,
                                 float &csqrt)
{
  cv::Point2f NP=px3Smooth(frame, ppoints[ppoints.size()-2], ppoints.back(), ppoints[0], ppoints[1], ppoints[2]); 
  xpoints.push_back(NP);
  for(unsigned int i=1;i<ppoints.size();i++)
  {
    if(ppoints[i-1]!=ppoints[i])
    {
      cv::Point2f NP;
      if(i==1)
        NP=px3Smooth(frame,ppoints.back(),ppoints[i-1],ppoints[i],ppoints[i+1],ppoints[i+2]);
      else if(i<ppoints.size()-2)
        NP=px3Smooth(frame,ppoints[i-2],ppoints[i-1],ppoints[i],ppoints[i+1],ppoints[i+2]);
      else if (i==ppoints.size()-1)
        NP=px3Smooth(frame,ppoints[i-2],ppoints[i-1],ppoints[i],ppoints[0],ppoints[1]);
      else if (i==ppoints.size()-2)
        NP=px3Smooth(frame,ppoints[i-2],ppoints[i-1],ppoints[i],ppoints[i+1],ppoints[0]);

      if(NP!=xpoints.back())
      {
        csqrt+=sqrt(p2fdist(NP,xpoints.back()));
        xpoints.push_back(NP);
      }
    }
  }
  csqrt+=sqrt(p2fdist(xpoints.back(),xpoints[0]));
  d.push_back(0);
  for(unsigned int i=1;i<xpoints.size();i++)
  {
    float newD = d.back() + sqrt(p2fdist(xpoints[i-1],xpoints[i]))/csqrt;
    d.push_back(newD);
  }
  d.push_back(d.back()+sqrt(p2fdist(xpoints.back(),xpoints[0]))/csqrt);
}

void getNextPoint(std::vector<cv::Point2f> &a,
                  unsigned int idx,
                  cv::Point2f &n)
{
  if(idx==a.size()-1)
    n=a[0];
  else
    n=a[idx+1];
}

void fixContour(
    cvb::CvBlob &blob,
    larvaDistanceMap &Distances,
    unsigned int RES,
    cv::Mat &frame)
{
  std::vector<cv::Point2f> cntPoints;
  std::vector<cv::Point2f> ppoints;
  std::vector<cv::Point2f> newPoints;
  std::vector<float> x2;
  std::vector<float> y2;
  std::vector<float> d;
  //blobToPointVector(blob,cntPoints);
  //cv::approxPolyDP(cntPoints,ppoints,0.9,true);
  std::vector<cv::Point2f> xpoints;
  int PAD=2;
  blobToContourVector(blob,frame,PAD,ppoints);
  //double perimeter=cv::arcLength(xpoints,true);
  float csqrt=0.0;
  //smoothPoints(cntPoints,ppoints,frame,1);
  smoothAndExtractCentripetal(frame,ppoints,xpoints,d,csqrt);
  /*csqrt+=sqrt(p2fdist(xpoints[1],xpoints[0]));
  csqrt+=sqrt(p2fdist(xpoints[2],xpoints[1]));*/

  /*d.push_back(d.back()+sqrt(p2fdist(xpoints[1],xpoints[0]))/csqrt);*/
  /*d.push_back(d.back()+sqrt(p2fdist(xpoints[2],xpoints[1]))/csqrt);*/
  xpoints.push_back(xpoints[0]);
  xpoints.push_back(xpoints[1]);
  /*xpoints.push_back(xpoints[2]);*/
  std::vector<float> w;
  float AvgGrVal=0.0;
  std::vector<float> dx;
  std::vector<float> dy;
  std::vector<float> d2x;
  std::vector<float> d2y;
  std::vector<float> cvt;
  d2x.push_back(0);
  d2y.push_back(0);
  cvt.push_back(0);
  for(unsigned int i=0;i<xpoints.size();i++)
  {
    float cval=frame.at<cv::Vec3b>(xpoints[i])[0];
    AvgGrVal+=cval;
    w.push_back(cval);
    unsigned int pi;
    if(i==0)
    {
      pi=xpoints.size()-1;
    }
    double ratio=(d[i]-d[pi]);
    dx.push_back((xpoints[i].x-xpoints[pi].x)/ratio);
    dy.push_back((xpoints[i].y-xpoints[pi].y)/ratio);
    if(i>0)
    {
      d2x.push_back((dx[i]-dx[i-1])/ratio);
      d2y.push_back((dy[i]-dy[i-1])/ratio);
      cvt.push_back((dx[i]-d2x[i])/pow(dx[i]*dx[i]+dy[i]*dy[i],1.5));
    }
  }
  double ratio=(d[0]-d.back());
  d2x[0]=(dx[0]-dx.back())/ratio;
  d2y[0]=(dy[0]-dy.back())/ratio;
  cvt[0]=(dx[0]-d2x[0])/pow(dx[0]*dx[0]+dy[0]*dy[0],1.5);
  AvgGrVal=AvgGrVal/w.size();
  for(unsigned int i=0;i<w.size();i++)
  {
    w[i]=w[i]/AvgGrVal;
  }


std::vector<unsigned int> dIdx(10,0);
try{
spline2(xpoints,
           d,
           w,
           xpoints.size()-2,
           RES,
           newPoints,
           dIdx);
}
catch(...)
{
  std::cerr << "Spline Error [" << blob.label << ", " << xpoints.size() << "]" << std::endl;
  return;
}

  /*cv::Mat ROI=cv::Mat(frame,
              cv::Rect(blob.minx-PAD,
                       blob.miny-PAD,
                       blob.maxx-blob.minx+1+2*PAD,
                       blob.maxy-blob.miny+1+2*PAD
                      )
             );
  cv::Mat cROI,final;
  ROI.copyTo(cROI);
  cv::Point2f bp(blob.minx-0.5,blob.miny-0.5);
  for(unsigned int i=0;i<xpoints.size()-3;i++)
  {
  cv::Vec3b val= cROI.at<cv::Vec3b>((xpoints[i]-bp)+cv::Point2f(PAD,PAD));
  cv::Scalar sval(0,0,val[0]+val[1]+val[2]);
    cv::circle(cROI, 
        (xpoints[i]-bp)+cv::Point2f(PAD,PAD),
        0,
        sval,-1);
  }*/
  /*static int MULT=16;
  cv::resize(cROI,cROI,cv::Size(),MULT,MULT,cv::INTER_NEAREST);

  for(unsigned int i=0;i<newPoints.size();i++)
  {
    cv::circle(cROI, 
        MULT*(newPoints[i]-bp)+cv::Point2f(MULT*PAD,MULT*PAD),
        0,
        cv::Scalar(0,255,0),-1);
  }
  cv::circle(cROI, 
      MULT*(newPoints[0]-bp)+cv::Point2f(MULT*PAD,MULT*PAD),
      3,
      cv::Scalar(255,255,0),-1);

  for(unsigned int i=0;i<dIdx.size();i++)
  {
    cv::circle(cROI, 
        MULT*(newPoints[dIdx[i]]-bp)+cv::Point2f(MULT*PAD,MULT*PAD),
        2,
        cv::Scalar(0,255,0),-1);
  }*/
  double maxDist;
  int k=0;
  k++;
  std::vector<cv::Point2f> spine;
  double maxLength=-1;
  for (unsigned int i=0; i<dIdx.size()-1;i++)
  {
    for(unsigned int j=i+1; j<dIdx.size();j++)
    {
      /*unsigned int MAX=newPoints.size();
      if(dIdx[i]>dIdx[j] && 
          dIdx[i]-dIdx[j]>(MAX/3) &&
          dIdx[i]-(MAX+dIdx[j])>(MAX/3)
        )
        continue;
      else if(dIdx[j]>dIdx[i] && 
          dIdx[j]-dIdx[i]>(MAX/3) &&
          dIdx[j]-(MAX+dIdx[i])>(MAX/3)
          )
        continue;*/
      std::vector<cv::Point2f> cSpine;
      std::vector<PointPair> pPairs;

      if(dIdx[i]==0 || dIdx[j]==0)
        continue;
      larvaDistanceMap candidateMap(newPoints);
      getCandidateSpine(dIdx[i],
          dIdx[j],
          newPoints,
          candidateMap.Spine,
          candidateMap.Widths,
          candidateMap.Angles,
          candidateMap.spinePairs,
          candidateMap.maxAngle,
          candidateMap.maxAngleLocation,
          candidateMap.firstHalfWidthsSum,
          candidateMap.secondHalfWidthsSum,
          candidateMap.WidthDist,
          candidateMap.MaxDist,
          100);

      if(maxLength<=candidateMap.MaxDist)
      {
        maxLength=candidateMap.MaxDist;
        Distances=candidateMap;
      }
    }
  }
  maxDist=maxLength;

  /*for(unsigned int i=0;i<spine.size();i++)
  {
    cv::circle(cROI, 
        MULT*(spine[i]-bp)+cv::Point2f(MULT*PAD,MULT*PAD),
        2,
        cv::Scalar(0,0,255),-1);
  }

  if(blob.label==1)
  {
    int i;
    //copyMakeBorder( cROI, final, (1600-cROI.rows)/2 , (1600-cROI.rows)/2, (1600-cROI.cols)/2, (1600-cROI.cols)/2, cv::BORDER_CONSTANT , 0 );
    cv::imshow("LRV1",cROI);
    cv::waitKey(1);
    i=0;
  }*/

  Distances.MaxDistPoints.first=Distances.Spine.front();
  Distances.MaxDistPoints.second=Distances.Spine.back();
  Distances.p20.x=Distances.Spine[2].x;
  Distances.p20.y=Distances.Spine[2].y;
  Distances.p80.x=Distances.Spine[8].x;
  Distances.p80.y=Distances.Spine[8].y;
  Distances.MidPoint.x=Distances.Spine[5].x;
  Distances.MidPoint.y=Distances.Spine[5].y;
  return;
}

void computeSpine(
    cvb::CvBlob &blob,
    larvaDistanceMap &Distances,
    cv::Mat &frame
    )
{
  cv::Mat contour,lcontour;
  //createLarvaContour(contour,blob,CV_8UC1,0,false);
  //cv::resize(contour,lcontour,cv::Size(),16,16,cv::INTER_CUBIC);
  //std::vector<cv::Point2f> &Dpoints=Distances.points;
  std::vector<cv::Point2f> points;
  std::vector<cv::Point2f> procPoints;
  std::vector<cv::Point2f> PlainPoints;
  std::vector<float> t;
  float csqrt;
  int PAD=2;
  blobToContourVector(blob,frame,PAD,points);
  //cv::approxPolyDP(Dpoints,points,0.4,true);
  //smoothVec(points,PlainPoints,3,cv::Point2f(0,0));
  smoothAndExtractCentripetal(frame,points,PlainPoints,t,csqrt);
  std::vector<float> scurvature;
  curvVec(PlainPoints,t,scurvature);
  /*for(unsigned int i=0;i<PlainPoints.size();i++)
  {
  cv::Vec3b val= frame.at<cv::Vec3b>((PlainPoints[i]));
  cv::Scalar sval(0,0,val[0]+val[1]+val[2]);
    cv::circle(frame, 
        PlainPoints[i],
        0,
        sval,-1);
  }
  cv::imshow("test",frame);*/
  std::vector<unsigned int> dIdx(10,0);
  std::vector<float> dmax(10,0);
  cv::Point2f bp(blob.minx,blob.miny);

  getBestCurvature(scurvature,
      dIdx,
      dmax
      );
  
  //cv::approxPolyDP(Dpoints,procPoints,0.,true);
  //cv::approxPolyDP(Dpoints,PlainPoints,0.8,true);

  //std::cerr << "computeSpine: points: " << points.size() << std::endl;
  //int PADDING=2;
  //double RES=0.5;
  //Distances.WidthDist=DBL_MAX;
  /*PointPair wpoints;
  double maxDist=0;
  std::vector<double> widths;

  cv::Mat foo;
  std::map<double,int> anglesPoints;
  //frame.copyTo(foo);
  for (unsigned int i=0; i<PlainPoints.size(); ++i)
  {
    double A;
    if(i==0)
    {
      A=(angleC(PlainPoints.back(),PlainPoints[i],PlainPoints[i+1]));
    }
    else if(i==PlainPoints.size()-1)
    {
      A=(angleC(PlainPoints[i-1],PlainPoints[i],PlainPoints.front()));
    }
    else
    {
      A=(angleC(PlainPoints[i-1],PlainPoints[i],PlainPoints[i+1]));
    }
    while(anglesPoints.find(A)!=anglesPoints.end())
    {
      A=A+0.000001;
    }
    anglesPoints[A]=i;*/
    /*cv::circle(foo, 
      PlainPoints[i],
      0,
      cv::Scalar(0,0,255),-1);*/
 // }

  std::vector<unsigned int> &candidatePoints=dIdx;
  /*std::map<double,int>::iterator f=anglesPoints.begin();
  //frame.copyTo(foo);
  unsigned int maxel=2<PlainPoints.size()?2:PlainPoints.size();
  for(unsigned int i=0;i<maxel;i++)
  {
    bool jump=false;
    for(unsigned int j=0;j<candidatePoints.size();j++)
    {
      if(abs(candidatePoints[j]-f->second)<2)
      {
        jump=true;
        break;
      }
    }
    if (jump)
    {
      i--;
      f++;
      continue;
    }
    candidatePoints.push_back(f->second);*/
  /*cv::circle(foo, 
      PlainPoints[candidatePoints.back()],
      0,
      cv::Scalar(0,250,0),-1);
    f++;
  }*/

  int k=0;
  k++;
  double maxDist;
  double MaxLength=-1;
  for (unsigned int i=0; i<candidatePoints.size()-1;i++)
  {
    for(unsigned int j=i+1; j<candidatePoints.size();j++)
    {
      larvaDistanceMap candidateMap(PlainPoints);
      std::vector<cv::Point2f> candidateSpine;
      getCandidateSpine(candidatePoints[i],
                        candidatePoints[j],
                        PlainPoints,
                        candidateMap.Spine,
                        candidateMap.Widths,
                        candidateMap.Angles,
                        candidateMap.spinePairs,
                        candidateMap.maxAngle,
                        candidateMap.maxAngleLocation,
                        candidateMap.firstHalfWidthsSum,
                        candidateMap.secondHalfWidthsSum,
                        candidateMap.WidthDist,
                        candidateMap.MaxDist,
                        30);
      if(MaxLength<candidateMap.MaxDist)
      {
        MaxLength=candidateMap.MaxDist;
        Distances=candidateMap;
      }
        
    }
  }
  maxDist=MaxLength;

  
  std::vector<PointPair>::iterator SPit;
  /*unsigned int i=0;
  for(SPit=Distances.spinePairs.begin();
      SPit!=Distances.spinePairs.end();
      SPit++)
  {
    cv::circle(foo, 
        SPit->first,
        0,
        cv::Scalar(250,20*i,0),-1);
    cv::circle(foo, 
        SPit->second,
        0,
        cv::Scalar(250,20*i,0),-1);
    i++;
  }

  cv::imshow("endpoints",foo);
  cv::waitKey(1);*/

  if(Distances.Spine.size()==0)
  {
    std::cerr << "Compute Spine: No spine constructed." << std::endl;
    return;
  }
  if(Distances.Spine.front()==Distances.Spine.back())
  {
    std::cerr << "Head and Tail are the same!!!" << std::endl;
  }

  //std::vector<cv::Point2f> spinepoints=spine;
  Distances.MaxDistPoints.first=Distances.Spine.front();
  Distances.MaxDistPoints.second=Distances.Spine.back();
  Distances.p20.x=Distances.Spine[2].x;
  Distances.p20.y=Distances.Spine[2].y;
  Distances.p80.x=Distances.Spine[8].x;
  Distances.p80.y=Distances.Spine[8].y;
  Distances.MidPoint.x=Distances.Spine[5].x;
  Distances.MidPoint.y=Distances.Spine[5].y;
}

void computeInnerDistances(cvb::CvBlob &blob,
                           larvaDistanceMap &Distances,
                           cv::Point2f &MidPoint)
{
  cv::Mat contour;
  createLarvaContour(contour,blob);
  cv::Mat workingContour;
  contour.copyTo(workingContour);
  std::vector<cv::Point2f> &points=Distances.points;
  std::vector<cv::Point2f> SimplePoints;
  cv::approxPolyDP(points,SimplePoints,0.8,true);
  cv::Point2f MP(MidPoint.x+blob.minx,MidPoint.y+blob.miny);
  points=SimplePoints;
  double mWidth=0;
  Distances.WidthDist=DBL_MAX;

  int origNZ=countNonZero(contour);
  double MAX=0;
  for (unsigned int i=0; i<points.size(); ++i)
    {
      cv::Point2f p1=cv::Point2f(points[i].x,points[i].y);

      Distances[i][i]=0;
      for (unsigned int j=i+1; j<points.size(); ++j)
        {
          cv::Point2f p2(points[j].x,points[j].y);
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
              if(abs(j-i)>1)
                {
                  std::vector<cv::Point2f> widthArc;
                  widthArc.push_back(p1);
                  widthArc.push_back(MP);
                  widthArc.push_back(p2);
                  mWidth=cv::arcLength(widthArc, false);
                  if (Distances.WidthDist > mWidth)
                    {
                      Distances.WidthDist=mWidth;
                      //Distances.WidthDistPoints=p1p2;
                    }
                }
            }
          contour.copyTo(workingContour);
        }
    }

  for (unsigned int i=0; i<points.size(); ++i)
    {
      cv::Point2f p1(points[i].x,points[i].y);
      lBFS(i,points,Distances);
    }
}

