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

bool checkFinish(
                unsigned int &pointsCrossed,
                std::vector<cv::Point2f > &points,
                std::vector<cv::Point2f > &curSpine,
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

void computeSpine(
    cvb::CvBlob &blob,
    larvaDistanceMap &Distances
    )
{
  cv::Mat contour,lcontour;
  std::vector<cv::Point2f> &Dpoints=Distances.points;
  std::vector<cv::Point2f> points;
  cv::approxPolyDP(Dpoints,points,0.4,true);
  //std::cerr << "computeSpine: points: " << points.size() << std::endl;
  //int PADDING=2;
  double RES=1.0;
  //Distances.WidthDist=DBL_MAX;
  PointPair wpoints;
  std::vector<cv::Point2f> spine;
  double maxDist=0;
  double maxWidth=0;
  std::vector<double> widths;

  for (unsigned int i=0; i<points.size()-1; ++i)
  {
    cv::Point2f fwdPoint;
    cv::Point2f bwdPoint;
    
    double curMaxWidth=0;
    std::vector<double> curWidths;
    PointPair curWidthPoints;

    std::vector<cv::Point2f> curSpine;  
    curSpine.push_back(points[i]);
    double curDist=0.0;

    unsigned int pfwdIdx=i;
    unsigned int fwdIdx=i+1;
    if(fwdIdx==points.size())
      fwdIdx=0;

    unsigned int pbwdIdx=i;
    unsigned int bwdIdx;
    if(i==0)
      bwdIdx=points.size()-1;
    else
      bwdIdx=i-1;

    double fwd_dst_fPoint_fwdIdx=p2fdist(points[i],points[fwdIdx]);
    double bwd_dst_bPoint_bwdIdx=p2fdist(points[i],points[bwdIdx]);

    double fwdRESLeft=RES;
    double bwdRESLeft=RES;

    unsigned int pointsCrossed=1;
    
    assignPoint(bwdPoint,points[i]);
    assignPoint(fwdPoint,points[i]);

    while(pointsCrossed<points.size())
    {
      if(distanceMatch(fwd_dst_fPoint_fwdIdx,fwdRESLeft))
      {
        assignPoint(fwdPoint,points[fwdIdx]);
        pfwdIdx=fwdIdx;
        pointsCrossed++;
        fwdIdx++;
        if(fwdIdx==points.size())
          fwdIdx=0;

        fwd_dst_fPoint_fwdIdx=p2fdist(fwdPoint,points[fwdIdx]);
      }
      else if(fwd_dst_fPoint_fwdIdx>fwdRESLeft) 
      {
        cv::Point2f lastP;
        assignPoint(lastP,fwdPoint);
        pointWithDist(lastP,points[fwdIdx],fwdPoint,fwdRESLeft);
        fwd_dst_fPoint_fwdIdx-=fwdRESLeft;
        fwdRESLeft=RES;
      }
      else if(fwd_dst_fPoint_fwdIdx<fwdRESLeft)
      {
        while(fwd_dst_fPoint_fwdIdx<fwdRESLeft && pointsCrossed<points.size()-2 )
        {
          fwdRESLeft-=fwd_dst_fPoint_fwdIdx;
          assignPoint(fwdPoint,points[fwdIdx]);
          pfwdIdx=fwdIdx;
          pointsCrossed++;
          fwdIdx++;
          if(fwdIdx==points.size())
            fwdIdx=0;
          
          fwd_dst_fPoint_fwdIdx=p2fdist(points[pfwdIdx],points[fwdIdx]);

        }
        if(checkFinish(pointsCrossed,
                       points,
                       curSpine,
                       i,
                       fwdIdx,
                       bwdIdx,
                       pfwdIdx,
                       pbwdIdx,
                       bwdPoint,
                       fwdPoint,
                       curDist))
        {
          break;
        }
        if(distanceMatch(fwd_dst_fPoint_fwdIdx,fwdRESLeft))
        {
          assignPoint(fwdPoint,points[fwdIdx]);
          pfwdIdx=fwdIdx;
          pointsCrossed++;
          fwdIdx++;
          if(fwdIdx==points.size())
            fwdIdx=0;

          fwd_dst_fPoint_fwdIdx=p2fdist(fwdPoint,points[fwdIdx]);
        }
        if(fwd_dst_fPoint_fwdIdx>fwdRESLeft) 
        {
          cv::Point2f lastP;
          assignPoint(lastP,fwdPoint);
          pointWithDist(lastP,points[fwdIdx],fwdPoint,fwdRESLeft);
          fwd_dst_fPoint_fwdIdx-=fwdRESLeft;
          fwdRESLeft=RES;
        }
      }

      if(distanceMatch(bwd_dst_bPoint_bwdIdx,bwdRESLeft))
      {
        assignPoint(bwdPoint,points[bwdIdx]);
        pbwdIdx=bwdIdx;
        pointsCrossed++;
        if(bwdIdx==0)
          bwdIdx=points.size();
        bwdIdx--;

        bwd_dst_bPoint_bwdIdx=p2fdist(bwdPoint,points[bwdIdx]);
      }
      else if(bwd_dst_bPoint_bwdIdx>bwdRESLeft) 
      {
        cv::Point2f lastP;
        assignPoint(lastP,bwdPoint);
        pointWithDist(lastP,points[bwdIdx],bwdPoint,bwdRESLeft);
        bwd_dst_bPoint_bwdIdx-=bwdRESLeft;
        bwdRESLeft=RES;
      }
      else if(bwd_dst_bPoint_bwdIdx<bwdRESLeft)
      {
        while(bwd_dst_bPoint_bwdIdx<bwdRESLeft && pointsCrossed<points.size()-2)
        {
          bwdRESLeft-=bwd_dst_bPoint_bwdIdx;
          assignPoint(bwdPoint,points[bwdIdx]);
          pbwdIdx=bwdIdx;
          pointsCrossed++;
          if(bwdIdx==0)
            bwdIdx=points.size();
          bwdIdx--;
          
          bwd_dst_bPoint_bwdIdx=p2fdist(points[pbwdIdx],points[bwdIdx]);

        }
        if(checkFinish(pointsCrossed,
              points,
              curSpine,
              i,
              fwdIdx,
              bwdIdx,
              pfwdIdx,
              pbwdIdx,
              bwdPoint,
              fwdPoint,
              curDist))
        {
          break;
        }
        if(distanceMatch(bwd_dst_bPoint_bwdIdx,bwdRESLeft))
        {
          assignPoint(bwdPoint,points[bwdIdx]);
          pbwdIdx=bwdIdx;
          pointsCrossed++;
          if(bwdIdx==0)
            bwdIdx=points.size();
          bwdIdx--;

          bwd_dst_bPoint_bwdIdx=p2fdist(bwdPoint,points[bwdIdx]);
        }
        if(bwd_dst_bPoint_bwdIdx>bwdRESLeft) 
        {
          cv::Point2f lastP;
          assignPoint(lastP,bwdPoint);
          pointWithDist(lastP,points[bwdIdx],bwdPoint,bwdRESLeft);
          bwd_dst_bPoint_bwdIdx-=bwdRESLeft;
          bwdRESLeft=RES;
        }
      }


      cv::Point2f newMP;
      newMP=(bwdPoint+fwdPoint)*0.5;
      double width=p2fdist(bwdPoint,fwdPoint);
      if(fwdIdx != i && 
          bwdIdx != i && 
          pfwdIdx != i && 
          pbwdIdx != i )
      {
        if (curMaxWidth<width)
        {
          curMaxWidth=width;
          curWidthPoints.first.x=bwdPoint.x;
          curWidthPoints.first.y=bwdPoint.y;
        }
        curWidths.push_back(width);
      }
      if(cv::pointPolygonTest(points,newMP,false)<=0)
      {
        curSpine.pop_back();
        break;
      }
      curDist+=p2fdist(newMP,curSpine.back());
      curSpine.push_back(newMP);

      if(checkFinish(pointsCrossed,
            points,
            curSpine,
            i,
            fwdIdx,
            bwdIdx,
            pfwdIdx,
            pbwdIdx,
            bwdPoint,
            fwdPoint,
            curDist))
      {
        break;
      }


    }
    if(curDist>maxDist)
    {
      spine=curSpine;
      maxDist=curDist;
      maxWidth=curMaxWidth;
      widths=curWidths;
      wpoints=curWidthPoints;
    }
  }


  //maxDist=arcLength(spine,false);
  //double dist80=0.8 * maxDist;
  //double dist20=0.2 * maxDist;
  //double dist50=0.5 * maxDist;

  cv::Point2f p20,p80,p50;
  
  if(spine.size()==0)
  {
    std::cerr << "Compute Spine: No spine constructed." << std::endl;
    return;
  }

  cv::Point2f aPoint;
  aPoint.x=spine.front().x;
  aPoint.y=spine.front().y;

  double cLength=0.0;
  double pLength=0.0;

  std::vector<cv::Point2f> &spinepoints=Distances.spinePoints;
  spinepoints[0].x=spine[0].x;
  spinepoints[0].y=spine[0].y;
  unsigned int j=1;

  for(unsigned i=1;i<spine.size();i++)
  {
    double xDist=aPoint.x-spine[i].x;
    double yDist=aPoint.y-spine[i].y;
    double lDist=sqrt(xDist*xDist+yDist*yDist);
    pLength=cLength;
    cLength+=lDist;
    double r=0.1;
    double splength=r*j*maxDist;

    while(cLength>splength && pLength<splength)
    {
      double xLength=r-pLength;
      double f=xLength/lDist;
      spinepoints[j].x=aPoint.x+f*(spine[i].x-aPoint.x);
      spinepoints[j].y=aPoint.y+f*(spine[i].y-aPoint.y);
      j++;
      splength=r*j*maxDist;
    }
    aPoint.x=spine[i].x;
    aPoint.y=spine[i].y;
  }
  spinepoints.back().x=spine.back().x;
  spinepoints.back().y=spine.back().y;

  //cv::Mat wMat(widths);
  //double wavg=cv::mean(wMat)[0];

  Distances.Spine=spine;
  Distances.MaxDist=maxDist;
  Distances.MaxDistPoints.first=Distances.Spine.front();
  Distances.MaxDistPoints.second=Distances.Spine.back();
  Distances.p20.x=spinepoints[2].x;
  Distances.p20.y=spinepoints[2].y;
  Distances.p80.x=spinepoints[8].x;
  Distances.p80.y=spinepoints[8].y;
  Distances.MidPoint.x=spinepoints[5].x;
  Distances.MidPoint.y=spinepoints[5].y;
  Distances.WidthDist=maxWidth;
  //Distances.WidthDist=(double) wavg;
  Distances.WidthDistPoints=wpoints;
  std::vector<std::vector<cv::Point2f> > cnts;
  cnts.push_back(spine);
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
                      Distances.WidthDistPoints=p1p2;
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

