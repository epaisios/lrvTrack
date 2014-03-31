#include "larvaDistanceMap.hpp"
#include "lrvTrackDebug.hpp"
#include <queue>
#include <utility>
#include <opencv2/highgui/highgui.hpp>
//#include "gnuplot_i.hpp"
#include "blobUtils.hpp"
#include "alglib/ap.h"

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
              float sqrtres;
              sqrtres=sqrt(multres);
              double newDst = Distances[cur][i]+Distances[p1][cur] +
                              2*sqrtres;
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
  //cv::imshow("Spine3",lcontour);
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
    ++i;

  if(i>0)
  {
    double d=distance-arcDistanceMap[i-1];
    double diip1=arcDistanceMap[i]-arcDistanceMap[i-1];
    double r=d/diip1;
    newPoint=arc[i-1]+r*(arc[i]-arc[i-1]);
    r=1;
  }

}

double distanceBetweenSpines(larvaDistanceMap &a,larvaDistanceMap &b)
{
  double adist=0.0;
  double rdist=0.0;
  for(unsigned int i=0;i<a.Spine.size();i++)
  {
    adist+=p2fdist(a.Spine[i],b.Spine[i]);
    rdist+=p2fdist(a.Spine[i],b.Spine[a.Spine.size()-1-i]);
  }
  if(rdist>adist)
    return adist;
  else
    return rdist;
}

unsigned int findContourPointClosestTo(
                        cv::Point2f &p,
                        std::vector<cv::Point2f> &v
                        )
{
  double min=65535;
  unsigned int retval;
  unsigned int step=v.size()/10;
  for(unsigned int i=0;i<v.size();i+=step)
  {
    double cdist=p2fdist(p,v[i]);
    if(cdist<min)
    {
      min=cdist;
      retval=i;
    }
  }

  for(int i=retval-step;i<retval+step;++i)
  {
    unsigned int j=i;
    if(i<0)
      j=v.size()+i;
    if(i>v.size()-1)
      j=i-v.size();
    double cdist=p2fdist(v[j],p);
    if(cdist<min)
    {
      min=cdist;
      retval=j;
    }
  }
  return retval;
  
}

void getCandidateSpine(cvb::CvBlob &blob,
                       unsigned int A,
                       unsigned int B, 
                       std::vector<cv::Point2f> &contour,
                       larvaDistanceMap &candidateMap,
                       std::map<float,unsigned int> &curvatures,
                       unsigned int RESOLUTION=SPINE_SEGMENTS,
                       std::vector<cv::Point2f> *heads=NULL,
                       std::vector<cv::Point2f> *tails=NULL,
                       std::vector<cvb::CvBlob> *blobs=NULL)
{
  std::vector<cv::Point2f> leftArc;
  std::vector<double> leftArcDistanceMap;
  std::vector<cv::Point2f> rightArc;
  std::vector<double> rightArcDistanceMap;
  std::vector<cv::Point2f> wholeSpine;
  //cv::Mat FOO(1390,1038,CV_8UC3,cv::Scalar(0,0,0));
  
  getLeftArc(A,B,contour,leftArc,leftArcDistanceMap);
  getRightArc(A,B,contour,rightArc,rightArcDistanceMap);
  double dleft=leftArcDistanceMap.back();
  double dright=rightArcDistanceMap.back();
  /*std::cout << "P: " << cv::arcLength(contour,true) << " l: " << dleft
            << " r: " << dright << std::endl;*/
  if(dleft > 4*dright || dright > 4*dleft) 
  {
    candidateMap.MaxDist=0;
    return;
  }
  double stepLeft=dleft/(RESOLUTION-1);
  double stepRight=dright/(RESOLUTION-1);
  candidateMap.Spine.clear();
  wholeSpine.clear();
  wholeSpine.push_back(contour[A]);
  candidateMap.MaxDist=0;
  candidateMap.WidthDist=0;
  candidateMap.maxAngle=0;
  candidateMap.WidthAvg=0;
  candidateMap.firstHalfWidthsSum=0;
  candidateMap.secondHalfWidthsSum=0;

  std::map<float,unsigned int>::reverse_iterator curvIT=curvatures.rbegin();
  unsigned int CURVTEST=1;
  int TOTAL=curvatures.size();
  candidateMap.curvatureBias=0.0;
  for(unsigned int i=0; i<CURVTEST;++i)
  {
    unsigned int dA=std::min(abs(curvIT->second-A),TOTAL-abs(curvIT->second-A));
    unsigned int dB=std::min(abs(curvIT->second-B),TOTAL-abs(curvIT->second-B));
    if(dA<dB)
    {
      candidateMap.curvatureBias+=1.0/CURVTEST;
    }
    curvIT++;
  }

  for(unsigned int i=1;i<RESOLUTION-1;++i)
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
    candidateMap.spinePairs.push_back(std::make_pair(cPointLeft,cPointRight));
    candidateMap.Widths.push_back(p2fdist(cPointLeft,cPointRight));
    candidateMap.WidthAvg+=candidateMap.Widths.back();
    //cv::line(FOO,cPointLeft*2,cPointRight*2,cv::Scalar(0,0,200));
    if(candidateMap.Widths.back()>candidateMap.WidthDist)
    {
      candidateMap.WidthDist=candidateMap.Widths.back();
    }
    if(i>0.85*RESOLUTION)
      candidateMap.firstHalfWidthsSum+=candidateMap.Widths.back();

    if(i<0.15*RESOLUTION)
      candidateMap.secondHalfWidthsSum+=candidateMap.Widths.back();

    wholeSpine.push_back(newPoint);

  }
  candidateMap.WidthAvg=candidateMap.WidthAvg/candidateMap.Widths.size();
  for(unsigned int i=0;i<candidateMap.Widths.size();++i)
  {
    double D=(candidateMap.Widths[i]-candidateMap.WidthAvg);
    candidateMap.WidthSD+=D*D;
  }
  candidateMap.WidthSD=sqrt(candidateMap.WidthSD/candidateMap.Widths.size());

  unsigned int i=RESOLUTION-1;
  wholeSpine.push_back(contour[B]);
  if(wholeSpine.size()==0)
  {
    candidateMap.MaxDist=-1;
    return;
  }
  candidateMap.MaxDist=cv::arcLength(wholeSpine,false);
  std::vector<double> spineDistanceMap;

  spineDistanceMap.push_back(0);
  spineDistanceMap.push_back(p2fdist(wholeSpine[1],wholeSpine[0]));
  for (unsigned int i=2;i<wholeSpine.size()-1;++i)
    spineDistanceMap.push_back(spineDistanceMap.back()+
                               p2fdist(wholeSpine[i],wholeSpine[i-1]));

  spineDistanceMap.push_back(candidateMap.MaxDist);

  double stepS=candidateMap.MaxDist/(SPINE_SEGMENTS-1);

  candidateMap.Spine.push_back(wholeSpine[0]);
  for(unsigned int i=1;i<SPINE_SEGMENTS-1;++i)
  {
    cv::Point2f cPoint;
    getNextPointInArc(wholeSpine,
                      spineDistanceMap,
                      i*stepS,
                      cPoint);
    candidateMap.Spine.push_back(cPoint);

    if(i>=2)
    {
      candidateMap.Angles.push_back(
          angleD(
            candidateMap.Spine[i-2],
            candidateMap.Spine[i-1],
            candidateMap.Spine[i]));
      if(candidateMap.maxAngle<candidateMap.Angles.back())
      {
        candidateMap.maxAngle=candidateMap.Angles.back();
        candidateMap.maxAngleLocation=i-1;
      }
    }
  }
  candidateMap.Spine.push_back(wholeSpine.back());
  candidateMap.Angles.push_back(angleD(
                    candidateMap.Spine[SPINE_SEGMENTS-3],
                    candidateMap.Spine[SPINE_SEGMENTS-2],
                    candidateMap.Spine.back()));
  if(candidateMap.maxAngle<candidateMap.Angles.back())
  {
    candidateMap.maxAngle=candidateMap.Angles.back();
    candidateMap.maxAngleLocation=i-1;
  }

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

void expandcontour(cv::Mat &frame,
                   std::vector<cv::Point2f> &points,
                   std::vector<cv::Point2f> &outpoints)
{
  for(unsigned int i=0;i<points.size();++i)
  {
    //check points lrud
    cv::Point2f up(points[i].x,points[i].y-1);
    cv::Point2f down(points[i].x,points[i].y+1);
    cv::Point2f left(points[i].x-1,points[i].y);
    cv::Point2f right(points[i].x+1,points[i].y);

    cv::Vec3b uval= frame.at<cv::Vec3b>(cv::Point2f(up));
    cv::Vec3b dval= frame.at<cv::Vec3b>(cv::Point2f(down));
    cv::Vec3b lval= frame.at<cv::Vec3b>(cv::Point2f(left));
    cv::Vec3b rval= frame.at<cv::Vec3b>(cv::Point2f(right));
    
    if(uval[0]<dval[0] && uval[0] < lval[0] && uval[0] <rval[0])
    {
      outpoints.push_back(up);
    }
    else if(dval[0]<uval[0] && dval[0] < lval[0] && dval[0] < rval[0])
    {
      outpoints.push_back(down);
    }
    else if(lval[0]<uval[0] && lval[0] < dval[0] && lval[0] <rval[0])
    {
      outpoints.push_back(left);
    }
    else if(rval[0]<uval[0] && rval[0] < uval[0] && rval[0] <lval[0])
    {
      outpoints.push_back(right);
    }
    else
    {
      outpoints.push_back(points[i]);
    }
  }
}

void extractCentripetal(
                       std::vector<cv::Point2f> &ppoints,
                       std::vector<cv::Point2f> &xpoints,
                       std::vector<float> &d,
                       float &csqrt)
{
  unsigned int step=1;
  xpoints.push_back(ppoints[0]);
  for(unsigned int i=step;i<ppoints.size();i+=step)
  {
    cv::Point2f NP=ppoints[i];
    if(NP!=xpoints.back() && NP!=xpoints[0])
    {
      csqrt+=sqrt(p2fdist(NP,xpoints.back()));
      xpoints.push_back(NP);
    }
    else
    {
      if( i<ppoints.size()-1 && 
          ppoints[i+1]!=xpoints.back() && 
          ppoints[i+1] != xpoints[0] )
      {
        NP=ppoints[i+1];
        csqrt+=sqrt(p2fdist(NP,xpoints.back()));
        xpoints.push_back(NP);
      }
    }
  }
  csqrt+=sqrt(p2fdist(xpoints.back(),xpoints[0]));
  d.push_back(0);
  for(unsigned int i=1;i<xpoints.size();++i)
  {
    float newD = d.back() + sqrt(p2fdist(xpoints[i-1],xpoints[i]))/csqrt;
    d.push_back(newD);
  }
  d.push_back(d.back()+sqrt(p2fdist(xpoints.back(),xpoints[0]))/csqrt);
}

void smoothAndExtractCentripetal(cv::Mat &frame,
                                 std::vector<cv::Point2f> &ppoints,
                                 std::vector<cv::Point2f> &xpoints,
                                 std::vector<float> &d,
                                 float &csqrt)
{
  cv::Point2f NP=px5Smooth(frame, 
                           ppoints[ppoints.size()-2],
                           ppoints.back(), 
                           ppoints[0],
                           ppoints[1],
                           ppoints[2]); 
  xpoints.push_back(NP);
  for(unsigned int i=1;i<ppoints.size();++i)
  {
    if(ppoints[i-1]!=ppoints[i])
    {
      cv::Point2f NP;
      if(i==1)
        NP=px5Smooth(frame,ppoints.back(),ppoints[i-1],ppoints[i],ppoints[i+1],ppoints[i+2]);
      else if(i<ppoints.size()-2)
        NP=px5Smooth(frame,ppoints[i-2],ppoints[i-1],ppoints[i],ppoints[i+1],ppoints[i+2]);
      else if (i==ppoints.size()-1)
        NP=px5Smooth(frame,ppoints[i-2],ppoints[i-1],ppoints[i],ppoints[0],ppoints[1]);
      else if (i==ppoints.size()-2)
        NP=px5Smooth(frame,ppoints[i-2],ppoints[i-1],ppoints[i],ppoints[i+1],ppoints[0]);

      if(NP!=xpoints.back())
      {
        if(NP==xpoints[0])
          break;
        csqrt+=sqrt(p2fdist(NP,xpoints.back()));
        xpoints.push_back(NP);
      }
    }
  }
  csqrt+=sqrt(p2fdist(xpoints.back(),xpoints[0]));
  d.push_back(0);
  for(unsigned int i=1;i<xpoints.size();++i)
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
    cv::Mat &frame,
    cv::Mat &previousFrame,
    std::vector<cv::Point2f> *heads,
    std::vector<cv::Point2f> *tails,
    std::vector<cvb::CvBlob> *blobs)
{
  int ID=blob.label;
  std::vector<cv::Point2f> cntPoints;
  std::vector<cv::Point2f> ppoints;
  std::vector<cv::Point2f> newPoints;
  std::vector<float> x2;
  std::vector<float> y2;
  std::vector<float> d;
  cv::Point2f bp(blob.minx,blob.miny);
  //blobToPointVector(blob,cntPoints);
  //cv::approxPolyDP(cntPoints,ppoints,0.9,true);
  std::vector<cv::Point2f> xpoints;
  int PAD=3;
  blobToContourVector(blob,frame,PAD,cntPoints);
  cv::Point2f nh(0,0),nt(0,0);
  if(heads==(std::vector<cv::Point2f> *) -14)
  {
    cvb::CvBlob &pblob=blobs->back();
    cv::Point2f ph=heads->back();
    cv::Point2f pt=tails->back();
    cv::Mat cROI=cv::Mat(frame,
        cv::Rect(pblob.minx-PAD,
          pblob.miny-PAD,
          pblob.maxx-pblob.minx+1+2*PAD,
          pblob.maxy-pblob.miny+1+2*PAD
          )
        );
  cv::Mat ROI;
  cvtColor(cROI,ROI,CV_BGR2GRAY);
    cv::Mat pROI=cv::Mat(previousFrame,
        cv::Rect(pblob.minx-PAD,
          pblob.miny-PAD,
          pblob.maxx-pblob.minx+1+2*PAD,
          pblob.maxy-pblob.miny+1+2*PAD
          )
        );
    cv::Mat d1,d2;
    double diff1,diff2;
    double total1, total2;
    d1=ROI-pROI;
    d2=pROI-ROI;
    cv::Mat pcROI;
    cvtColor(pROI,pcROI,CV_GRAY2BGR);
    cv::circle(pcROI, 
      (ph)+cv::Point2f(PAD,PAD),
      0,
      cv::Scalar(0,255,0),-1);
    cv::circle(pcROI, 
      (pt)+cv::Point2f(PAD,PAD),
      0,
      cv::Scalar(0,255,0),-1);
    std::vector<uchar> status;
    std::vector<float> err;
    std::vector<cv::Point2f> pp;
    std::vector<cv::Point2f> pn;
    pp.push_back(ph);
    pp.push_back(pt);

    cv::TermCriteria termcrit(CV_TERMCRIT_ITER|CV_TERMCRIT_EPS, 200, 0.000000003);
    cv::Mat flow;
    calcOpticalFlowSF(pROI, ROI, flow,3,8,3);
    //calcOpticalFlowFarneback(pROI, ROI, flow, 0.8, 4, 4, 100,3,0.6,cv::OPTFLOW_FARNEBACK_GAUSSIAN);
    const cv::Point2f& fhxy = flow.at<cv::Point2f>(ph.y, ph.x);
    const cv::Point2f& ftxy = flow.at<cv::Point2f>(pt.y, pt.x);
    nh=cv::Point2f(ph.x+fhxy.x,ph.y+fhxy.y);
    nt=cv::Point2f(pt.x+ftxy.x,pt.y+ftxy.y);
    pn.push_back(nh);
    pn.push_back(nt);
    //calcOpticalFlowPyrLK(pROI, ROI, pp, pn, status, err, cv::Size(41,41),
    //    29,termcrit,cv::OPTFLOW_LK_GET_MIN_EIGENVALS);
    //nh=pn[0];
    //nt=pn[1];
    /*cv::circle(cROI, 
      (pn[0])+cv::Point2f(PAD,PAD),
      0,
      cv::Scalar(0,255,0),-1);
    cv::circle(cROI, 
      (pn[1])+cv::Point2f(PAD,PAD),
      0,
      cv::Scalar(0,255,0),-1);*/
  }
  //expandcontour(frame,ppoints,cntPoints);
  //ppoints.clear();
  pointsToContourVector(blob,cntPoints,frame,PAD,ppoints);

  //double perimeter=cv::arcLength(xpoints,true);
  float csqrt=0.0;
  //smoothPoints(cntPoints,ppoints,frame,1);
  //smoothAndExtractCentripetal(frame,ppoints,xpoints,d,csqrt);
  extractCentripetal(ppoints,xpoints,d,csqrt);
  /*csqrt+=sqrt(p2fdist(xpoints[1],xpoints[0]));
  csqrt+=sqrt(p2fdist(xpoints[2],xpoints[1]));*/

  /*d.push_back(d.back()+sqrt(p2fdist(xpoints[1],xpoints[0]))/csqrt);*/
  /*d.push_back(d.back()+sqrt(p2fdist(xpoints[2],xpoints[1]))/csqrt);*/
  //xpoints.push_back(xpoints[0]);
  //xpoints.push_back(xpoints[1]);
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
std::map<float,unsigned int> curvatures;
std::vector<unsigned int> dIdx(3,0);

  cv::Mat ROI=cv::Mat(frame,
              cv::Rect(blob.minx-PAD,
                       blob.miny-PAD,
                       blob.maxx-blob.minx+1+2*PAD,
                       blob.maxy-blob.miny+1+2*PAD
                      )
             );
  cv::Mat cROI,final;
  ROI.copyTo(cROI);
  for(unsigned int i=0;i<ppoints.size();++i)
  {
    cv::Vec3b val= cROI.at<cv::Vec3b>((ppoints[i]-bp)+cv::Point2f(PAD,PAD));
    cv::Scalar sval(0,0,val[0]+val[1]+val[2]);
    cv::circle(cROI, 
      (ppoints[i]-bp)+cv::Point2f(PAD,PAD),
      0,
      sval,-1);
  }
  static int MULT=24;
  //cv::resize(cROI,cROI,cv::Size(),MULT,MULT,cv::INTER_NEAREST);
  //cv::imshow("LRV1",cROI);

//if(heads==NULL)
  dIdx.resize(2);

std::vector<float> vcurv;
try{
spline4(xpoints,
           d,
           w,
           xpoints.size(),
           RES,
           newPoints,
           dIdx,
           curvatures,
           vcurv);
}
catch(alglib::ap_error &e)
{
  std::cerr << "Spline Error:" << e.msg << std::endl;
  std::cerr << printVector(xpoints) << std::endl;
  std::cerr << printVector(d) << std::endl;
  return;
}
catch(...)
{
  std::cerr << "Spline Error [" << blob.label << ", " << xpoints.size() << "]" << std::endl;
  return;
}


  bp=cv::Point2f(blob.minx-0.5,blob.miny-0.5);
  for(unsigned int i=0;i<newPoints.size();++i)
  {
    cv::circle(cROI, 
        MULT*(newPoints[i]-bp)+cv::Point2f(MULT*PAD,MULT*PAD),
        1,
        cv::Scalar(0,255,0),-1);
  }
  cv::circle(cROI, 
      MULT*(newPoints[0]-bp)+cv::Point2f(MULT*PAD,MULT*PAD),
      2,
      cv::Scalar(0,0,255),-1);

  double maxDist;
  larvaDistanceMap bestDist(newPoints);

  int k=0;
  k++;
  std::vector<cv::Point2f> spine;
  double maxLength=-1;
  if(!(nh.x==0 && nh.y==0 && nt.x==0 && nt.y == 0))
  {
    int hi,ti;
    cv::Point2f gnh=nh+bp;
    cv::Point2f gnt=nt+bp;
    hi=findContourPointClosestTo(gnh,newPoints);
    ti=findContourPointClosestTo(gnt,newPoints);
    dIdx[0]=hi;
    dIdx[1]=ti;
  }
  for (unsigned int i=0; i<dIdx.size()-1;++i)
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
      size_t di=dIdx[i];
      size_t dj=dIdx[j];
      if(blobs!=NULL && blobs->size()>0)
      {
        cv::Point2f bp((*blobs)[blobs->size()-1].minx,
            (*blobs)[blobs->size()-1].miny);
        cv::Point2f gH=heads->back()+bp;
        if(p2fdist(newPoints[dIdx[i]],gH)
            >
            p2fdist(newPoints[dIdx[j]],gH)
          )
        {
          di=dIdx[j];
          dj=dIdx[i];
        }
      }
      if(di == dj)
        continue;
      larvaDistanceMap candidateMap(newPoints);
      getCandidateSpine(blob,
          di,
          dj,
          newPoints,
          candidateMap,
          curvatures,
          50,
          heads,
          tails,
          blobs);

      if(maxLength<=candidateMap.MaxDist)
      {
        maxLength=candidateMap.MaxDist;
        bestDist=candidateMap;
      }
    }
  }

/*unsigned int A;
unsigned int B;
if( heads!=NULL && heads->size()>2)
{
  cv::Point2f pA=bestDist.Spine[0];
  cv::Point2f pB=bestDist.Spine.back();
  cv::Point2f bp((*blobs)[blobs->size()-1].minx,(*blobs)[blobs->size()-1].miny);
  cv::Point2f bp1((*blobs)[blobs->size()-2].minx,(*blobs)[blobs->size()-2].miny);
  cv::Point2f bp2((*blobs)[blobs->size()-3].minx,(*blobs)[blobs->size()-3].miny);
  cv::Point2f h1=(*heads)[heads->size()-1]+bp;
  cv::Point2f h2=(*heads)[heads->size()-2]+bp;
  cv::Point2f h3=(*heads)[heads->size()-3]+bp;
  cv::Point2f t1=(*tails)[tails->size()-1]+bp;
  cv::Point2f t2=(*tails)[tails->size()-2]+bp;
  cv::Point2f t3=(*tails)[tails->size()-3]+bp;

  if(fabs((pA-h1).x)+fabs((pA-h1).y)>fabs((pB-h1).x)+fabs((pB-h1).y) &&
      fabs((pA-t1).x)+fabs((pA-t1).y)<fabs((pB-t1).x)+fabs((pB-t1).y))
  {
    //pA=(pA+t1)*0.5;
    //pB=(pB+h1)*0.5;
    pA=(pA+t1+t2)*0.33333333333333;
    pB=(pB+h1+h2)*0.33333333333333;
    //pA=(pA+t1+t2+t3)*0.25;
    //pB=(pB+h1+h2+h3)*0.25;
    A=findContourPointClosestTo(pA,newPoints);
    B=findContourPointClosestTo(pB,newPoints);
  }
  else if(fabs((pA-h1).x)+fabs((pA-h1).y)<fabs((pB-h1).x)+fabs((pB-h1).y) &&
      fabs((pA-t1).x)+fabs((pA-t1).y)>fabs((pB-t1).x)+fabs((pB-t1).y))
  {
    //pA=(pA+h1)*0.5;
    //pB=(pB+t1)*0.5;
    pA=(pA+h1+h2)*0.33333333333333;
    pB=(pB+t1+t2)*0.33333333333333;
    //pA=(pB+h1+h2+h3)*0.25;
    //pB=(pA+t1+t2+t3)*0.25;
    A=findContourPointClosestTo(pA,newPoints);
    B=findContourPointClosestTo(pB,newPoints);
  }
  if(A==B)
  {
    Distances=bestDist;
    maxLength=Distances.MaxDist;
  }
  else{
    getCandidateSpine(blob,
        A,
        B,
        newPoints,
        Distances.Spine,
        Distances.Widths,
        Distances.Angles,
        Distances.spinePairs,
        Distances.maxAngle,
        Distances.maxAngleLocation,
        Distances.firstHalfWidthsSum,
        Distances.secondHalfWidthsSum,
        Distances.WidthDist,
        Distances.MaxDist,
        Distances.curvatureBias,
        curvatures,
        48,
        heads,
        tails,
        blobs);

    maxLength=Distances.MaxDist;
  }
}
else
{*/
  Distances=bestDist;
  maxLength=Distances.MaxDist;
//}

  maxDist=maxLength;

  for(unsigned int i=1;i<Distances.Spine.size();++i)
  {
    cv::circle(cROI, 
        MULT*(Distances.Spine[i]-bp)+cv::Point2f(MULT*PAD,MULT*PAD),
        3,
        cv::Scalar(0,0,255),-1);
    cv::line(cROI, 
        MULT*(Distances.Spine[i]-bp)+cv::Point2f(MULT*PAD,MULT*PAD),
        MULT*(Distances.Spine[i-1]-bp)+cv::Point2f(MULT*PAD,MULT*PAD),
        cv::Scalar(255,255,0),
        2
        );
  }
  cv::circle(cROI, 
      MULT*(Distances.Spine.back()-bp)+cv::Point2f(MULT*PAD,MULT*PAD),
      4,
      cv::Scalar(255,0,255),-1);

  for(unsigned int i=0;i<dIdx.size();++i)
  {
    cv::circle(cROI, 
        MULT*(newPoints[dIdx[i]]-bp)+cv::Point2f(MULT*PAD,MULT*PAD),
        3,
        cv::Scalar(0,255-(i*4),0),-1);
  }
  cv::circle(cROI, 
    MULT*(newPoints[dIdx.front()]-bp)+cv::Point2f(MULT*PAD,MULT*PAD),
    2,
    cv::Scalar(255,50,50),-1);
  
  cv::circle(cROI, 
      MULT*(newPoints[0]-bp)+cv::Point2f(MULT*PAD,MULT*PAD),
      5,
      cv::Scalar(0,0,255),-1);

  if(blob.label==0)
  {
    unsigned int i;

    cv::imshow("LRV1",cROI);
    for(i=0;i<Distances.Spine.size()-1;i++)
      std::cout << Distances.Spine[i] << ", " ;
    std::cout << Distances.Spine[i] << std::endl;

    for(i=0;i<Distances.Angles.size()-1;i++)
      std::cerr << Distances.Angles[i] << ", " ;
    std::cerr << Distances.Angles.back() << std::endl;
    cv::waitKey(1);
    //g1.remove_tmpfiles();
    //std::cerr << "dIdx: " <<  dIdx[0] << " , " << dIdx[1] << std::endl;
    i=0;
  }

  Distances.MaxDistPoints.first=Distances.Spine.front();
  Distances.MaxDistPoints.second=Distances.Spine.back();
  Distances.p20.x=Distances.Spine[2].x;
  Distances.p20.y=Distances.Spine[2].y;
  Distances.p80.x=Distances.Spine[9].x;
  Distances.p80.y=Distances.Spine[9].y;
  Distances.MidPoint=0.5*(Distances.Spine[5]+Distances.Spine[6]);
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
      /*getCandidateSpine(candidatePoints[i],
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
                        */
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

