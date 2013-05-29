#include "larvaDistanceMap.hpp"
#include <queue>
#include "blobUtils.hpp"

void lBFS(int p1, 
          std::vector<cv::Point> &Points ,
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
    for(int i=0; i<Points.size() ; i++)
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

void computeInnerDistances(cvb::CvBlob &blob,
                           larvaDistanceMap &Distances,
                           cv::Point &MidPoint)
{
  cv::Mat contour;
  createLarvaContour(contour,blob);
  cv::Mat workingContour;
  contour.copyTo(workingContour);
  std::vector<cv::Point> &points=Distances.points;
  std::vector<cv::Point> SimplePoints;
  cv::approxPolyDP(points,SimplePoints,0.9,true);
  cv::Point MP(MidPoint.x+blob.minx,MidPoint.y+blob.miny);
  points=SimplePoints;
  double mWidth=0;
  Distances.WidthDist=9999;

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
	if(abs(j-i)>1)
	{
	  std::vector<cv::Point> widthArc={p1,MP,p2};
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

  for (int i=0;i<points.size();i++)
  {
    cv::Point p1(points[i].x,points[i].y);
    lBFS(i,points,Distances);
  }
}

