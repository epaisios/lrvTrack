#include "lrvTrackDebug.hpp"
#include "lrvTrackFit.hpp"
#include "blobUtils.hpp"
#include <tbb/concurrent_hash_map.h>
#include <tbb/parallel_for.h>
#include <tbb/parallel_reduce.h>


inline double w(double x)
{
  double d=(x-1)*(x-1)*(x-1);
  if(x>=1.3)
    return sqrt(1-d*d);
  if(x<1.3)
    return sqrt(1-pow((1-x/1.3),1.7));
  return 0;
}

/*double optimize(cvb::CvBlob &blob, std::vector<larvaObject> &l)
{
  for(auto i=0;i<l.size();i++)
  {
    larvaFit f(l[i]);
    f.optimize(blob);
  }
}*/

ostream& operator<<(ostream& os,larvaFit &l)
{
  os << "[l:" << l.larvaFitData.length << ", w:" <<
        l.larvaFitData.width << ", mt:" <<
        l.larvaFitData.midtail << ", ga:" <<
        l.larvaFitData.global_angle << ", as:" <<
        printVector(l.larvaFitData.angles) << "]";
  return os;
}

void larvaFit::paintPoly(Mat &ROI, Point* f,size_t fsize)
{
  /*Point vi1[4];
  for(auto i=0;i<f.size();i++)
  {
    vi1[i]=Point(f[i].x,f[i].y);
  }*/

  stringstream pps;
  unsigned long long atime=getmsofday();
  paintPolyCalls++;
  fillConvexPoly(ROI,f,fsize,Scalar(255));
  unsigned long long diff=(getmsofday()-atime);
  paintPolyTime+=diff;
  if(diff>ppdiffmax)
  {
    ppdiffmax=diff;
    pps<<"MaxPP: Time: " << diff << " Points: [";
    for(auto i=0;i<fsize-1;i++)
      pps << f[i] << ", ";
    pps << f[fsize-1] << "]";
    BOOST_LOG_TRIVIAL(debug) << pps.str();
  }
}

larvaFit::larvaFit(const larvaFit &p)
{
  ID=p.ID;
  ipol=p.ipol;
  orig_width=p.orig_width;
  orig_length=p.orig_length;
  width_lookup=p.width_lookup;
  spine=p.spine;
  intSpine=p.intSpine;
  cpoints_size=p.cpoints_size;
  cpoints=(Point *) malloc(sizeof(Point)*(cpoints_size));
  degree_step=p.degree_step;
  wstep=p.wstep;
  ga=p.ga;
  a1=p.a1;
  a2=p.a2;
  a3=p.a3;
  a4=p.a4;
  wl=p.wl;
  PAD=p.PAD;
  ppdiffmax=p.ppdiffmax;
  createMatFromFitTime=p.createMatFromFitTime;
  createMatFromFitCalls=p.createMatFromFitCalls;
  paintPolyTime=p.paintPolyTime;
  paintPolyCalls=p.paintPolyCalls;
  setupSpineTime=p.setupSpineTime;
  setupSpineCalls=p.setupSpineCalls;
  calculateContourPointsTime=p.calculateContourPointsTime;
  calculateContourPointsCalls=p.calculateContourPointsCalls;
  larvaFitData=p.larvaFitData;
  //p.completeFrame.copyTo(completeFrame);
}

double larvaFit::errorFunc(Mat &b1C,
                           size_t minx,
                           size_t maxx,
                           size_t miny,
                           size_t maxy
                           )
{
  Mat l1C,Diff;
  
  createMatfromFit(l1C,minx,maxx,miny,maxy,true);
  bitwise_xor(b1C,l1C*2.5,Diff);
  size_t nz=countNonZero(Diff);
  //size_t bnz=countNonZero(b1C);
  return nz;
}

double larvaFit::errorFunc(cvb::CvBlob &blob)
{
  Mat b1C,l1C,Diff;
  
  createMatfromFit(l1C,blob.minx,blob.maxx,blob.miny,blob.maxy,true);
  createLarvaContour(b1C,blob,CV_8UC1,PAD,true,
      Scalar(255),8);
  bitwise_xor(b1C,l1C*2.5,Diff);
  size_t nz=countNonZero(Diff);
  //size_t lnz=countNonZero(l1C);
  return nz;
}

void larvaFit::pointToPoint(Point2f &p1, double angle, double d, Point2f &p2)
{
  //For Point spine[6]
  double k=tan(angle);
  double qpi=0.5*CV_PI;
  double dpi=1.5*CV_PI;

  if(angle<qpi || angle>dpi)
  {
    double xt=sqrt((d*d)/(1+k*k));
    p2.x=p1.x+xt;
    double yt=-k*xt;
    p2.y=p1.y+yt;
  }
  else if(angle>qpi || angle<dpi)
  {
    double xt=-sqrt((d*d)/(1+k*k));
    p2.x=p1.x+xt;
    double yt=-k*xt;
    p2.y=p1.y+yt;
  }
  else
  {
    p2.x=p1.x;
    p2.y=p1.y+(((angle-qpi)/qpi)-1)*d;
  }
}

void larvaFit::setloops(int ds,double ws)
{
  degree_step=ds*0.0174532925199432957692369076848;
  wstep=ws;
  setloops();
}

void larvaFit::setloops()
{
    ga.clear();
    //ga.push_back(-2*degree_step);
    ga.push_back(-degree_step);
    ga.push_back(0);
    ga.push_back(degree_step);
    //ga.push_back(2*degree_step);

    a1.clear();
    a1.push_back(-2*degree_step);
    a1.push_back(-degree_step);
    a1.push_back(0);
    a1.push_back(degree_step);
    a1.push_back(2*degree_step);
/*
    a2.clear();
    a2.push_back(-2*degree_step);
    a2.push_back(-degree_step);
    a2.push_back(0);
    a2.push_back(degree_step);
    a2.push_back(2*degree_step);
*/
    //a1=ga;
    a2=ga;
    a3=ga;
    a4=ga;


    wl.clear();
    //wl.push_back(-2*wstep);
    wl.push_back(-wstep);
    wl.push_back(0);
    wl.push_back(wstep);
    //wl.push_back(2*wstep);

}

void larvaFit::filterAngle(std::vector<double> &a, 
                           double &val,
                           double lim,
                           double add)
{
  if(val+a[0]<=lim)
    transform(a.begin(), a.end(), a.begin(),
        bind2nd(std::plus<double>(), add)); 
  else if(val+a.back()>=2*CV_PI-lim)
    transform(a.begin(), a.end(), a.begin(),
        bind2nd(std::minus<double>(), add));
}

void larvaFit::generate(std::vector<fitData> &fitSpace)
{
  /*BOOST_LOG_TRIVIAL(debug) << printVector(ga);
  BOOST_LOG_TRIVIAL(debug) << printVector(a1);
  BOOST_LOG_TRIVIAL(debug) << printVector(a2);
  BOOST_LOG_TRIVIAL(debug) << printVector(a3);
  BOOST_LOG_TRIVIAL(debug) << printVector(a4);
  BOOST_LOG_TRIVIAL(debug) << printVector(wl);*/

  //fitSpace=std::vector<fitData>(84376);
  //fitSpace=std::vector<fitData>(30376);
  //fitSpace=std::vector<fitData>(18226);

  filterAngle(a1,larvaFitData.angles[0],2*degree_step,2*degree_step);
  filterAngle(a2,larvaFitData.angles[1],2*degree_step,2*degree_step);
  filterAngle(a3,larvaFitData.angles[2],2*degree_step,degree_step);
  filterAngle(a4,larvaFitData.angles[3],2*degree_step,degree_step);
  //filterAngle(ga,larvaFitData.global_angle,degree_step);

  fitSpace=std::vector<fitData>((a1.size()*
                                a2.size()*
                                a3.size()*
                                a4.size()*
                                ga.size()*wl.size()*9)+1);
  //fitSpace=std::vector<fitData>(6076);
  
  fitData &l=larvaFitData;
  //BOOST_LOG_TRIVIAL(debug) << printVector(l.angles);

  fitSpace[0]=larvaFitData;
  size_t i=1;
  for(auto iag=0;iag<ga.size();++iag)
  {
    for(auto ia1=0; ia1<a1.size() ; ++ia1)
    {
      for(auto ia2=0; ia2<a2.size() ; ++ia2)
      {
        for(auto ia3=0; ia3<a3.size() ; ++ia3)
        {
          for(auto ia4=0; ia4<a4.size() ; ++ia4)
          {
            for(auto wp=0;wp<wl.size();wp++)
            {
              for(auto mpx=-1;mpx<2;mpx++)
              {
                for(auto mpy=-1;mpy<2;mpy++)
                {
                  Point2f mp=Point2f(larvaFitData.midtail.x+mpx,
                              larvaFitData.midtail.y+mpy);
                  /*BOOST_LOG_TRIVIAL(debug) << i << " ,"
                                           << iag << ", "
                                           << ia1 << ", "
                                           << ia2 << ", "
                                           << ia3 << ", "
                                           << ia4 << ", "
                                           << wp << ", "
                                           << mpx << " ,"
                                           << mpy << ", "
                                           << orig_length << ", "
                                           << orig_width << ", "
                                           << l.length << ", "
                                           << l.width;*/
                  fitSpace[i++]=fitData(mp,
                                        orig_length*(1.0+wl[wp]),
                                        orig_width*(1.0-wl[wp]),
                                        l.angles[0]+a1[ia1],
                                        l.angles[1]+a2[ia2],
                                        l.angles[2]+a3[ia3],
                                        l.angles[3]+a4[ia4],
                                        l.global_angle+ga[iag]
                                        );
                }
              }
            }
          }
        }
      }
    }
  }
}

class optimizeBody
{
  private:
    std::vector<larvaFit::fitData> &fitSpace;
    const larvaFit &lf;
    cvb::CvBlob &blob;
    Mat &b1C;
  public:
    double minerr;
    int mini;
    optimizeBody(std::vector<larvaFit::fitData> &fs, 
        larvaFit &l,
        cvb::CvBlob &b,
        Mat &b1
        ): 
      fitSpace(fs),
      lf(l),
      blob(b),
      b1C(b1),
      minerr(DBL_MAX),
      mini(-1)
  {}

    optimizeBody(optimizeBody &x, tbb::split) :
    fitSpace(x.fitSpace),
    lf(x.lf),
    blob(x.blob),
    b1C(x.b1C),
    minerr(DBL_MAX),
    mini(-1)
  {}
    void join(const optimizeBody &y)
    {
      if(y.minerr<=minerr)
      {
        minerr=y.minerr;
        mini=y.mini;
      }
    }

    void operator ()(const tbb::blocked_range<size_t>& range)
    {
      for (int i = range.begin(); i < range.end(); ++i)
      {
        larvaFit nl=lf;
        nl.larvaFitData=fitSpace[i];
        double r=nl.errorFunc(b1C,
            blob.minx,
            blob.maxx,
            blob.miny,
            blob.maxy
            );
        if(r<=minerr)
        {
          minerr=r;
          mini=i;
        }
      }
    }
};

double larvaFit::optimize(vector<cvb::CvBlob> &vec, cvb::CvBlob &res)
{
  innerCountersReset();
  std::vector<fitData> fitSpace;
  generate(fitSpace);
  std::vector<fitData>::iterator min=fitSpace.begin()+1;
  double minerr=DBL_MAX;
  Mat b1C;

  size_t minx=INT_MAX,maxx=0,miny=INT_MAX,maxy=0;
  for(auto &v:vec)
  {
    if(minx>v.minx)
      minx=v.minx;
    if(maxx<v.maxx)
      maxx=v.maxx;
    if(miny>v.miny)
      miny=v.miny;
    if(maxy<v.maxy)
      maxy=v.maxy;
  }
  Mat all;
  //Create contour bitmap of the resulting blob (the one detected)
  createLarvaContour_custom(b1C,res,CV_8UC1,
      minx,
      maxx,
      miny,
      maxy,
      PAD);

  //Create contour bitmap of the blobs of the corresponding larvae
  //in the previous step
  for(auto &v:vec)
  {
    if(v.label==ID)
      continue;
    Mat vmat;
    createLarvaContour_custom(vmat,
        v,
        CV_8UC1,
        minx,
        maxx,
        miny,
        maxy,
        PAD);
    if(all.empty())
      vmat.copyTo(all);
    else
      all=vmat | all;
  }

  unsigned long long atime=getmsofday();
  for(auto i=fitSpace.begin();i!=fitSpace.end();++i)
  {
    larvaFitData=*i;
    //cout << *this << endl;
    double r=errorFunc(b1C,
                       minx,
                       maxx,
                       miny,
                       maxy
                       );

    if(r<=minerr)
    {
      minerr=r;
      min=i;
    }
  }
  larvaFitData=*min;

  BOOST_LOG_TRIVIAL(debug) << "Computed Angles: " << printVector(larvaFitData.angles) << endl;

  Mat l1C;
  /*BOOST_LOG_TRIVIAL(debug) << "LarvaFitData: " << larvaFitData.width <<" ," 
                           << larvaFitData.length << ", " 
                           << larvaFitData.midtail << ", "  
                           << printVector(larvaFitData.angles) << ", "
                           << larvaFitData.global_angle 
                           << endl;*/
  createMatfromFit(l1C,all,minx,maxx,miny,maxy,true);
  cv::imshow("Fit Contour", l1C);
  //cv::imshow("Old Frame", greyFrame);
  cv::imshow("Old Contour", b1C);
  cv::moveWindow("Old Contour",200,200);
  cv::moveWindow("Fit Contour",200,130);
  cv::waitKey(1);
  return minerr;
  //setupSpine();
}

double larvaFit::optimize(cvb::CvBlob &blob)
{
  double minerr=DBL_MAX;
  std::vector<fitData> fitSpace;
  generate(fitSpace);
  Mat b1C;
  createLarvaContour(b1C,blob,CV_8UC1,PAD,true,
      Scalar(255),8);
  optimizeBody body(fitSpace,*this,blob,b1C);
  tbb::parallel_reduce(tbb::blocked_range<size_t>(0, fitSpace.size()), body);
  minerr=body.minerr;
  if(body.mini!=-1)
    larvaFitData=fitSpace[body.mini];
  else
    larvaFitData=fitSpace[0];
  Mat l1C;
  BOOST_LOG_TRIVIAL(debug) << "LarvaFitData: " << larvaFitData.width <<" ," 
                           << larvaFitData.length << ", " 
                           << larvaFitData.midtail << ", "  
                           << printVector(larvaFitData.angles) << ", "
                           << larvaFitData.global_angle;
  createMatfromFit(l1C,blob.minx,blob.maxx,blob.miny,blob.maxy,true);
  cv::Mat comb(l1C.rows,2*l1C.cols+2,l1C.type(),cv::Scalar(255));
  cv::Mat left(comb, Rect(0, 0, l1C.cols, l1C.rows));
  cv::Mat right(comb, Rect(l1C.cols+2, 0, l1C.cols, l1C.rows));
  l1C.copyTo(left);
  b1C.copyTo(right);
  //cv::imshow("Fit Contour", l1C);
  //cv::imshow("Old Frame", greyFrame);
  //cv::imshow("Old Contour", b1C);
  cv::imshow("combined",comb);
  stringstream n;
  vector<int> compression_params;
  compression_params.push_back(CV_IMWRITE_PNG_COMPRESSION);
  compression_params.push_back(0);
  n << "img/cc" << CURRENT_FRAME << ".png";
  imwrite(n.str(), comb, compression_params);
  //cv::moveWindow("Old Contour",200,200);
  //cv::moveWindow("Fit Contour",200,130);
  cv::waitKey(1);
  return minerr;
}

double larvaFit::optimize2(cvb::CvBlob &blob)
{
  innerCountersReset();
  std::vector<fitData> fitSpace;
  generate(fitSpace);
  std::vector<fitData>::iterator min=fitSpace.begin()+1;
  double minerr=DBL_MAX;
  Mat b1C;
  createLarvaContour(b1C,blob,CV_8UC1,PAD,true,
      Scalar(255),8);

  unsigned long long atime=getmsofday();
  for(auto i=fitSpace.begin();i!=fitSpace.end();++i)
  {
    larvaFitData=*i;
    //cout << *this << endl;
    double r=errorFunc(b1C,
                       blob.minx,
                       blob.maxx,
                       blob.miny,
                       blob.maxy
                       );

    if(r<=minerr)
    {
      minerr=r;
      min=i;
    }
  }
  larvaFitData=*min;

  BOOST_LOG_TRIVIAL(debug) << "Computed Angles: " << printVector(larvaFitData.angles) << endl;

  Mat l1C;
  /*BOOST_LOG_TRIVIAL(debug) << "LarvaFitData: " << larvaFitData.width <<" ," 
                           << larvaFitData.length << ", " 
                           << larvaFitData.midtail << ", "  
                           << printVector(larvaFitData.angles) << ", "
                           << larvaFitData.global_angle 
                           << endl;*/
  createMatfromFit(l1C,blob.minx,blob.maxx,blob.miny,blob.maxy,true);
  cv::imshow("Fit Contour", l1C);
  //cv::imshow("Old Frame", greyFrame);
  cv::imshow("Old Contour", b1C);
  cv::moveWindow("Old Contour",200,200);
  cv::moveWindow("Fit Contour",200,130);
  cv::waitKey(1);
  return minerr;
  //setupSpine();
}

void larvaFit::showContour()
{
  Mat l1C;
  /*BOOST_LOG_TRIVIAL(debug) << "LarvaFitData: " << larvaFitData.width <<" ," 
                           << larvaFitData.length << ", " 
                           << larvaFitData.midtail << ", "  
                           << printVector(larvaFitData.angles) << ", "
                           << larvaFitData.global_angle 
                           << endl;*/
  createMatfromFit(l1C,0,800,0,800,true);
  cv::imshow("Contour", l1C);
  cv::imshow("Picture", greyFrame);
  cv::waitKey(1);
}

larvaFit::larvaFit()
{
  cpoints_size=0;
}

void larvaFit::setup(larvaObject &l,int dstep, double wstep, size_t FRAME)
{
  setup(l,FRAME);
  setloops(dstep,wstep);
}  

larvaFit::larvaFit(larvaObject &l,int dstep, double wstep, size_t FRAME)
{
  larvaFit(l,FRAME);
  setloops(dstep,wstep);
}  

larvaFit::~larvaFit()
{
  if(cpoints_size>0)
    free(cpoints);
}

larvaFit::larvaFit(larvaObject &l,size_t FRAME)
{
  setup(l,FRAME);
}

void larvaFit::setup(larvaObject &l,size_t FRAME)
{
  ID=l.larva_ID;
  innerCountersReset();
  spine=l.lrvDistances.back().Spine;
  intSpine=vector<Point2f>(ipol*(spine.size()-1)+1,Point2f(-1,-1));
  cpoints=(Point *) malloc(sizeof(Point)*(2*(intSpine.size()-1)));
  cpoints_size=2*(intSpine.size()-1);
  //completeFrame=cv::Mat(FRAME_ROWS,FRAME_COLS,
  //                        CV_8UC1,Scalar(0));
  double width=l.width.back();
  orig_width=width/2;
  double length=l.length.back();
  orig_length=length;
  Point2f bp(l.blobs.back().minx,l.blobs.back().miny);
  for(auto i=0;i<intSpine.size();i++)
  {
    width_lookup.push_back(w(2.0*i/(intSpine.size()-1)));
  }
  
  //BOOST_LOG_TRIVIAL(debug) << "lrvfit: Spine: " << printVector(p) << endl;
  //BOOST_LOG_TRIVIAL(debug) << "lrvfit: Head: " << l.heads.back()+bp << endl;
  //BOOST_LOG_TRIVIAL(debug) << "lrvfit: Tail: " << l.tails.back()+bp << endl;

  Point2f V1=spine[8]-spine[6];
  double angleg=atan2(V1.y,-V1.x);
  if(angleg<0)
    angleg=2*CV_PI+angleg;
  vector<double> angles(4,0.0);
  angles[0]=((l.lrvDistances.back().Angles[1]));
  angles[1]=((l.lrvDistances.back().Angles[3]));
  angles[2]=((l.lrvDistances.back().Angles[5]));
  angles[3]=((l.lrvDistances.back().Angles[7]));
  larvaFitData=fitData(spine[8],
                        length,
                        width/2,
                        angles,
                        angleg);
  setloops();
  /*BOOST_LOG_TRIVIAL(debug) << "lrvfit: Width: " << width << endl;
  BOOST_LOG_TRIVIAL(debug) << "lrvfit: Length: " << length << endl;
  BOOST_LOG_TRIVIAL(debug) << "lrvfit: Angleg: " << angleg*180/CV_PI << endl;
  BOOST_LOG_TRIVIAL(debug) << "lrvfit: MidTail: " << spine[8] << endl;
  BOOST_LOG_TRIVIAL(debug) << "lrvfit: Angles: " << printVector(angles) << endl;*/

}

void larvaFit::setupSpine()
{
  unsigned long long atime=getmsofday();
  setupSpineCalls++;
  double seglen=larvaFitData.length/11;
  spine[8].x=larvaFitData.midtail.x;
  spine[8].y=larvaFitData.midtail.y;
  /*BOOST_LOG_TRIVIAL(debug) << "LarvaFitData: " << larvaFitData.width <<" ," 
                           << larvaFitData.length << ", " 
                           << larvaFitData.midtail << ", "  
                           << printVector(larvaFitData.angles) << ", "
                           << larvaFitData.global_angle 
                           << endl;*/

  //For point spine[6]
  pointToPoint(spine[8],larvaFitData.global_angle,seglen*2,spine[6]);
  //For point spine[11]
  double phi=larvaFitData.global_angle+(2*CV_PI-larvaFitData.angles[3]);
  if(phi>=2*CV_PI)
    phi-=2*CV_PI;
  if(phi<0)
    phi+=2*CV_PI;
  pointToPoint(spine[8],phi,seglen*3,spine[11]);
  //For point spine[10] we interpolate
  spine[9]=((spine[11]-spine[8])*(1.0/3.0))+spine[8];
  spine[10]=(spine[9]+spine[11])*0.5;

  //For point spine[4]
  phi=(2*CV_PI-larvaFitData.angles[2])-(CV_PI-larvaFitData.global_angle);
  if(phi>=2*CV_PI)
    phi-=2*CV_PI;
  if(phi<0)
    phi+=2*CV_PI;
  pointToPoint(spine[6],phi,seglen*2,spine[4]);

  //For point spine[2]
  phi=(2*CV_PI-larvaFitData.angles[1])-(CV_PI-phi);
  if(phi>=2*CV_PI)
    phi-=2*CV_PI;
  if(phi<0)
    phi+=2*CV_PI;
  pointToPoint(spine[4],phi,seglen*2,spine[2]);

  //For point spine[0]
  phi=(2*CV_PI-larvaFitData.angles[0])-(CV_PI-phi);
  if(phi>2*CV_PI)
    phi-=2*CV_PI;
  if(phi<0)
    phi+=2*CV_PI;
  pointToPoint(spine[2],phi,seglen*2,spine[0]);

  spine[7]=(spine[6]+spine[8])*0.5;
  spine[5]=(spine[4]+spine[6])*0.5;
  spine[3]=(spine[2]+spine[4])*0.5;
  spine[1]=(spine[0]+spine[2])*0.5;
 
  //BOOST_LOG_TRIVIAL(debug) << "Setup Spine: Computed Spine: " << printVector(spine) << endl;
  setupSpineTime+=(getmsofday()-atime);
}
void larvaFit::innerCountersReset()
{
    ppdiffmax=0;
    createMatFromFitTime=0;
    createMatFromFitCalls=0;
    paintPolyTime=0;
    paintPolyCalls=0;
    setupSpineTime=0;
    setupSpineCalls=0;
    calculateContourPointsTime=0;
    calculateContourPointsCalls=0;
}

void larvaFit::createMatfromFit(Mat &larvaFitContour,
                                Mat &fitBase,
                                size_t minx,
                                size_t maxx,
                                size_t miny,
                                size_t maxy,
                                bool verbose
                                )
{
  unsigned long long atime=getmsofday();
  createMatFromFitCalls++;
  setupSpine();
  Mat tmp=cv::Mat(maxy-miny+1+(2*PAD),
        maxx-minx+1+(2*PAD),
        CV_8UC1,Scalar(0));
  /*tmp=Mat(completeFrame,
      Rect(minx-PAD,
        miny-PAD,
        maxx-minx+1+2*PAD,
        maxy-miny+1+2*PAD));*/
  Point2f bp(minx-PAD,miny-PAD);
  
  if(fitBase.rows!=tmp.rows || fitBase.cols!=tmp.cols)
    return;

  //tmp.setTo(Scalar::all(0));
  tmp=Scalar(0);
  //cv::Mat completeFrame=cv::Mat(FRAME_ROWS,FRAME_COLS,
  //                        CV_8UC1,Scalar(0));
  for(size_t i=0; i<spine.size()-1;i++)
  {
    size_t v=i*ipol;
    for(auto j=0;j<ipol;j++)
    {
      intSpine[v+j]=(spine[i+1]-spine[i])*((double)j/ipol)+spine[i];
    }
  }
  intSpine.back()=spine.back();
  //if(verbose)
    //BOOST_LOG_TRIVIAL(debug) << "Interpolated Spine: " << printVector(intSpine) << endl;
  //vector<Point2f> cpoints(2*(intSpine.size()-1));
  cpoints[0]=intSpine[0]-bp;
  cpoints[cpoints_size/2]=intSpine.back()-bp;
  for(auto i=1;i<intSpine.size()-1;i++)
  {
    calculateContourPoints(intSpine[i-1],
                           intSpine[i],
                           intSpine[i+1],
                           bp,
                           //2.0*(i+1.0)/intSpine.size(),
                           i,
                           larvaFitData.width,
                           cpoints[i],
                           cpoints[cpoints_size-i]);
    /*if(i==1)
    {
      Point t[3]={cpoints[i-1],
                      cpoints[i],
                      cpoints[cpoints.size()-i]};
      //fillConvexPoly(tmp,t,3,Scalar(255));
      paintPoly(completeFrame,t,3);
    }
    else if(i==intSpine.size()-2)
    {
      Point2f t1[4]={cpoints[i-1],
                      cpoints[i],
                      cpoints[cpoints.size()-i+1],
                      cpoints[cpoints.size()-i]};
      vector<Point2f> tv1(&t1[0],&t1[0]+4);
      paintPoly(completeFrame,tv1);

      Point t1[4]={cpoints[i-1],
                      cpoints[i],
                      cpoints[cpoints.size()-i],
                      cpoints[cpoints.size()-i+1]};
      //fillConvexPoly(tmp,t1,4,Scalar(255));
      paintPoly(completeFrame,t1,4);

      Point t2[3]={cpoints[cpoints.size()/2],
                      cpoints[i],
                      cpoints[cpoints.size()-i]};
      //fillConvexPoly(tmp,t2,3,Scalar(255));
      paintPoly(completeFrame,t2,3);
    }
    else if(i>1)
    {
      Point2f t1[4]={cpoints[i-1],
                      cpoints[i],
                      cpoints[cpoints.size()-i+1],
                      cpoints[cpoints.size()-i]};
      vector<Point2f> tv1(&t1[0],&t1[0]+4);
      paintPoly(completeFrame,tv1);
      Point t1[4]={cpoints[i-1],
                      cpoints[i],
                      cpoints[cpoints.size()-i],
                      cpoints[cpoints.size()-i+1]};
      //fillConvexPoly(tmp,t1,4,Scalar(255));
      paintPoly(completeFrame,t1,4);
    }*/
    /*for(size_t i=0;i<intSpine.size();i++)
    {
      circle(completeFrame,intSpine[i],0,Scalar(100),-1);
    }*/
  }
  //int npts=cpoints.size();
  fillPoly(tmp,(const Point**) &cpoints,(int *) &cpoints_size,1,Scalar(255));
  for(size_t i=1;i<spine.size();i++)
  {
    line(tmp,spine[i-1]-bp,spine[i]-bp,Scalar(100));
  }
  /*for(size_t i=0;i<spine.size();i++)
  {
    circle(completeFrame,spine[i],0,Scalar(50),-1);
  }*/
  //if(verbose)
    //BOOST_LOG_TRIVIAL(debug) << "CPOINTS: " << printVector(cpoints);
  larvaFitContour = tmp | fitBase;
  createMatFromFitTime+=(getmsofday()-atime);
  //free(t1);
}

void larvaFit::createMatfromFit(Mat &larvaFitContour,
                                size_t minx,
                                size_t maxx,
                                size_t miny,
                                size_t maxy,
                                bool verbose
                                )
{
  unsigned long long atime=getmsofday();
  createMatFromFitCalls++;
  setupSpine();
  Mat tmp=cv::Mat(maxy-miny+1+(2*PAD),
        maxx-minx+1+(2*PAD),
        CV_8UC1,Scalar(0));
  /*tmp=Mat(completeFrame,
      Rect(minx-PAD,
        miny-PAD,
        maxx-minx+1+2*PAD,
        maxy-miny+1+2*PAD));*/
  Point2f bp(minx-PAD,miny-PAD);
  //tmp.setTo(Scalar::all(0));
  tmp=Scalar(0);
  //cv::Mat completeFrame=cv::Mat(FRAME_ROWS,FRAME_COLS,
  //                        CV_8UC1,Scalar(0));
  for(size_t i=0; i<spine.size()-1;i++)
  {
    size_t v=i*ipol;
    for(auto j=0;j<ipol;j++)
    {
      intSpine[v+j]=(spine[i+1]-spine[i])*((double)j/ipol)+spine[i];
    }
  }
  intSpine.back()=spine.back();
  //if(verbose)
    //BOOST_LOG_TRIVIAL(debug) << "Interpolated Spine: " << printVector(intSpine) << endl;
  //vector<Point2f> cpoints(2*(intSpine.size()-1));
  cpoints[0]=intSpine[0]-bp;
  cpoints[cpoints_size/2]=intSpine.back()-bp;
  for(auto i=1;i<intSpine.size()-1;i++)
  {
    calculateContourPoints(intSpine[i-1],
                           intSpine[i],
                           intSpine[i+1],
                           bp,
                           //2.0*(i+1.0)/intSpine.size(),
                           i,
                           larvaFitData.width,
                           cpoints[i],
                           cpoints[cpoints_size-i]);
    /*if(i==1)
    {
      Point t[3]={cpoints[i-1],
                      cpoints[i],
                      cpoints[cpoints.size()-i]};
      //fillConvexPoly(tmp,t,3,Scalar(255));
      paintPoly(completeFrame,t,3);
    }
    else if(i==intSpine.size()-2)
    {
      Point2f t1[4]={cpoints[i-1],
                      cpoints[i],
                      cpoints[cpoints.size()-i+1],
                      cpoints[cpoints.size()-i]};
      vector<Point2f> tv1(&t1[0],&t1[0]+4);
      paintPoly(completeFrame,tv1);

      Point t1[4]={cpoints[i-1],
                      cpoints[i],
                      cpoints[cpoints.size()-i],
                      cpoints[cpoints.size()-i+1]};
      //fillConvexPoly(tmp,t1,4,Scalar(255));
      paintPoly(completeFrame,t1,4);

      Point t2[3]={cpoints[cpoints.size()/2],
                      cpoints[i],
                      cpoints[cpoints.size()-i]};
      //fillConvexPoly(tmp,t2,3,Scalar(255));
      paintPoly(completeFrame,t2,3);
    }
    else if(i>1)
    {
      Point2f t1[4]={cpoints[i-1],
                      cpoints[i],
                      cpoints[cpoints.size()-i+1],
                      cpoints[cpoints.size()-i]};
      vector<Point2f> tv1(&t1[0],&t1[0]+4);
      paintPoly(completeFrame,tv1);
      Point t1[4]={cpoints[i-1],
                      cpoints[i],
                      cpoints[cpoints.size()-i],
                      cpoints[cpoints.size()-i+1]};
      //fillConvexPoly(tmp,t1,4,Scalar(255));
      paintPoly(completeFrame,t1,4);
    }*/
    /*for(size_t i=0;i<intSpine.size();i++)
    {
      circle(completeFrame,intSpine[i],0,Scalar(100),-1);
    }*/
  }
  //int npts=cpoints.size();
  fillPoly(tmp,(const Point**) &cpoints,(int *) &cpoints_size,1,Scalar(255));
  for(size_t i=1;i<spine.size();i++)
  {
    line(tmp,spine[i-1]-bp,spine[i]-bp,Scalar(100));
  }
  /*for(size_t i=0;i<spine.size();i++)
  {
    circle(completeFrame,spine[i],0,Scalar(50),-1);
  }*/
  //if(verbose)
    //BOOST_LOG_TRIVIAL(debug) << "CPOINTS: " << printVector(cpoints);
  tmp.copyTo(larvaFitContour);
  createMatFromFitTime+=(getmsofday()-atime);
  //free(t1);
}

void larvaFit::calculateContourPoints(Point2f &a,
                            Point2f &b,
                            Point2f &c,
                            double b_index,
                            double width,
                            Point &cl,
                            Point &cr)
{
  unsigned long long atime=getmsofday();
  calculateContourPointsCalls++;
  //Make b our reference point
  Point2f r=a-b;
  Point2f f=c-b;
  Point2f sp=r+f;
  Point2f lPoint,rPoint;
  width=width_lookup[b_index]*width;
  if(fabs(sp.x)<0.0001 && fabs(sp.y)<0.0001) // On the same line
  {
    //We need the perpendicular vector
    if (r.x!=0 && r.y!=0)
    {
      double pslope=-r.x/r.y;
      double xc1=sqrt((width*width)/(1+pslope*pslope));
      lPoint.x=xc1;
      lPoint.y=pslope*xc1;
      rPoint.x=-xc1;
      rPoint.y=-pslope*xc1;
    }
    else if(r.x==0 && r.y==0)
    {
      cerr << "Weird case" << endl;
      exit(0);
    }
    else if(r.x==0)
    {
      lPoint.x=width;
      lPoint.y=0;
      rPoint.x=-width;
      rPoint.y=0;
    }
    else if(r.y==0)
    {
      lPoint.x=0;
      lPoint.y=width;
      rPoint.x=0;
      rPoint.y=-width;
    }
  
  }
  else
  {
    double ratio=width/sqrt(sp.x*sp.x+sp.y*sp.y);
    lPoint.x=ratio*sp.x;
    lPoint.y=ratio*sp.y;
    rPoint.x=-lPoint.x;
    rPoint.y=-lPoint.y;
  }

  double avals[4]={atan2f(r.y,r.x),
                    atan2f(-r.y,-r.x),
                    atan2f(rPoint.y,rPoint.x),
                    atan2f(lPoint.y,lPoint.x)};
  if(avals[0]<0)
    avals[0]=2*CV_PI+avals[0];
  if(avals[1]<0)
    avals[1]=2*CV_PI+avals[1];
  if(avals[2]<0)
    avals[2]=2*CV_PI+avals[2];
  if(avals[3]<0)
    avals[3]=2*CV_PI+avals[3];
  if((avals[2]>avals[0] || avals[2]<avals[1]) &&
     (avals[3]<avals[0] || avals[3]>avals[1]))
  {
    //points are correct lPoint is to the left, rPoint is to the right
    cl=lPoint+b;
    cr=rPoint+b;
  }
  else if((avals[2]<=avals[0] || avals[2]>=avals[1]) &&
          (avals[3]>=avals[0] || avals[3]<=avals[1]))
  {
    //points are reverse lPoint is to the left, rPoint is to the right
    cl=rPoint+b;
    cr=lPoint+b;
  }
  calculateContourPointsTime+=(getmsofday()-atime);
}

void larvaFit::calculateContourPoints(Point2f &a,
                            Point2f &b,
                            Point2f &c,
                            Point2f &bp,
                            double b_index,
                            double width,
                            Point &cl,
                            Point &cr)
{
  unsigned long long atime=getmsofday();
  calculateContourPointsCalls++;
  //Make b our reference point
  Point2f r=a-b;
  Point2f f=c-b;
  Point2f sp=r+f;
  Point2f lPoint,rPoint;
  width=width_lookup[b_index]*width;
  if(fabs(sp.x)<0.0001 && fabs(sp.y)<0.0001) // On the same line
  {
    //We need the perpendicular vector
    if (r.x!=0 && r.y!=0)
    {
      double pslope=-r.x/r.y;
      double xc1=sqrt((width*width)/(1+pslope*pslope));
      lPoint.x=xc1;
      lPoint.y=pslope*xc1;
      rPoint.x=-xc1;
      rPoint.y=-pslope*xc1;
    }
    else if(r.x==0 && r.y==0)
    {
      cerr << "Weird case" << endl;
      exit(0);
    }
    else if(r.x==0)
    {
      lPoint.x=width;
      lPoint.y=0;
      rPoint.x=-width;
      rPoint.y=0;
    }
    else if(r.y==0)
    {
      lPoint.x=0;
      lPoint.y=width;
      rPoint.x=0;
      rPoint.y=-width;
    }
  
  }
  else
  {
    double ratio=width/sqrt(sp.x*sp.x+sp.y*sp.y);
    lPoint.x=ratio*sp.x;
    lPoint.y=ratio*sp.y;
    rPoint.x=-lPoint.x;
    rPoint.y=-lPoint.y;
  }

  double avals[4]={atan2f(r.y,r.x),
                    atan2f(-r.y,-r.x),
                    atan2f(rPoint.y,rPoint.x),
                    atan2f(lPoint.y,lPoint.x)};
  if(avals[0]<0)
    avals[0]=2*CV_PI+avals[0];
  if(avals[1]<0)
    avals[1]=2*CV_PI+avals[1];
  if(avals[2]<0)
    avals[2]=2*CV_PI+avals[2];
  if(avals[3]<0)
    avals[3]=2*CV_PI+avals[3];
  if((avals[2]>avals[0] || avals[2]<avals[1]) &&
     (avals[3]<avals[0] || avals[3]>avals[1]))
  {
    //points are correct lPoint is to the left, rPoint is to the right
    cl=lPoint+b-bp;
    cr=rPoint+b-bp;
  }
  else if((avals[2]<=avals[0] || avals[2]>=avals[1]) &&
          (avals[3]>=avals[0] || avals[3]<=avals[1]))
  {
    //points are reverse lPoint is to the left, rPoint is to the right
    cl=rPoint+b-bp;
    cr=lPoint+b-bp;
  }
  calculateContourPointsTime+=(getmsofday()-atime);
}

