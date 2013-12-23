#ifndef __LRVTRACK_LRVTRACKBASE_HPP__
#define __LRVTRACK_LRVTRACKBASE_HPP__
#include <utility>
#include <opencv2/core/core.hpp>
#include <opencv2/opencv.hpp>
#include <opencv2/imgproc/imgproc.hpp>

#ifdef LRVTRACK_WITH_CUDA
#include "opencv2/gpu/gpu.hpp"
#endif

#ifdef LRVTRACK_WITH_OPENCL
#include "opencv2/ocl/ocl.hpp"
#endif

#define MIN_DISCARD_DISTANCE 30
#define ROI_PADDING 0
#define MAX_HORIZ_RESOLUTION 22000 //Pixels
#define ROUNDNESS_THRESHOLD 230000

typedef std::pair<cv::Point2f,cv::Point2f > PointPair;
#ifndef NO_SSE
#define ltsqrt SSESqrt_Recip_Times_X

// SSE Optimized reciprocal sqrt and mutliplies it by pIn to get the sqrt
// Optimized sqrt from http://stackoverflow.com/questions/1528727/why-is-sse-scalar-sqrtx-slower-than-rsqrtx-x
inline void SSESqrt_Recip_Times_X( float* pOut, float* pIn )
{
  __m128 in = _mm_load_ss( pIn );
  _mm_store_ss( pOut, _mm_mul_ss( in, _mm_rsqrt_ss( in ) ) );
  // compiles to movss, movaps, rsqrtss, mulss, movss
}
#else
#define ltsqrt mysqrt
inline void mysqrt(float* pOut, float* pIn)
{
  *pOut=sqrt(*pIn);
}
#endif

#ifdef LRVTRACK_WITH_CUDA
inline void lrvTrackNormalize(cv::Mat &src,
                              cv::Mat &dst,
                              double alpha,
                              double beta,
                              int norm_type)
{
  cv::gpu::GpuMat dst_host, src_host;
  src_host.upload(src);
  cv::gpu::normalize(src_host,dst_host,alpha,beta,norm_type);
  dst_host.download(dst);
}

#elif defined(LRVTRACK_WITH_OPENCL)
inline void lrvTrackNormalize(cv::Mat &_src , cv::Mat &_dst, double a, double b,
                              int norm_type)
{
  cv::ocl::oclMat src(_src);
  std::cout << "Using OPENCL normalization" << std::endl;
  int rtype;
  double scale = 1, shift = 0;
  if( norm_type == CV_MINMAX )
    {
      double smin = 0, smax = 0;
      double dmin = MIN( a, b ), dmax = MAX( a, b );
      cv::ocl::minMaxLoc( src, &smin, &smax, 0, 0);
      scale = (dmax - dmin)*(smax - smin > DBL_EPSILON ? 1./(smax - smin) : 0);
      shift = dmin - smin*scale;
    }
  else if( norm_type == CV_L2 || norm_type == CV_L1 || norm_type == CV_C )
    {
      scale = cv::ocl::norm( src, norm_type );
      scale = scale > DBL_EPSILON ? a/scale : 0.;
      shift = 0;
    }
  else
    CV_Error( CV_StsBadArg, "Unknown/unsupported norm type" );

  rtype =  src.depth();

  _dst.create(_src.dims , _src.size, CV_MAKETYPE(rtype, _src.channels()));
  cv::ocl::oclMat dst(_dst);

  src.convertTo( dst, rtype, scale, shift );
  dst.upload(_dst);
}

#else
inline void lrvTrackNormalize(cv::InputArray src,
                              cv::OutputArray dst,
                              double alpha,
                              double beta,
                              int norm_type)
{
  cv::normalize(src,dst,alpha,beta,norm_type);
}
#endif

namespace std
{
  template <typename number >
    std::string printVector(std::vector<number> vec,int position=0)
    {
      if (vec.size()==0)
        return "";
      std::stringstream sstm;
      //bool const is_number= std::is_arithmetic<number>::value;
      //static_assert( is_number, "Provided type is not an arithmetic type");
      sstm << "[" ;
      typename std::vector<number>::const_iterator i=vec.begin()+position;
      sstm << *i ;
      ++i;
      for( ; i != vec.end(); ++i)
        sstm << ","<< *i;
      sstm << "]";
      return sstm.str();
    }
}

#endif

