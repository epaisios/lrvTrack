#ifndef __LRVTRACK_LRVTRACKBASE_HPP__
#define __LRVTRACK_LRVTRACKBASE_HPP__
#include <utility>
#include <opencv2/core/core.hpp>

#ifdef LRVTRACK_WITH_CUDA
#include "opencv2/opencv.hpp"
#include "opencv2/gpu/gpu.hpp"
#endif

#define MIN_DISCARD_DISTANCE 30
#define ROI_PADDING 0
#define MAX_HORIZ_RESOLUTION 22000 //Pixels
#define ROUNDNESS_THRESHOLD 230000

typedef std::pair<cv::Point_<int>,cv::Point_<int> > PointPair;

#ifndef NO_SSE
#define ltsqrt SSESqrt_Recip_Times_X

// SSE Optimized reciprocal sqrt and mutliplies it by pIn to get the sqrt
// Optimized sqrt from http://stackoverflow.com/questions/1528727/why-is-sse-scalar-sqrtx-slower-than-rsqrtx-x
inline void SSESqrt_Recip_Times_X( float * pOut, float * pIn )
{
  __m128 in = _mm_load_ss( pIn );
  _mm_store_ss( pOut, _mm_mul_ss( in, _mm_rsqrt_ss( in ) ) );
  // compiles to movss, movaps, rsqrtss, mulss, movss
}
#else
#define ltsqrt mysqrt
inline void mysqrt(float * pOut, float * pIn)
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

#endif

