/* Skeletonization based on the algorithm of Zhang and Suen
 *
 * Implementation copied from here by ZachTM:
 * http://answers.opencv.org/question/3207/what-is-a-good-thinning-algorithm-for-getting-the/
 * It contains three functions:
 *  ThinSubiteration1, ThinSubiteration2 and larvae_skel (original: normalizeLetter)
 * */
#ifndef __LRVTRACK_BLOBUTILS_HPP__
#define __LRVTRACK_BLOBUTILS_HPP__

#include <opencv2/core/core.hpp>
#include <vector>

class larvaSkel
{
private:
  cv::Mat skelImg;
  std::vector<cv::Point> skelPoints;
  void ThinSubiteration1(cv::Mat & pSrc, cv::Mat & pDst);
  void ThinSubiteration2(cv::Mat & pSrc, cv::Mat & pDst);

public:
  std::vector <cv::Point> startPoints;
  cv::Point MidPoint;
  cv::Point Point20;
  cv::Point Point80;
  bool emptySkeleton;
  larvaSkel(bool emptySkel):emptySkeleton(true) {};
  larvaSkel(cv::Mat &inputarray, cv::Point &centroid);
  void drawSkeleton(cv::Mat &img,cv::Scalar col=cv::Scalar(0,0,0));
};
#endif
