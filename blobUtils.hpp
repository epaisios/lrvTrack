//#ifndef __LRVTRACK_BLOBUTILS_HPP__
//#define __LRVTRACK_BLOBUTILS_HPP__
#include <opencv2/core/core.hpp>
#include "cvblob.h"

float diff(cv::Point &a, cv::Point &b);
void blobToPointVector(cvb::CvBlob &p,std::vector<cv::Point> &points);
void createLarvaROI(cv::Mat &frame, cv::Mat &ROI, cvb::CvBlob &blob);

void createLarvaContour(cv::Mat &lrvROI,
                        cvb::CvBlob &blob,
                        int type=CV_8UC1);

void createLarvaContourPoints(cv::Mat &lrvROI,
                              cvb::CvBlob &blob,
                              int type=CV_8UC1);

float getGreyValue(cv::Mat &larvaROI, cvb::CvBlob &blob,cv::Mat &grey_frame);

float getPerimeter(cvb::CvBlob &blob);

float getSurroundingSize(cv::Point &point, cvb::CvBlob &blob,cv::Mat &grey_frame);
float plotAngle(cvb::CvBlob *blob,cv::Mat &ROIimg,int PAD=0);
float angle( cv::Point &pt1, cv::Point &pt0, cv::Point &pt2 );
//#endif
