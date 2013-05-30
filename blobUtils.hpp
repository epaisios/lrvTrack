//#ifndef __LRVTRACK_BLOBUTILS_HPP__
//#define __LRVTRACK_BLOBUTILS_HPP__
#include <opencv2/core/core.hpp>
#include "cvblob.h"

double diff(cv::Point &a, cv::Point &b);
void blobToPointVector(cvb::CvBlob &p,std::vector<cv::Point> &points);
void createLarvaROI(cv::Mat &frame, cv::Mat &ROI, cvb::CvBlob &blob);

void createLarvaContour(cv::Mat &lrvROI, 
                        cvb::CvBlob &blob,
                        int type=CV_8UC1);

void createLarvaContourPoints(cv::Mat &lrvROI, 
                              cvb::CvBlob &blob,
                              int type=CV_8UC1);

double getGreyValue(cv::Mat &larvaROI, cvb::CvBlob &blob,cv::Mat &grey_frame);

double getPerimeter(cvb::CvBlob &blob);

double getSurroundingSize(cv::Point &point, cvb::CvBlob &blob,cv::Mat &grey_frame);
double plotAngle(cvb::CvBlob *blob,cv::Mat &ROIimg,int PAD=0);
//#endif
