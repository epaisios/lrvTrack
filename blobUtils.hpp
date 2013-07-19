//#ifndef __LRVTRACK_BLOBUTILS_HPP__
//#define __LRVTRACK_BLOBUTILS_HPP__
#include <opencv2/core/core.hpp>
#include "cvblob.h"

double diff(cv::Point2f &a, cv::Point2f &b);
void blobToPointVector(cvb::CvBlob &p,std::vector<cv::Point2f> &points);
void createLarvaROI(cv::Mat &frame, cv::Mat &ROI, cvb::CvBlob &blob);

void createLarvaContour(cv::Mat &lrvROI,
                        cvb::CvBlob &blob,
                        int type=CV_8UC1,
                        int PADDING=0,
                        bool FILL=true);

void createLarvaContourPoints(cv::Mat &lrvROI,
                              cvb::CvBlob &blob,
                              int type=CV_8UC1,
                              int PADDING=0);

double getGreyValue(cv::Mat &larvaROI, cvb::CvBlob &blob,cv::Mat &grey_frame);

double getPerimeter(cvb::CvBlob &blob);

double getSurroundingSize(cv::Point2f &point, cvb::CvBlob &blob,cv::Mat &grey_frame);
double plotAngle(cvb::CvBlob *blob,cv::Mat &ROIimg,int PAD=0);
double angle( cv::Point2f &pt1, cv::Point2f &pt0, cv::Point2f &pt2 );
double angleD( cv::Point2f &pt1, cv::Point2f &pt0, cv::Point2f &pt2 );
//#endif
