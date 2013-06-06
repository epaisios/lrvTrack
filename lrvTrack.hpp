#ifndef __LRVTRACK_ALL_HPP__
#define __LRVTRACK_ALL_HPP__
#include <map>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include "larvaObject.hpp"

static double VIDEO_FPS=24.0;
static unsigned int CURRENT_FRAME=0;
unsigned int LARVAE_COUNT;

std::map<unsigned int, std::vector<unsigned int> > detected_clusters;
std::map<unsigned int,larvaObject> detected_larvae;
std::map<unsigned int, std::vector<unsigned int> > current_clusters;
std::map<unsigned int, std::vector<unsigned int> > current_diverged;
std::vector<unsigned int> current_new;
std::vector<unsigned int> current_gone;

struct timeval tS;
struct timeval tP;
struct timeval tC;
double FrameEllapsedTime;
cv::Mat frame;
cv::Mat thresholded_frame;
cv::Mat grey_frame;
cv::Mat bgFrame;
cv::Mat previousFrame;
IplImage *labelImg;
int DEBUG_INFO=0;

double Wlength=1.0;
double Wsize=0.2;

#endif
