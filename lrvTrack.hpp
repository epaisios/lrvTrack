#ifndef __LRVTRACK_ALL_HPP__
#define __LRVTRACK_ALL_HPP__
#include <map>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <fstream>
#include "larvaObject.hpp"
#include <boost/timer/timer.hpp>
#ifdef _WIN32
#include <winsock.h>
//#include "lrvTrackConfig.h"
#endif

using boost::timer::cpu_timer;
using boost::timer::cpu_times;
using boost::timer::nanosecond_type;

std::string VIDEO_TYPE=".avi";
//int VIDEO_CODEC=CV_FOURCC('X','2','6','4');
//int VIDEO_CODEC=CV_FOURCC('S','V','Q','3');
int VIDEO_CODEC=CV_FOURCC('F','M','P','4');
//int VIDEO_CODEC=CV_FOURCC('L','M','P','4');
//int VIDEO_CODEC=CV_FOURCC('M','P','4','2');
//int VIDEO_CODEC=CV_FOURCC('A','V','C','1');
//int VIDEO_CODEC=CV_FOURCC('M','P','G','4');
//int VIDEO_CODEC=CV_FOURCC('D','A','V','C');

std::ofstream summary;
int START_FRAME=0;

// *** FLAGS ***
int LRVTRACK_VERBOSE_LEVEL=1;
std::string LRVTRACK_RESULTS_FOLDER;
std::string LRVTRACK_DATE;
std::string LRVTRACK_NAME;
std::string LRVTRACK_SAVE_VIDEO;
std::string LRVTRACK_SAVE_PROCESSED_VIDEO;
std::string LRVTRACK_FILE_INPUT;
int  LRVTRACK_CAMERA_INPUT;
int  LRVTRACK_ODOUR_CUPS=0;
bool LRVTRACK_NORMALIZE=true;
bool LRVTRACK_SHOW_SKELETON=false;
bool LRVTRACK_SHOW_CONTOUR=false;
bool LRVTRACK_SHOW_ORIENTATION=false;
bool LRVTRACK_SHOW_CENTROID=true;
bool LRVTRACK_SHOW_HEAD_TAIL=true;
bool LRVTRACK_SHOW_TAGS=true;
// *** FLAGS ***

//VALUES FOR VARIOUS COMPARISONS
double LARVA_SIZE_COMPARISON_FACTOR=1.3;
double LARVA_CENTRE_COMPARISON_FACTOR=1.05;
double LARVA_MAHALANOBIS_THRESHOLD=11.0;
double IS_LARVA_THRESHOLD=1300;
double COLLISION_DURATION_THRESHOLD=48;
unsigned int HISTORY_SIZE=10;

double VIDEO_FPS=24.1;
static unsigned int CURRENT_FRAME=0;
unsigned int LARVAE_COUNT;

std::map<unsigned int, std::vector<unsigned int> > detected_clusters;
std::map<unsigned int,larvaObject> detected_larvae;
std::map<unsigned int, std::vector<unsigned int> > current_clusters;
std::map<unsigned int, std::vector<unsigned int> > current_diverged;

//********* Used by the new tracking algorithm ************************
// map of assignments of the blobs in the previous frame.
// Assignments are such that:
//    [ IDP of Previous Frame, [ IDN1, ..., IDNN]] IDs of
//    new Frame assigned to IDP
std::map<unsigned int, std::vector<unsigned int> > assignedPrevious;
std::map<unsigned int, unsigned int> assignedPreMap;

// map of assignments of the blobs in the current frame.
// Assignments are such that:
//    [ IDP of Current Frame, [ID_NN, IDN1, ..., IDNN]] IDs of
//    [ IDP of Current Frame, [ID_NN, IDN1, ..., IDNN]] IDs of
//    Previous Frame assigned to IDP
//    ID NN
std::map<unsigned int, std::vector<unsigned int> > assignedNew;
//**********************************************************************
std::vector<unsigned int> lost_larvae;
std::vector<unsigned int> current_new;
std::vector<unsigned int> current_gone;

cpu_timer tS,tP;
cpu_times FrameEllapsedTime;
//struct timeval tP;
//struct timeval tC;
//double FrameEllapsedTime;
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
