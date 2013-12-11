#include <lrvTrack.hpp>
#include <fstream>
#ifndef _WIN32
#include <sys/time.h>
#endif
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/timer/timer.hpp>
#include <iomanip>
#include <string>

#include "cvblob.h"
#include "lrvTrackBase.hpp"
#include "blobUtils.hpp"
#include "larvaDistanceMap.hpp"
#include "larvaSkel.hpp"

namespace po = boost::program_options;
namespace fs = boost::filesystem;
using namespace cv;
using namespace std;


//The function removes the background from the new frame
//and applies thresholding to derive the image in which we look
//for blobs.
// Input:
//    * newFrame: the new frame from the get_next_frame function
//    * greyBgFrame: a grey image of the background
//    * showFrame: the BGR image returned to the user
// Output:
//    * processedFrame: the frame containing only the greyscale shapes
//                      of the larvae
void process_frame(Mat &newFrame, //origFrame, 
                  Mat &greyBgFrame,
                  Mat &showFrame,
                  Mat &processedFrame, //outFrame, 
                  int mode)
{
    int THRESH=THRESH_BINARY;
    unsigned int CMAX=255;
    unsigned int NB=137;
    int THRESHOLD=-30;
    int METHOD=ADAPTIVE_THRESH_MEAN_C;
    Mat fgFrame; //Foreground frame
    Mat fgROI; // Foreground region of interest (only the part we want)
    Mat fgImage; //Final image where the foreground is grey and (
                  // hopefully everything else is black

    // fgFrame is created by the absolute difference of the
    // newFrame and the background. This gives us the best
    // flexibility also when the larva are over areas that we 
    // consider background.
    // TODO: Check if this is true:
    // The addWeighted would provide us with a fuzzier
    // removal which could also remove larvae. 
    absdiff(newFrame,greyBgFrame,fgFrame);
    //addWeighted(fgFrame, 1.0, test, -2.0, 0.0, test);

    fgROI=Mat(fgFrame.rows , fgFrame.cols,fgFrame.depth(),Scalar(0));

    //Use the registered bestCircle to get construct the ROI as a circle
    //where all things outside the circle (petri-dish) are black.
    if(circles.size()>0)
    {
      circle(fgROI, 
          Point2f(circles[bestCircle][0],circles[bestCircle][1]),
          int(circles[bestCircle][2]),Scalar(255),-1);

      // Add the circle as a marker to the showFrame
      circle(showFrame,
          Point2f(circles[bestCircle][0],circles[bestCircle][1]),
          int(circles[bestCircle][2]),
          Scalar(0,255,0),1);
    }
    else
    {
      // This is a case where there is no circle found
      fgROI=Mat(greyBgFrame.rows ,
          greyBgFrame.cols,
          greyBgFrame.depth(),
          Scalar(255));
    }

    // Same process for odor cups
    if(cups.size()>0)
    {
      for(unsigned int i=0; i<cups.size(); ++i)
      {
        circle(fgROI,
            Point2f(cups[i][0],cups[i][1]),
            int(cups[i][2]),
            Scalar(0),
            -1);
        circle(showFrame, 
            Point2f(cups[i][0],cups[i][1]),
            int(cups[i][2]),
            Scalar(0,0,255),
            1);
      }
    }

    // Construct a complete image of the BW fgROI and the whole frame
    fgImage=fgFrame&fgROI;
    // Apply a thresholding to extract those elements that escape.
    adaptiveThreshold(fgImage,
        processedFrame,
        CMAX,
        METHOD,
        THRESH,
        NB,
        THRESHOLD);
    Mat element = getStructuringElement(MORPH_CROSS, Size(3, 3));
    // Here we dilate since our thresholding effort will almost always
    // remove outer parts of the contour
    dilate(processedFrame,processedFrame,element);

    // The processedFrame is the outcome of the good image and filtered for
    // noise by the processedFrame
    processedFrame=fgImage&processedFrame;
}

// Function to extract and process each frame so that the background is dark
// and the forground white.
//   Input: 
//       * capture: the capture device
//   Output:
//       * output: the processed frame
//       * origFrame: the RGB version of the image

bool get_next_frame(VideoCapture &capture, Mat &output, Mat &colorOutput)
{
  bool retval=capture.read(output);
  if(output.channels()==3)
  {
    output.copyTo(colorOutput);
    cvtColor(output,output,CV_BGR2GRAY);
  }
  else
  {
    cvtColor(output,colorOutput,CV_GRAY2BGR);
  }
  if(retval==false)
    return retval;

  if(LRVTRACK_INVERT==true)
  {
    Mat ctout;
    double MF=1.3;
    unsigned int PF=0;
    int THRESHOLD=255;
    double cpow=1.25;
    output=output*MF+PF;
    addWeighted(output,1.0,output,-2.0,THRESHOLD,output);
    output.convertTo(output,CV_32F);
    pow(output,cpow,output);
    convertScaleAbs(output,output,1,0);
  }

  return retval;
}

//Function to extract the background
// The function retrieves a frame from the capture device and 
// returns the background into the bgFrame image.
//   Input:
//     * capture: the video capture device
//   Output:
//     * greyBgFrame: the Suggested background
void extractBackground(VideoCapture &capture,
                       Mat &greyBgFrame)
{
  Mat origFrame;
  // Grab a frame from the stream
  if(!get_next_frame(capture,greyBgFrame,origFrame))
  {
    //TODO: Error handling
    exit(1);
  }
  
  int votes=60; //For the HoughCircles call
  int THRESH=THRESH_BINARY;

  Mat ctout;
  unsigned int CMAX=255;
  unsigned int NB=275;
  int THRESHOLD=0;
  int METHOD=ADAPTIVE_THRESH_MEAN_C;
  // Thresholding here helps circle detection (trial and error)
  // TODO: Check if it's true.
  adaptiveThreshold(greyBgFrame,
      ctout,
      CMAX,
      METHOD,
      THRESH,
      NB,
      THRESHOLD);

  //Loop until we find a set of circles we can work with :)
  while (circles.size()==0 && votes >= 10)
  {
    HoughCircles(ctout, circles, CV_HOUGH_GRADIENT,
        4,   // accumulator resolution (size of the image / 2)
        2,  // minimum distance between two circles
        2, // Canny high threshold
        votes, // minimum number of votes
        greyBgFrame.rows/2.1, greyBgFrame.rows/1.95); // min and max radiusV
    votes-=10;
  }

  // Once we have a set of circles we try to get the one that will give
  // us the best ratio brigthness/size. Assign its ID to bestCircle
  double cutlim;
  cutlim=DBL_MAX;
  unsigned int mysz=circles.size();
  for(unsigned int i=0; i<circles.size();i++)
  {
    Mat cutout=Mat::zeros(ctout.rows,ctout.cols, ctout.type());
    circle(cutout,
               Point2f(circles[i][0],circles[i][1]),
               circles[i][2],
               Scalar(255),-1);
    cutout=cutout&greyBgFrame;
    double val=norm(cutout,NORM_L1)/
      (CV_PI*circles[i][2]*circles[i][2]);
    if(val<cutlim)
    {
      cutlim=val;
      bestCircle=i;
    }
  }
  // Find the odour cups:
  //  TODO: More testing on this needed.
  if (LRVTRACK_ODOUR_CUPS>0 || LRVTRACK_ODOUR_CUPS==-1)
    {
      Mat fgROI;
      fgROI=Mat::zeros(greyBgFrame.rows , 
                           greyBgFrame.cols,
                           greyBgFrame.depth());

      if(circles.size()>0)
        {
          circle(fgROI, 
                     Point2f(circles[bestCircle][0],circles[bestCircle][1]),
                     int(circles[bestCircle][2]),
                     Scalar(255),
                     -1);
        }

      Mat thr;
      thr=greyBgFrame&fgROI;

      morphologyEx(thr,thr,MORPH_OPEN,Mat(),Point2f(-1,-1),9);
      threshold(thr,
                    thr,
                    thresholdlow,
                    thresholdhigh,
                    THRESH_BINARY|THRESH_OTSU);
      HoughCircles(thr, cups, CV_HOUGH_GRADIENT,
                       2,   // accumulator resolution (size of the image / 2)
                       1,  // minimum distance between two circles
                       240, // Canny high threshold
                       50, // minimum number of votes
                       3, 40); // min and max radiusV
    }

  //Initialize the background frame. Everything that is in the circle
  //is black (i.e. will be excluded from the background). 
  //Outside the circle is white.
  if(circles.size()>0)
  {
    circle(greyBgFrame,
        Point2f(circles[bestCircle][0], circles[bestCircle][1]),
        circles[bestCircle][2],
        Scalar(),
        -1);
  }
  else
  {
    greyBgFrame=Mat(greyBgFrame.rows,
        greyBgFrame.cols,
        greyBgFrame.depth(),
        Scalar(0));
  }
  Mat tempFrame;
  process_frame(greyBgFrame,
               tempFrame,
               greyBgFrame,
               origFrame,
               cups,
               thresholdlow,
               thresholdhigh,
               THRESH_BINARY|THRESH_OTSU);

  Vec3f &mc=circles[bestCircle];
  double radius;
  findFurthestLarva(tempFrame,
                    Point2f(circles[bestCircle][0],circles[bestCircle][1]),
                    radius);
  if(radius>circles[bestCircle][2]*0.9)
    radius=circles[bestCircle][2]*0.9;
  if(bgFrame.channels()>1)
    {
      cvtColor(bgFrame,greyBgFrame,CV_BGR2GRAY);
    }
  else
    {
      bgFrame.copyTo(greyBgFrame);
    }
  if(circles.size()>0)
  {
    circle(greyBgFrame,
        Point2f(circles[bestCircle][0], circles[bestCircle][1]), // circle centre
        radius,       // circle radius
        Scalar(), // color
        -1);              // thickness
  }
  else
  {
    circle(greyBgFrame,
        Point2f(greyBgFrame.rows/2, greyBgFrame.cols/2), // circle centre
        radius,       // circle radius
        Scalar(), // color
        -1);              // thickness

  }
  Mat greyBgFrameDefault;
  greyBgFrame.copyTo(greyBgFrameDefault);
}

// Function to take care of the various input possibilities. 
// Set up the parameters in the capture device for direct camera 
// input/file input and correct FPS.
int setup_capture_input(VideoCapture &capture)
{
  // Check whether we have camera input or file
  if(LRVTRACK_CAMERA_INPUT != -2)
    {
      //capture.open(CV_CAP_DC1394);
      // These settings are for our setup.
      // TODO: Autodetection is easy for windows but still need to look into this.
      capture.open(300); 
      capture.set(CV_CAP_PROP_FRAME_WIDTH,1280); 
      capture.set(CV_CAP_PROP_FRAME_HEIGHT,1024);
      capture.set(CV_CAP_PROP_FPS,24);
    }
  else if (LRVTRACK_FILE_INPUT != "")
    {
      capture.open(LRVTRACK_FILE_INPUT);
    }
  capture.set(CV_CAP_PROP_FORMAT,CV_8U);

  if (capture.get(CV_CAP_PROP_FPS)==0)
    {
      VIDEO_FPS=24.1;
    }
  else
    {
      // Funny case where the FPS are mentioned to be double what they actually are.
      if ( capture.get(CV_CAP_PROP_FPS) > 30 )
        {
          VIDEO_FPS=capture.get(CV_CAP_PROP_FPS)/2;
          //cout << VIDEO_FPS << endl;
        }
      else
        {
          VIDEO_FPS=capture.get(CV_CAP_PROP_FPS);
          //cout << VIDEO_FPS << endl;
        }
    }

  if (!capture.isOpened())
    return -1;

  return 0;
}

// Function to handle command-line arguments
int handle_args(int argc, char* argv[])
{
  char cDate[16];
  time_t curtime;
  struct tm *loctime;

  /* Get the current time. */
  curtime = time (NULL);

  /* Convert it to local time representation. */
  loctime = localtime (&curtime);

  strftime(cDate,16,"%Y%m%d_%H%M%S",loctime);

  string DATE=string(cDate);
  LRVTRACK_DATE=DATE;
  string VIDEO_NAME=DATE + VIDEO_TYPE;
  string PROCESSED_VIDEO_NAME=DATE + "_PROC"+VIDEO_TYPE;

  try
    {
      // allowed only on command line
      po::options_description generic("Generic options");
      generic.add_options()
      ("version,V", "Print version string")
      ("help,h", "Produce help message")
      ("list-cameras,l","List available cameras.")
      ("verbose,v",
       po::value<int>(&LRVTRACK_VERBOSE_LEVEL),
       "Verbosity (0-3).")
      ;

      po::options_description IOopts("Input Output Options");

      IOopts.add_options()

      ("results-folder,r",
       po::value<string>(&LRVTRACK_RESULTS_FOLDER)->implicit_value("."),
       "Main folder where the results will be recorded. The experiment folder for each experiment (based on the date and time) will be created under this folder.(Default: Current folder)"
      )

      ("file-suffix,f",
       po::value<string>(&LRVTRACK_NAME)->implicit_value(""),
       "Suffix to use for the experiment output files under the experiment folder. The experiment folder name will be DATE-<file-suffix>.(Default: empty)"
      )

      ("save-video,s",
       po::value<string>
       (&LRVTRACK_SAVE_VIDEO)->implicit_value(DATE+VIDEO_TYPE),
       "Save original input as video. \
                         The option is disabled when the input is from a video \
                         source.(Default: <date-filesuffix>)"
      )

      ("save-tracked-video,p",
       po::value<string>
       (&LRVTRACK_SAVE_PROCESSED_VIDEO)->implicit_value(
         DATE+"_PROC"+VIDEO_TYPE),
       "Save resulting output as video with the filename \
                         specified.(Default: <date-filesuffix>_PROC)"
      )
      ("file-input,i",
       po::value<string>(&LRVTRACK_FILE_INPUT)->implicit_value(""),
       "Filename of video to be used as input for tracker."
      )

      ("camera-input,c",
       po::value<int>(&LRVTRACK_CAMERA_INPUT)->implicit_value(0),
       "Camera feed to be used as input for tracker. For now unfortunately \
                         one can only use numbers and no listing is possible."
      )

      ("csv-output,t",
       "Output will be csv file named as \"date_time@<suffix>.csv\". The file will be \
       saved under the experiment folder."
      )

      ("choreography-output,j",
       "Output summary and blob files to be processed by choreography."
      );

      po::options_description setupOpts("Setup Options");
      setupOpts.add_options()

      ("odour-cups,o",
       po::value<int> (&LRVTRACK_ODOUR_CUPS)->implicit_value(10),
       "When the flag is specified lrvTrack will look for odour cups. If specified\
                         with no arguments (e.g. -o) lrvTrack will look for as many as it can find.\
                         The number of odour cups can be specified as a value for this option."
      );

      po::options_description procOpts("Processing Options");
      procOpts.add_options()

      ("no-normalize,n",
       "Setting this flag will disable normalization of the brightness values of \
                         each frame. This degrades the accuracy of results but reduces \
                         cpu utilization. (Default: Not set)"
      );

      po::options_description displayOpts("Display Options");
      displayOpts.add_options()

      ("show-skeleton,S",
       po::value<bool> (&LRVTRACK_SHOW_SKELETON),
       "Show skeleton for detected larvae (Default: False)"
      )

      ("show-contour,C",
       po::value<bool> (&LRVTRACK_SHOW_CONTOUR),
       "Show contour for detected larvae (Default: False)"
      )

      ("show-orientation,O",
       po::value<bool> (&LRVTRACK_SHOW_ORIENTATION),
       "Show orientation for detected larvae (Default: False)"
      )

      ("show-centroid,Z",
       po::value<bool> (&LRVTRACK_SHOW_CENTROID),
       "Show centroid for detected larvae (Default: True)"
      )

      ("show-head-tail,H",
       po::value<bool> (&LRVTRACK_SHOW_HEAD_TAIL),
       "Show head and tail for detected larvae (Default: True)"
      )

      ("show-tags,T",
       po::value<bool> (&LRVTRACK_SHOW_TAGS),
       "Show larvae tags (Default: True)"
      )

      ("batch-testing,B",
       po::value<bool> (&LRVTRACK_SHOW_HEAD_TAIL),
       "Show no output except from the matchings. Used for testing internal parameters."
      );

      po::options_description cmdline_options;
      cmdline_options.add(generic).add(IOopts).add(setupOpts).add(procOpts).add(displayOpts);

      po::variables_map vm;
      po::store(po::parse_command_line(argc, argv, cmdline_options), vm);
      po::notify(vm);

      if (vm.count("help"))
        {
          cout << cmdline_options << "\n";
          exit(1);
        }
      if (vm.count("version"))
        {
          cout << "NO VERSION YET" << "\n";
          exit(1);
        }
      if (vm.count("results-folder"))
        {
          LRVTRACK_RESULTS_FOLDER=vm["results-folder"].as<string>();
        }
      else
        {
          LRVTRACK_RESULTS_FOLDER=".";
        }

      if(vm.count("csv-output"))
        {
          LRVTRACK_CSV_OUTPUT=true;
        }
      if(vm.count("choreography-output"))
        {
          LRVTRACK_CHOREOGRAPHY_OUTPUT=true;
        }
      if(vm.count("no-normalize"))
        {
          LRVTRACK_NORMALIZE=false;
        }
      if(vm.count("odour-cups"))
        {
          LRVTRACK_ODOUR_CUPS=vm["odour-cups"].as<int>();
        }

      if (vm.count("file-suffix"))
        {
          LRVTRACK_NAME=LRVTRACK_NAME;
        }
      else
        {
          LRVTRACK_NAME=DATE;
        }

      if (vm.count("save-video"))
        {
          LRVTRACK_SAVE_VIDEO=LRVTRACK_NAME+VIDEO_TYPE;
        }
      else
        {
          LRVTRACK_SAVE_VIDEO="";
        }

      if (vm.count("save-tracked-video"))
        {
          LRVTRACK_SAVE_PROCESSED_VIDEO=LRVTRACK_NAME+"_PROC"+VIDEO_TYPE;
        }
      else
        {
          LRVTRACK_SAVE_PROCESSED_VIDEO="";
        }

      if (vm.count("list-cameras"))
        {
          cerr << "Camera listing not yet implemented." <<  endl;
          exit(1);
        }

      if (vm.count("file-input")<=0 && vm.count("camera-input")<=0)
        {
          cerr << "Error: No Input Specified. Exiting..." <<  endl;
          exit(2);
        }
      else if (vm.count("file-input")>0 && vm.count("camera-input")>0)
        {
          cerr << "Error: Ambiguous Input Specified. Exiting..." <<  endl;
          exit(2);
        }

      if (vm.count("file-input")>0)
        {
          LRVTRACK_FILE_INPUT=vm["file-input"].as<string>();
          if (LRVTRACK_FILE_INPUT=="")
            {
              cerr << "Error: Input file flag given but no file specified. Exiting..." << endl;
              exit(3);
            }
        }

      if (vm.count("camera-input")>0)
        {
          LRVTRACK_CAMERA_INPUT=vm["camera-input"].as<int>();
        }
      else
        {
          LRVTRACK_CAMERA_INPUT=-2;
        }
    }
  catch(po::error &e)
    {
      cerr << "Problem parsing options: " << e.what() << endl;
      exit(1);
    }

  return 0;
}

int main(int argc, char* argv[])
{
  bool SHOWTAGS=true;
  bool TRACK=false;
  bool STEP=true;

#ifdef LRVTRACK_WITH_CUDA
  gpu::printShortCudaDeviceInfo(gpu::getDevice());
#elif defined(LRVTRACK_WITH_OPENCL)
  vector<ocl::Info> oclinfo;
  ocl::getDevice(oclinfo);
  cout << "Found " << oclinfo.size() << " OpenCL devices. Using: ";
  cout << oclinfo[0].DeviceName[0] << endl;
#endif
  // Initialize video capture device
  VideoCapture capture;
  // Handle command line arguments
  handle_args(argc,argv);
  
  if(setup_capture_input(capture) == -1)
  {
    cerr << "Error setting up the capture device (camera or video file)" << endl;
    exit(1);
  }

}
