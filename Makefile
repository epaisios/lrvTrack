OPENCV_VER=2.4.7
INCLUDES=-I/opt/opencv/${OPENCV_VER}/include/opencv -I/opt/opencv/${OPENCV_VER}/include -I.
INCLUDES_GPU=-I/opt/opencv/${OPENCV_VER}/include/opencv -I/opt/opencv/${OPENCV_VER}/include -I.
INCLUDES_QT=-I/opt/opencv/${OPENCV_VER}-qt/include/opencv -I/opt/opencv/${OPENCV_VER}-qt/include -I.
#GCCFLAGS=-fopenmp -DAE_CPU=AE_INTEL
#GCCFLAGS=-DAE_CPU=AE_INTEL -std=c++11 -Wno-c++11-extensions
GCCFLAGS=-std=c++0x -pipe -march=core-avx-i -Wno-c++11-extensions -DBOOST_LOG_DYN_LINK
GCCFLAGS_GPU=$(GCCFLAGS) -DLRVTRACK_WITH_CUDA
GCCFLAGS_OPENCL=$(GCCFLAGS) -DLRVTRACK_WITH_OPENCL
OPTFLAGS=-O2
DEBUGFLAGS=-g
OPENCL_LIBS=-lopencv_ocl
LIBS=-lm\
		  -lopencv_contrib\
		 	-lopencv_core\
		 	-lopencv_highgui\
		 	-lopencv_video\
		 	-lopencv_imgproc\
			-lopencv_photo\
			-lopencv_ml\
		 	-lcvblob\
			-ltbb\
		 	-lboost_log-mt\
			-lboost_log_setup-mt\
		 	-lboost_program_options-mt\
		 	-lboost_filesystem-mt\
		 	-lboost_system-mt\
			-lboost_timer-mt

ALG_OBJ=alglib/alglibinternal.o\
 alglib/alglibmisc.o\
 alglib/ap.o\
 alglib/integration.o\
 alglib/interpolation.o\
 alglib/linalg.o\
 alglib/optimization.o\
 alglib/solvers.o\
 alglib/specialfunctions.o\
 alglib/statistics.o

LT_OBJ=lrvTrackBase.o\
 blobUtils.o\
 larvaDistanceMap.o\
 lrvTrackDebug.o\
 larvaObject.o\
 lrvTrackFit.o\
 lrvTrackOL.o

LT_DEPS=lrvTrackBase.hpp\
blobUtils.hpp\
larvaDistanceMap.hpp\
lrvTrackDebug.hpp\
larvaObject.hpp\
lrvTrackFit.hpp\
lrvTrackPartitionGenerator.hpp\
lrvTrack.hpp

alglib/%.o: alglib/%.cpp
	g++ $(GCCFLAGS) $(DEBUGFLAGS) -o $@ -c $?
	
%.o: %.cpp $(LT_DEPS)
	g++ $(GCCFLAGS) $(DEBUGFLAGS) $(INCLUDES) -c $<

lrvTrack: $(ALG_OBJ) $(LT_OBJ)
	g++ $(GCCFLAGS) $(DEBUGFLAGS) $(INCLUDES) -L/opt/opencv/${OPENCV_VER}/lib -L/usr/local/lib $(LIBS) -o $@ $^

all: lrvTrack

opt: GCCFLAGS+=$(OPTFLAGS) 
opt: DEBUGFLAGS= 
opt: all
	

clean: 
	@rm -rf alglib/*.o *.o lrvTrack *.dSYM

