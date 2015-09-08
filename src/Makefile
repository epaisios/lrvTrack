OPENCV_VER=3.0.0
#INCLUDES=-I/opt/opencv/${OPENCV_VER}/include/opencv -I/opt/opencv/${OPENCV_VER}/include -I/usr/local/Cellar/swi-prolog/HEAD/libexec/lib/swipl-7.1.28/include -I.
INCLUDES=-I/opt/opencv/${OPENCV_VER}/include/opencv -I/opt/opencv/${OPENCV_VER}/include -I.
#INCLUDES_GPU=-I/opt/opencv/${OPENCV_VER}/include/opencv -I/opt/opencv/${OPENCV_VER}/include -I/usr/local/Cellar/swi-prolog/HEAD/libexec/lib/swipl-7.1.28/include -I.
INCLUDES_GPU=-I/opt/opencv/${OPENCV_VER}/include/opencv -I/opt/opencv/${OPENCV_VER}/include -I.
INCLUDES_QT=-I/opt/opencv/${OPENCV_VER}-qt/include/opencv -I/opt/opencv/${OPENCV_VER}-qt/include -I.
#GCCFLAGS=-fopenmp -DAE_CPU=AE_INTEL
#GCCFLAGS=-DAE_CPU=AE_INTEL -std=c++11 -Wno-c++11-extensions
GCCFLAGS=-std=c++0x -Wno-c++11-extensions -DBOOST_LOG_DYN_LINK
GCCFLAGS_GPU=-DLRVTRACK_WITH_CUDA
GCCFLAGS_OCL=-DLRVTRACK_WITH_OPENCL
OPTFLAGS=-O3 -msse2 -DAE_CPU=AE_INTEL -march=native 
DEBUGFLAGS=-g -O0
OPENCL_LIBS=-lopencv_ocl
     #-lopencv_contrib#
LIBS=-lm\
     -lopencv_core\
     -lopencv_highgui\
     -lopencv_videoio\
     -lopencv_video\
     -lopencv_imgcodecs\
     -lopencv_imgproc\
     -lopencv_photo\
     -lopencv_ml\
     -lcvblob\
     -ltbb\
     -llpsolve55\
     -lboost_log-mt\
     -lboost_log_setup-mt\
     -lboost_program_options-mt\
     -lboost_filesystem-mt\
     -lboost_system-mt\
     -lboost_timer-mt
     #-lswipl

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
 lpsolver.o\
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
	#g++ $(GCCFLAGS) $(DEBUGFLAGS) $(INCLUDES) -L/usr/local/Cellar/swi-prolog/HEAD/libexec/lib/swipl-7.1.28/lib/x86_64-darwin13.4.0 -L/opt/opencv/${OPENCV_VER}/lib -L/usr/local/lib $(LIBS) -o $@ $^

all: lrvTrack

opt: GCCFLAGS+=$(OPTFLAGS) 
opt: DEBUGFLAGS= 
opt: all

gpu_debug: GCCFLAGS+=$(GCCFLAGS_GPU) 
gpu_debug: LIBS+=-lopencv_gpu
gpu_debug: all

ocl_debug: GCCFLAGS+=$(GCCFLAGS_OCL) 
#ocl_debug: LIBS+=-lopencv_ocl
ocl_debug: all
	
gpu: GCCFLAGS+=$(GCCFLAGS_GPU) 
gpu: LIBS+=-lopencv_gpu
gpu: opt

ocl: GCCFLAGS+=$(GCCFLAGS_OCL) 
#ocl: LIBS+=-lopencv_ocl
ocl: opt

clean: 
	@rm -rf alglib/*.o *.o lrvTrack *.dSYM *.clb

