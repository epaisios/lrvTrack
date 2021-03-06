OPENCV_VER=2.4.7
INCLUDES=-I/opt/opencv/${OPENCV_VER}/include/opencv -I/opt/opencv/${OPENCV_VER}/include -I/usr/local/include -I.
INCLUDES_GPU=-I/opt/opencv/${OPENCV_VER}/include/opencv -I/opt/opencv/${OPENCV_VER}/include -I/usr/local/include -I.
INCLUDES_QT=-I/opt/opencv/${OPENCV_VER}-qt/include/opencv -I/opt/opencv/${OPENCV_VER}-qt/include -I/usr/local/include -I.
#GCCFLAGS=-fopenmp -DAE_CPU=AE_INTEL
GCCFLAGS=-DAE_CPU=AE_INTEL -std=c++11 -Wno-c++11-extensions
GCCFLAGS_GPU=$(GCCFLAGS) -DLRVTRACK_WITH_CUDA
GCCFLAGS_OPENCL=$(GCCFLAGS) -DLRVTRACK_WITH_OPENCL
OPTFLAGS=-Ofast -O4
LIBS=-lm\
		  -lopencv_contrib\
		 	-lopencv_core\
		 	-lopencv_highgui\
		 	-lopencv_video\
		 	-lopencv_imgproc\
			-lopencv_ml\
		 	-lcvblob\
			-ltbb\
		 	-lboost_program_options-mt\
		 	-lboost_filesystem-mt\
		 	-lboost_system-mt\
			-lboost_timer-mt
all:
	@g++ $(GCCFLAGS) -g -O0 -ansi $(INCLUDES) -c alglib/alglibinternal.cpp -o alglib/alglibinternal.o
	@g++ $(GCCFLAGS) -g -O0 -ansi $(INCLUDES) -c alglib/alglibmisc.cpp -o alglib/alglibmisc.o
	@g++ $(GCCFLAGS) -g -O0 -ansi $(INCLUDES) -c alglib/ap.cpp -o alglib/ap.o
	@g++ $(GCCFLAGS) -g -O0 -ansi $(INCLUDES) -c alglib/integration.cpp -o alglib/integration.o
	@g++ $(GCCFLAGS) -g -O0 -ansi $(INCLUDES) -c alglib/interpolation.cpp -o alglib/interpolation.o
	@g++ $(GCCFLAGS) -g -O0 -ansi $(INCLUDES) -c alglib/linalg.cpp -o alglib/linalg.o
	@g++ $(GCCFLAGS) -g -O0 -ansi $(INCLUDES) -c alglib/optimization.cpp -o alglib/optimization.o
	@g++ $(GCCFLAGS) -g -O0 -ansi $(INCLUDES) -c alglib/solvers.cpp -o alglib/solvers.o
	@g++ $(GCCFLAGS) -g -O0 -ansi $(INCLUDES) -c alglib/specialfunctions.cpp -o alglib/specialfunctions.o
	@g++ $(GCCFLAGS) -g -O0 -ansi $(INCLUDES) -c alglib/statistics.cpp -o alglib/statistics.o
	@g++ $(GCCFLAGS) -g -O0 -ansi $(INCLUDES) -c blobUtils.cpp -o blobUtils.o
	@g++ $(GCCFLAGS) -g -O0 -ansi $(INCLUDES) -c larvaSkel.cpp -o larvaSkel.o
	@g++ $(GCCFLAGS) -g -O0 -ansi $(INCLUDES) -c larvaDistanceMap.cpp -o larvaDistanceMap.o
	@g++ $(GCCFLAGS) -g -O0 -ansi $(INCLUDES) -c larvaObject.cpp -o larvaObject.o
	@g++ $(GCCFLAGS) -g -O0 -ansi $(INCLUDES) -L/opt/opencv/${OPENCV_VER}/lib -L/usr/local/lib main.cpp $(LIBS) -o lrvTrack *.o alglib/*.o

gpu:
	g++ $(GCCFLAGS_GPU) -g -Wall -ansi $(INCLUDES) -c alglib/alglibinternal.cpp -o alglibinternal.o
	g++ $(GCCFLAGS_GPU) -g -Wall -ansi $(INCLUDES) -c alglib/interpolation.cpp -o interpolation.o
	g++ $(GCCFLAGS_GPU) -g $(INCLUDES_GPU) -L/opt/opencv/${OPENCV_VER}/lib -L/usr/local/lib -c blobUtils.cpp -o blobUtils.o
	g++ $(GCCFLAGS_GPU) -g $(INCLUDES_GPU) -L/opt/opencv/${OPENCV_VER}/lib -L/usr/local/lib -c larvaSkel.cpp -o larvaSkel.o
	g++ $(GCCFLAGS_GPU) -g $(INCLUDES_GPU) -L/opt/opencv/${OPENCV_VER}/lib -L/usr/local/lib -c larvaDistanceMap.cpp -o larvaDistanceMap.o
	g++ $(GCCFLAGS_GPU) -g $(INCLUDES_GPU) -L/opt/opencv/${OPENCV_VER}/lib -L/usr/local/lib -c larvaObject.cpp -o larvaObject.o
	g++ $(GCCFLAGS_GPU) -g $(INCLUDES_GPU) -framework CUDA -L/opt/opencv/${OPENCV_VER}/lib -L/usr/local/lib main.cpp $(LIBS) -lopencv_gpu  -o lrvTrack *.o

ocl:
	g++ $(GCCFLAGS_OPENCL) -g $(INCLUDES_GPU) -L/opt/opencv/${OPENCV_VER}/lib -L/usr/local/lib -c blobUtils.cpp -o blobUtils.o
	g++ $(GCCFLAGS_OPENCL) -g $(INCLUDES_GPU) -L/opt/opencv/${OPENCV_VER}/lib -L/usr/local/lib -c larvaSkel.cpp -o larvaSkel.o
	g++ $(GCCFLAGS_OPENCL) -g $(INCLUDES_GPU) -L/opt/opencv/${OPENCV_VER}/lib -L/usr/local/lib -c larvaDistanceMap.cpp -o larvaDistanceMap.o
	g++ $(GCCFLAGS_OPENCL) -g $(INCLUDES_GPU) -L/opt/opencv/${OPENCV_VER}/lib -L/usr/local/lib -c larvaObject.cpp -o larvaObject.o
	g++ $(GCCFLAGS_OPENCL) -g $(INCLUDES_GPU) -framework OpenCL -L/opt/opencv/${OPENCV_VER}/lib -L/usr/local/lib main.cpp $(LIBS) -lopencv_gpu -lopencv_ocl -o lrvTrack *.o

clean: 
	@rm -rf alglib/*.o *.o lrvTrack *.dSYM

qt:
	@g++-4.9 $(GCCFLAGS) -g $(INCLUDES_QT) -F/usr/local/lib -framework QtCore -framework QtGui -DQT_GUI_LIB -DQT_CORE_LIB -DUNICODE -L/opt/opencv/${OPENCV_VER}-qt/lib -L/usr/local/lib -c blobUtils.cpp -o blobUtils.o
	@g++-4.9 $(GCCFLAGS) -g $(INCLUDES_QT) -F/usr/local/lib -framework QtCore -framework QtGui -DQT_GUI_LIB -DQT_CORE_LIB -DUNICODE -L/opt/opencv/${OPENCV_VER}-qt/lib -L/usr/local/lib -c larvaSkel.cpp -o larvaSkel.o
	@g++-4.9 $(GCCFLAGS) -g $(INCLUDES_QT) -F/usr/local/lib -framework QtCore -framework QtGui -DQT_GUI_LIB -DQT_CORE_LIB -DUNICODE -L/opt/opencv/${OPENCV_VER}-qt/lib -L/usr/local/lib -c larvaDistanceMap.cpp -o larvaDistanceMap.o
	@g++-4.9 $(GCCFLAGS) -g $(INCLUDES_QT) -F/usr/local/lib -framework QtCore -framework QtGui  -DQT_GUI_LIB -DQT_CORE_LIB -DUNICODE -L/opt/opencv/${OPENCV_VER}-qt/lib -L/usr/local/lib -c larvaObject.cpp -o larvaObject.o
	@g++-4.9 $(GCCFLAGS) -g $(INCLUDES_QT) -F/usr/local/lib -framework QtCore -framework QtGui  -DQT_GUI_LIB -DQT_CORE_LIB -DUNICODE -L/opt/opencv/${OPENCV_VER}-qt/lib -L/usr/local/lib main.cpp -lm -lopencv_contrib -lopencv_core -lopencv_highgui -lopencv_video -lopencv_imgproc -lcvblob -o lrvTrack *.o

opt:
	@g++ $(GCCFLAGS) ${OPTFLAGS} -ansi $(INCLUDES) -c alglib/alglibinternal.cpp -o alglib/alglibinternal.o
	@g++ $(GCCFLAGS) ${OPTFLAGS} -ansi $(INCLUDES) -c alglib/alglibmisc.cpp -o alglib/alglibmisc.o
	@g++ $(GCCFLAGS) ${OPTFLAGS} -ansi $(INCLUDES) -c alglib/ap.cpp -o alglib/ap.o
	@g++ $(GCCFLAGS) ${OPTFLAGS} -ansi $(INCLUDES) -c alglib/integration.cpp -o alglib/integration.o
	@g++ $(GCCFLAGS) ${OPTFLAGS} -ansi $(INCLUDES) -c alglib/interpolation.cpp -o alglib/interpolation.o
	@g++ $(GCCFLAGS) ${OPTFLAGS} -ansi $(INCLUDES) -c alglib/linalg.cpp -o alglib/linalg.o
	@g++ $(GCCFLAGS) ${OPTFLAGS} -ansi $(INCLUDES) -c alglib/optimization.cpp -o alglib/optimization.o
	@g++ $(GCCFLAGS) ${OPTFLAGS} -ansi $(INCLUDES) -c alglib/solvers.cpp -o alglib/solvers.o
	@g++ $(GCCFLAGS) ${OPTFLAGS} -ansi $(INCLUDES) -c alglib/specialfunctions.cpp -o alglib/specialfunctions.o
	@g++ $(GCCFLAGS) ${OPTFLAGS} -ansi $(INCLUDES) -c alglib/statistics.cpp -o alglib/statistics.o
	@g++ $(GCCFLAGS) ${OPTFLAGS} -ansi $(INCLUDES) -c blobUtils.cpp -o blobUtils.o
	@g++ $(GCCFLAGS) ${OPTFLAGS} -ansi $(INCLUDES) -c larvaSkel.cpp -o larvaSkel.o
	@g++ $(GCCFLAGS) ${OPTFLAGS} -ansi $(INCLUDES) -c larvaDistanceMap.cpp -o larvaDistanceMap.o
	@g++ $(GCCFLAGS) ${OPTFLAGS} -ansi $(INCLUDES) -c larvaObject.cpp -o larvaObject.o
	@g++ $(GCCFLAGS) ${OPTFLAGS} -ansi $(INCLUDES) -L/opt/opencv/${OPENCV_VER}/lib -L/usr/local/lib main.cpp $(LIBS) -o lrvTrack *.o alglib/*.o
