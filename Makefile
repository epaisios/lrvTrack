INCLUDES=-I/opt/opencv/2.4.6.1-gpu/include/opencv -I/opt/opencv/2.4.6.1-gpu/include -I/usr/local/include -I.
INCLUDES_GPU=-I/opt/opencv/2.4.6.1-gpu/include/opencv -I/opt/opencv/2.4.6.1-gpu/include -I/usr/local/include -I.
INCLUDES_QT=-I/opt/opencv/2.4.5-qt/include/opencv -I/opt/opencv/2.4.5-qt/include -I/usr/local/include -I.
GCCFLAGS=-fopenmp
GCCFLAGS_GPU=-fopenmp -DLRVTRACK_WITH_CUDA
GCCFLAGS_OPENCL=-fopenmp -DLRVTRACK_WITH_OPENCL
OPTFLAGS=-Ofast
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
	@g++ $(GCCFLAGS) -g -Wall -ansi $(INCLUDES) -L/opt/opencv/2.4.6.1-gpu/lib -L/usr/local/lib -c blobUtils.cpp -o blobUtils.o
	@g++ $(GCCFLAGS) -g -Wall -ansi $(INCLUDES) -L/opt/opencv/2.4.6.1-gpu/lib -L/usr/local/lib -c larvaSkel.cpp -o larvaSkel.o
	@g++ $(GCCFLAGS) -g -Wall -ansi $(INCLUDES) -L/opt/opencv/2.4.6.1-gpu/lib -L/usr/local/lib -c larvaDistanceMap.cpp -o larvaDistanceMap.o
	@g++ $(GCCFLAGS) -g -Wall -ansi $(INCLUDES) -L/opt/opencv/2.4.6.1-gpu/lib -L/usr/local/lib -c larvaObject.cpp -o larvaObject.o
	@g++ $(GCCFLAGS) -g -Wall -ansi $(INCLUDES) -L/opt/opencv/2.4.6.1-gpu/lib -L/usr/local/lib main.cpp $(LIBS) -o lrvTrack *.o

gpu:
	g++ $(GCCFLAGS_GPU) -g $(INCLUDES_GPU) -L/opt/opencv/2.4.6.1-gpu/lib -L/usr/local/lib -c blobUtils.cpp -o blobUtils.o
	g++ $(GCCFLAGS_GPU) -g $(INCLUDES_GPU) -L/opt/opencv/2.4.6.1-gpu/lib -L/usr/local/lib -c larvaSkel.cpp -o larvaSkel.o
	g++ $(GCCFLAGS_GPU) -g $(INCLUDES_GPU) -L/opt/opencv/2.4.6.1-gpu/lib -L/usr/local/lib -c larvaDistanceMap.cpp -o larvaDistanceMap.o
	g++ $(GCCFLAGS_GPU) -g $(INCLUDES_GPU) -L/opt/opencv/2.4.6.1-gpu/lib -L/usr/local/lib -c larvaObject.cpp -o larvaObject.o
	g++ $(GCCFLAGS_GPU) -g $(INCLUDES_GPU) -framework CUDA -L/opt/opencv/2.4.6.1-gpu/lib -L/usr/local/lib main.cpp $(LIBS) -lopencv_gpu  -o lrvTrack *.o

ocl:
	g++ $(GCCFLAGS_OPENCL) -g $(INCLUDES_GPU) -L/opt/opencv/2.4.6.1-gpu/lib -L/usr/local/lib -c blobUtils.cpp -o blobUtils.o
	g++ $(GCCFLAGS_OPENCL) -g $(INCLUDES_GPU) -L/opt/opencv/2.4.6.1-gpu/lib -L/usr/local/lib -c larvaSkel.cpp -o larvaSkel.o
	g++ $(GCCFLAGS_OPENCL) -g $(INCLUDES_GPU) -L/opt/opencv/2.4.6.1-gpu/lib -L/usr/local/lib -c larvaDistanceMap.cpp -o larvaDistanceMap.o
	g++ $(GCCFLAGS_OPENCL) -g $(INCLUDES_GPU) -L/opt/opencv/2.4.6.1-gpu/lib -L/usr/local/lib -c larvaObject.cpp -o larvaObject.o
	g++ $(GCCFLAGS_OPENCL) -g $(INCLUDES_GPU) -framework OpenCL -L/opt/opencv/2.4.6.1-gpu/lib -L/usr/local/lib main.cpp $(LIBS) -lopencv_gpu -lopencv_ocl -o lrvTrack *.o

clean: 
	@rm -rf *.o lrvTrack *.dSYM

qt:
	@g++-4.9 $(GCCFLAGS) -g $(INCLUDES_QT) -F/usr/local/lib -framework QtCore -framework QtGui -DQT_GUI_LIB -DQT_CORE_LIB -DUNICODE -L/opt/opencv/2.4.5-qt/lib -L/usr/local/lib -c blobUtils.cpp -o blobUtils.o
	@g++-4.9 $(GCCFLAGS) -g $(INCLUDES_QT) -F/usr/local/lib -framework QtCore -framework QtGui -DQT_GUI_LIB -DQT_CORE_LIB -DUNICODE -L/opt/opencv/2.4.5-qt/lib -L/usr/local/lib -c larvaSkel.cpp -o larvaSkel.o
	@g++-4.9 $(GCCFLAGS) -g $(INCLUDES_QT) -F/usr/local/lib -framework QtCore -framework QtGui -DQT_GUI_LIB -DQT_CORE_LIB -DUNICODE -L/opt/opencv/2.4.5-qt/lib -L/usr/local/lib -c larvaDistanceMap.cpp -o larvaDistanceMap.o
	@g++-4.9 $(GCCFLAGS) -g $(INCLUDES_QT) -F/usr/local/lib -framework QtCore -framework QtGui  -DQT_GUI_LIB -DQT_CORE_LIB -DUNICODE -L/opt/opencv/2.4.5-qt/lib -L/usr/local/lib -c larvaObject.cpp -o larvaObject.o
	@g++-4.9 $(GCCFLAGS) -g $(INCLUDES_QT) -F/usr/local/lib -framework QtCore -framework QtGui  -DQT_GUI_LIB -DQT_CORE_LIB -DUNICODE -L/opt/opencv/2.4.5-qt/lib -L/usr/local/lib main.cpp -lm -lopencv_contrib -lopencv_core -lopencv_highgui -lopencv_video -lopencv_imgproc -lcvblob -o lrvTrack *.o

opt:
	@g++ $(GCCFLAGS) -fast -Wall -ansi $(INCLUDES) -L/opt/opencv/2.4.6.1-gpu/lib -L/usr/local/lib -c blobUtils.cpp -o blobUtils.o
	@g++ $(GCCFLAGS) -fast -Wall -ansi $(INCLUDES) -L/opt/opencv/2.4.6.1-gpu/lib -L/usr/local/lib -c larvaSkel.cpp -o larvaSkel.o
	@g++ $(GCCFLAGS) -fast -Wall -ansi $(INCLUDES) -L/opt/opencv/2.4.6.1-gpu/lib -L/usr/local/lib -c larvaDistanceMap.cpp -o larvaDistanceMap.o
	@g++ $(GCCFLAGS) -fast -Wall -ansi $(INCLUDES) -L/opt/opencv/2.4.6.1-gpu/lib -L/usr/local/lib -c larvaObject.cpp -o larvaObject.o
	@g++ $(GCCFLAGS) -fast -Wall -ansi $(INCLUDES) -L/opt/opencv/2.4.6.1-gpu/lib -L/usr/local/lib main.cpp $(LIBS) -o lrvTrack *.o
