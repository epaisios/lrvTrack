INCLUDES=-I/opt/opencv/2.4.5/include/opencv -I/opt/opencv/2.4.5/include -I/usr/local/include -I.
INCLUDES_QT=-I/opt/opencv/2.4.5-qt/include/opencv -I/opt/opencv/2.4.5-qt/include -I/usr/local/include -I.
GCCFLAGS=-std=c++0x -fopenmp
OPTFLAGS=-Ofast
all:
	@g++-4.9 $(GCCFLAGS) -g $(INCLUDES) -L/opt/opencv/2.4.5/lib -L/usr/local/lib -c blobUtils.cpp -o blobUtils.o
	@g++-4.9 $(GCCFLAGS) -g $(INCLUDES) -L/opt/opencv/2.4.5/lib -L/usr/local/lib -c larvaSkel.cpp -o larvaSkel.o
	@g++-4.9 $(GCCFLAGS) -g $(INCLUDES) -L/opt/opencv/2.4.5/lib -L/usr/local/lib -c larvaDistanceMap.cpp -o larvaDistanceMap.o
	@g++-4.9 $(GCCFLAGS) -g $(INCLUDES) -L/opt/opencv/2.4.5/lib -L/usr/local/lib -c larvaObject.cpp -o larvaObject.o
	@g++-4.9 $(GCCFLAGS) -g $(INCLUDES) -L/opt/opencv/2.4.5/lib -L/usr/local/lib main.cpp -lm -lopencv_contrib -lopencv_core -lopencv_highgui -lopencv_video -lopencv_imgproc -lcvblob -o lrvTrack *.o

clean: 
	@rm -f *.o lrvTrack

qt:
	@g++-4.9 $(GCCFLAGS) -g $(INCLUDES_QT) -F/usr/local/lib -framework QtCore -framework QtGui -DQT_GUI_LIB -DQT_CORE_LIB -DUNICODE -L/opt/opencv/2.4.5-qt/lib -L/usr/local/lib -c blobUtils.cpp -o blobUtils.o
	@g++-4.9 $(GCCFLAGS) -g $(INCLUDES_QT) -F/usr/local/lib -framework QtCore -framework QtGui -DQT_GUI_LIB -DQT_CORE_LIB -DUNICODE -L/opt/opencv/2.4.5-qt/lib -L/usr/local/lib -c larvaSkel.cpp -o larvaSkel.o
	@g++-4.9 $(GCCFLAGS) -g $(INCLUDES_QT) -F/usr/local/lib -framework QtCore -framework QtGui -DQT_GUI_LIB -DQT_CORE_LIB -DUNICODE -L/opt/opencv/2.4.5-qt/lib -L/usr/local/lib -c larvaDistanceMap.cpp -o larvaDistanceMap.o
	@g++-4.9 $(GCCFLAGS) -g $(INCLUDES_QT) -F/usr/local/lib -framework QtCore -framework QtGui  -DQT_GUI_LIB -DQT_CORE_LIB -DUNICODE -L/opt/opencv/2.4.5-qt/lib -L/usr/local/lib -c larvaObject.cpp -o larvaObject.o
	@g++-4.9 $(GCCFLAGS) -g $(INCLUDES_QT) -F/usr/local/lib -framework QtCore -framework QtGui  -DQT_GUI_LIB -DQT_CORE_LIB -DUNICODE -L/opt/opencv/2.4.5-qt/lib -L/usr/local/lib main.cpp -lm -lopencv_contrib -lopencv_core -lopencv_highgui -lopencv_video -lopencv_imgproc -lcvblob -o lrvTrack *.o

opt: 
	@g++-4.9 $(GCCFLAGS) $(OPTFLAGS) -g $(INCLUDES) -L/opt/opencv/2.4.5/lib -L/usr/local/lib -c blobUtils.cpp -o blobUtils.o
	@g++-4.9 $(GCCFLAGS) $(OPTFLAGS) -g $(INCLUDES) -L/opt/opencv/2.4.5/lib -L/usr/local/lib -c larvaSkel.cpp -o larvaSkel.o
	@g++-4.9 $(GCCFLAGS) $(OPTFLAGS) -g $(INCLUDES) -L/opt/opencv/2.4.5/lib -L/usr/local/lib -c larvaDistanceMap.cpp -o larvaDistanceMap.o
	@g++-4.9 $(GCCFLAGS) $(OPTFLAGS) -g $(INCLUDES) -L/opt/opencv/2.4.5/lib -L/usr/local/lib -c larvaObject.cpp -o larvaObject.o
	@g++-4.9 $(GCCFLAGS) $(OPTFLAGS) -g $(INCLUDES) -L/opt/opencv/2.4.5/lib -L/usr/local/lib main.cpp -lm -lopencv_contrib -lopencv_core -lopencv_highgui -lopencv_video -lopencv_imgproc -lcvblob -o lrvTrack *.o

