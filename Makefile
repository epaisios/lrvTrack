all:
	@g++-4.9 -std=c++0x -fopenmp -g -I/opt/opencv/2.4.5/include/opencv -I/opt/opencv/2.4.5/include -I/usr/local/include -I. -L/opt/opencv/2.4.5/lib -L/usr/local/lib -c blobUtils.cpp -o blobUtils.o
	@g++-4.9 -std=c++0x -fopenmp -g -I/opt/opencv/2.4.5/include/opencv -I/opt/opencv/2.4.5/include -I/usr/local/include -I. -L/opt/opencv/2.4.5/lib -L/usr/local/lib -c larvaSkel.cpp -o larvaSkel.o
	@g++-4.9 -std=c++0x -fopenmp -g -I/opt/opencv/2.4.5/include/opencv -I/opt/opencv/2.4.5/include -I/usr/local/include -I. -L/opt/opencv/2.4.5/lib -L/usr/local/lib -c larvaDistanceMap.cpp -o larvaDistanceMap.o
	@g++-4.9 -std=c++0x -fopenmp -g -I/opt/opencv/2.4.5/include/opencv -I/opt/opencv/2.4.5/include -I/usr/local/include -I. -L/opt/opencv/2.4.5/lib -L/usr/local/lib -c larvaObject.cpp -o larvaObject.o
	@g++-4.9 -std=c++0x -fopenmp -g -I/opt/opencv/2.4.5/include/opencv -I/opt/opencv/2.4.5/include -I/usr/local/include -I. -L/opt/opencv/2.4.5/lib -L/usr/local/lib main.cpp -lm -lopencv_contrib -lopencv_core -lopencv_highgui -lopencv_video -lopencv_imgproc -lcvblob -o lrvTrack *.o

clean: 
	@rm -f *.o lrvTrack

qt:
	@g++-4.9 -std=c++0x -fopenmp -g -I/opt/opencv/2.4.5-qt/include/opencv -I/opt/opencv/2.4.5-qt/include -I/usr/local/include -I. -L/opt/opencv/2.4.5-qt/lib -L/usr/local/lib -c blobUtils.cpp -o blobUtils.o
	@g++-4.9 -std=c++0x -fopenmp -g -I/opt/opencv/2.4.5-qt/include/opencv -I/opt/opencv/2.4.5-qt/include -I/usr/local/include -I. -L/opt/opencv/2.4.5-qt/lib -L/usr/local/lib -c larvaSkel.cpp -o larvaSkel.o
	@g++-4.9 -std=c++0x -fopenmp -g -I/opt/opencv/2.4.5-qt/include/opencv -I/opt/opencv/2.4.5-qt/include -I/usr/local/include -I. -L/opt/opencv/2.4.5-qt/lib -L/usr/local/lib -c larvaDistanceMap.cpp -o larvaDistanceMap.o
	@g++-4.9 -std=c++0x -fopenmp -g -I/opt/opencv/2.4.5-qt/include/opencv -I/opt/opencv/2.4.5-qt/include -I/usr/local/include -I. -L/opt/opencv/2.4.5-qt/lib -L/usr/local/lib -c larvaObject.cpp -o larvaObject.o
	@g++-4.9 -std=c++0x -fopenmp -g -I/opt/opencv/2.4.5-qt/include/opencv -I/opt/opencv/2.4.5-qt/include -I/usr/local/include -I. -L/opt/opencv/2.4.5-qt/lib -L/usr/local/lib main.cpp -lm -lopencv_contrib -lopencv_core -lopencv_highgui -lopencv_video -lopencv_imgproc -lcvblob -o lrvTrack *.o

opt: 
	@g++-4.9 -std=c++0x -fopenmp -Ofast -I/opt/opencv/2.4.5/include/opencv -I/opt/opencv/2.4.5/include -I. -L/opt/opencv/2.4.5/lib main.cpp -lm -lopencv_contrib -lopencv_core -lopencv_highgui -lopencv_video -lopencv_imgproc -lopencv_objdetect -lopencv_contrib -lopencv_legacy -lopencv_features2d -lcvblob -o exec2

win: 
	@/usr/local/gcc-4.8.0-qt-4.8.4-for-mingw32/win32-gcc/bin/i586-mingw32-c++ -std=c++0x -Ofast -mwindows -DNO_SSE=1 -I./cvBlob/ -I/usr/local/gcc-4.8.0-qt-4.8.4-for-mingw32/win32-gcc/include -I/opt/opencv/WIN/2.4.5/include/opencv -I/opt/opencv/WIN/2.4.5/modules/core/include -I/opt/opencv/WIN/2.4.5/modules/calib3d/include -I/opt/opencv/WIN/2.4.5/modules/flann/include -I/opt/opencv/WIN/2.4.5/modules/video/include -I/opt/opencv/WIN/2.4.5/modules/legacy/include -I/opt/opencv/WIN/2.4.5/modules/objdetect/include -I/opt/opencv/WIN/2.4.5/modules/features2d/include -I/opt/opencv/WIN/2.4.5/modules/highgui/include -I/opt/opencv/WIN/2.4.5/modules/imgproc/include -I/opt/opencv/WIN/2.4.5/modules/nonfree/include -I/opt/opencv/WIN/2.4.5/include -I. -L/usr/local/gcc-4.8.0-qt-4.8.4-for-mingw32/win32-gcc/lib -L/opt/opencv/WIN/2.4.5/build/x86/mingw/lib main.cpp cvBlob/cvblob.cpp cvBlob/cvcolor.cpp cvBlob/cvcontour.cpp cvBlob/cvlabel.cpp cvBlob/cvaux.cpp -lm -lm -lopencv_contrib245 -lopencv_core245 -lopencv_highgui245 -lopencv_video245 -lopencv_imgproc245 -lopencv_objdetect245 -lopencv_contrib245 -lopencv_legacy245 -lopencv_features2d245 -o exec2.exe


