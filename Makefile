all:
	@g++-4.8 -std=c++0x -fopenmp -g -I/opt/opencv/2.4.5/include/opencv -I/opt/opencv/2.4.5/include -I. -L/opt/opencv/2.4.5/lib main.cpp -lm -lopencv_contrib -lopencv_core -lopencv_highgui -lopencv_video -lopencv_imgproc -lopencv_objdetect -lopencv_contrib -lopencv_legacy -lopencv_features2d -lcvblob -o exec2

opt: 
	@g++-4.8 -std=c++0x -fopenmp -Ofast -I/opt/opencv/2.4.5/include/opencv -I/opt/opencv/2.4.5/include -I. -L/opt/opencv/2.4.5/lib main.cpp -lm -lopencv_contrib -lopencv_core -lopencv_highgui -lopencv_video -lopencv_imgproc -lopencv_objdetect -lopencv_contrib -lopencv_legacy -lopencv_features2d -lcvblob -o exec2
