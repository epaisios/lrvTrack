all:
	@g++-4.8 -std=c++0x -fopenmp -g -I/opt/opencv/2.4.4/include/opencv -I/opt/opencv/2.4.4/include -L/opt/opencv/2.4.4/lib  main.cpp -lgsl -lgslcblas -lm -lopencv_contrib -lopencv_core -lopencv_highgui -lopencv_video -lopencv_imgproc -lopencv_objdetect -lopencv_contrib -lopencv_legacy -lcvblob -o exec2

opt: 
	@g++-4.8 -std=c++0x -fopenmp -Os -I/opt/opencv/2.4.4/include/opencv -I/opt/opencv/2.4.4/include -L/opt/opencv/2.4.4/lib  main.cpp -lgsl -lgslcblas -lm -lopencv_contrib -lopencv_core -lopencv_highgui -lopencv_video -lopencv_imgproc -lopencv_objdetect -lopencv_contrib -lopencv_legacy -lcvblob -o exec2
