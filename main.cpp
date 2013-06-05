#include <lrvTrack.hpp>
#include <fstream>
#include <sys/time.h>
#include <iomanip>
#include "cvblob.h"
#include "lrvTrackBase.hpp"
#include "blobUtils.hpp"
#include "larvaDistanceMap.hpp"

void findHeadTail(std::vector<cv::Point> &startPoints, 
		larvaObject &lrv,
		cv::Point &Head, 
		cv::Point &Tail,
		bool force_SurroundingValSearch=false)
{
	int max=0,min=65535;
	int i=0;
	if (startPoints.size()<2)
	{
		Head.x=0.0;
		Head.y=0.0;
		Tail.x=0.0;
		Tail.y=0.0;
	}

	if( lrv.start_frame==CURRENT_FRAME || 
			lrv.inCluster.back()>0 || 
			force_SurroundingValSearch ||
			lrv.roundness.back()<=3.0 ||
			(lrv.heads.back().x == 0 && lrv.heads.back().y==0 && 
			 lrv.tails.back().x == 0 && lrv.tails.back().y==0)
		)
	{
		for (i=0; i<startPoints.size(); i++)
		{
			cv::Point p=startPoints[i];
			cvb::CvBlob blob=lrv.blobs.back();
			double tmpsz=getSurroundingSize(p,blob,grey_frame);
			if (tmpsz>max)
			{
				Tail=startPoints[i];
				max=tmpsz;
			}
			if (tmpsz<min)
			{
				Head=startPoints[i];
				min=tmpsz;
			}
		}
		if (max==min)
			std::cerr << "PROBLEM AREAS AROUND HEAD AND TAIL ARE THE SAME" << std::endl;
		if (Head==Tail && i!=2)
		{
			std::cerr << "PROBLEM HEAD AND TAIL ARE THE SAME" << std::endl;
		}
	}
	else
	{
		double hmin=65535;
		double tmin=65535;
		for (i=0; i<startPoints.size(); i++)
		{
			double hdiff=diff(lrv.heads.back(),startPoints[i]);
			double tdiff=diff(lrv.tails.back(),startPoints[i]);
			if (hdiff<hmin)
			{
				Head=startPoints[i];
				hmin=hdiff;
			}
			if (tdiff<tmin)
			{
				Tail=startPoints[i];
				tmin=tdiff;
			}

		}
	}
}

void updateLarvae(cvb::CvBlobs &In, cvb::CvBlobs &Prev)
{
	cvb::CvBlobs::iterator it=In.begin();
	std::vector<unsigned int> larvaeInClusters;

	while (it!=In.end())
	{
		unsigned int ID=(*it).first;
		cvb::CvBlob blob=*((*it).second);
		cv::Mat larvaROI,cntPoints;
		createLarvaContour(larvaROI,blob);
		createLarvaContourPoints(cntPoints,blob);

		//cv::imshow("Current Larva",larvaROI);
		//cv::imshow("Contour Points",cntPoints);
		//cv::waitKey(1);
		std::map<unsigned int,larvaObject>::iterator curLarva;
		// NEW LARVA OBJECT!
		if ((curLarva=detected_larvae.find(ID))==detected_larvae.end())
		{
			// Create and allocate the new object
			larvaObject newLarva;
			//
			// Initialize the speed to 0
			newLarva.centroid_speed_x.push_back(0);
			newLarva.centroid_speed_y.push_back(0);

			// State that the larva is not in a blob
			newLarva.inCluster.push_back(false);

			cv::Point centroid=cv::Point(
					static_cast<int>(blob.centroid.x-blob.minx+ROI_PADDING+0.5),
					static_cast<int>(blob.centroid.y-blob.miny+ROI_PADDING+0.5)
					);
			newLarva.centroids.push_back(centroid);
			newLarva.lifetimeWithStats++;
			newLarva.lastIdxWithStats=CURRENT_FRAME;
			// Set the frame of it's existence
			newLarva.start_frame=CURRENT_FRAME;

			// Add the blob of the larva to its blob history
			newLarva.blobs.push_back(blob);
			// Give the larva the necessary ID
			newLarva.larva_ID=ID;

			if (detected_clusters.find(ID)==detected_clusters.end())
			{
				newLarva.isCluster=false;

				//Initialize the area values
				newLarva.area.push_back(blob.area);
				newLarva.area_mean=blob.area;
				newLarva.area_sum=blob.area;
				newLarva.area_max=newLarva.area_min=blob.area;
				newLarva.area_min=newLarva.area_min=blob.area;

				double greyVal=getGreyValue(larvaROI,blob,grey_frame);
				newLarva.grey_value.push_back(greyVal);
				newLarva.grey_value_mean = greyVal; 
				newLarva.grey_value_sum= greyVal; 
				newLarva.grey_value_max = greyVal; 
				newLarva.grey_value_min = greyVal; 

				double perimeter=getPerimeter(blob);
				newLarva.perimeter.push_back(perimeter);
				newLarva.perimeter_mean=perimeter;
				newLarva.perimeter_sum=perimeter;
				newLarva.perimeter_max=perimeter;
				newLarva.perimeter_min=perimeter;

				// Initialize the speed to 0
				newLarva.centroid_speed_x.push_back(0);
				newLarva.centroid_speed_y.push_back(0);

				newLarva.roundness.push_back((perimeter*perimeter)/(2*CV_PI*blob.area));

				// State that the larva is not in a blob
				newLarva.inCluster.push_back(false);

				cv::Point centroid=cv::Point(
						static_cast<int>(blob.centroid.x-blob.minx+ROI_PADDING+0.5),
						static_cast<int>(blob.centroid.y-blob.miny+ROI_PADDING+0.5)
						);
				newLarva.centroids.push_back(centroid);
				larvaSkel newLarvaSkel(larvaROI,centroid);
				newLarva.lrvskels.push_back(newLarvaSkel);

				// In this block we compute the inner distances for the larva
				std::vector<cv::Point> cntPoints;
				blobToPointVector(blob,cntPoints);
				larvaDistanceMap Distances(cntPoints);
				computeInnerDistances(blob,Distances,newLarvaSkel.MidPoint);
				newLarva.length.push_back(Distances.MaxDist);
				newLarva.length_mean = Distances.MaxDist;
				newLarva.length_sum= Distances.MaxDist;
				newLarva.length_max = Distances.MaxDist;
				newLarva.length_min = Distances.MaxDist;
				PointPair MAXPair=Distances.MaxDistPoints;

				cv::Point Head,Tail;
				std::vector<cv::Point> startPoints;
				startPoints.push_back(cv::Point(
							MAXPair.first.x-blob.minx+ROI_PADDING,
							MAXPair.first.y-blob.miny+ROI_PADDING)
						);
				startPoints.push_back(cv::Point(
							MAXPair.second.x-blob.minx+ROI_PADDING,
							MAXPair.second.y-blob.miny+ROI_PADDING)
						);
        newLarva.angular_speed.push_back(0);

				findHeadTail(startPoints,newLarva,Head,Tail);
				newLarva.heads.push_back(Head);
				newLarva.tails.push_back(Tail);

				newLarva.width.push_back(Distances.WidthDist);
				newLarva.width_mean = Distances.WidthDist;
				newLarva.width_sum= Distances.WidthDist;
				newLarva.width_max = Distances.WidthDist;
				newLarva.width_min = Distances.WidthDist;


				if(DEBUG_INFO!=0)
				{
					std::cout << CURRENT_FRAME <<
						" , " << newLarva.larva_ID << 
						" , " << blob.area << 
						" , " << Distances.MaxDist << 
						" , " << greyVal << 
						" , " << perimeter << 
						" , " << Distances.WidthDist << 
						" , " << newLarva.roundness.back() << 
						" , " << newLarva.inCluster.back() <<
						std::endl;
				}
			}
			else
			{
				// Larva is a blob
				// IMPORTANT: When the larva IS a blob NOTHING is updated
				//            The only fields that get updated are:
				//             - isCluster
				//             - centroid
				//             - centroid speeds
				//
				// IMPORTANT: When the larva is INSIDE a blob NOTHING is updated
				//            The only field that gets updated is:
				//             - inCluster
				//             - blob
				//
				newLarva.isCluster=true;
				std::vector<unsigned int>::iterator cl;
				for ( cl=detected_clusters[ID].begin() ;cl!=detected_clusters[ID].end() ; cl++)
				{
					larvaObject &clusterLarva=detected_larvae[*cl];
					clusterLarva.inCluster.push_back(ID);
					clusterLarva.blobs.push_back(blob);
					if(DEBUG_INFO!=0)
					{
						std::cout << CURRENT_FRAME <<
							" , " << clusterLarva.larva_ID << 
							" , " << blob.area << 
							" , " << clusterLarva.length.back() << 
							" , " << clusterLarva.grey_value.back() << 
							" , " << clusterLarva.perimeter.back() << 
							" , " << clusterLarva.width.back() << 
							" , " << clusterLarva.roundness.back() << 
							" , " << clusterLarva.inCluster.back() << 
							std::endl;
					}
				}

			}
			detected_larvae[ID]=newLarva;
		}
		// UPDATED LARVA OBJECT!
		else
		{
			//Reference to current larva
			larvaObject &cur_larva=(*curLarva).second;
			//Pointer for the previous blob
			cvb::CvBlob preBlob = cur_larva.blobs.back();

			// Set the ID of the larvaObject to the ID found TODO:Probably unnecessary
			cur_larva.larva_ID=ID;

			// Add the current blob to the blobs history of the larva
			cur_larva.blobs.push_back(blob);

			// Create the skeleton of the larva and add it to the skeletons history
			cv::Point centroid=cv::Point(
					static_cast<int>(blob.centroid.x-blob.minx+ROI_PADDING+0.5),
					static_cast<int>(blob.centroid.y-blob.miny+ROI_PADDING+0.5)
					);

			cur_larva.centroids.push_back(centroid);

			// Update centroid_speeds (in pixel per second per axis)
			cur_larva.centroid_speed_x.push_back(
					(blob.centroid.x - preBlob.centroid.x)/FrameEllapsedTime);
			cur_larva.centroid_speed_y.push_back(
					(blob.centroid.y - preBlob.centroid.y)/FrameEllapsedTime
					);
      
      double curAngle=cvb::cvAngle(&blob);
      double preAngle=cvb::cvAngle(&preBlob);

      cur_larva.angular_speed.push_back(fabs(curAngle-preAngle)/FrameEllapsedTime);

			// Look if larva is a blob
			if ( (!cur_larva.isCluster) &&
					(detected_clusters.find(ID)==detected_clusters.end()))
			{

				cur_larva.lifetimeWithStats++;
				cur_larva.lastIdxWithStats=cur_larva.area.size();

				larvaSkel newLarvaSkel(larvaROI,centroid);
				cur_larva.lrvskels.push_back(newLarvaSkel);

				// If not then:
				//  Update area values for larva.
				cur_larva.area.push_back(blob.area);

				cur_larva.area_mean=(cur_larva.area_mean+blob.area)/2;
				cur_larva.area_sum = cur_larva.area_sum + blob.area;
				if (cur_larva.area_max < blob.area)
				{
					cur_larva.area_max=blob.area;
				}
				if (cur_larva.area_min > blob.area)
				{
					cur_larva.area_min=blob.area;
				}


				// Try to find Head and Tail

				// Point coordinates for head/tail
				cv::Point Head,Tail;

				// Map to keep the distances of each point to all the others
				// Pair of points to keep the points with the Maximum distance (i.e. head and tail :) )
				std::vector<cv::Point> cntPoints;
				blobToPointVector(blob,cntPoints);
				larvaDistanceMap Distances(cntPoints);
				// Compute all the inner distances for the larva
				computeInnerDistances(blob,Distances,newLarvaSkel.MidPoint);
				cur_larva.length.push_back(Distances.MaxDist);
				cur_larva.length_mean=(cur_larva.length_mean+Distances.MaxDist)/2;
				cur_larva.length_sum=cur_larva.length_sum+Distances.MaxDist;
				if (cur_larva.area_max < Distances.MaxDist)
				{
					cur_larva.area_max=Distances.MaxDist;
				}
				if (cur_larva.area_min > Distances.MaxDist)
				{
					cur_larva.area_min=Distances.MaxDist;
				}
				PointPair MAXPair=Distances.MaxDistPoints;

				cur_larva.width.push_back(Distances.WidthDist);
				cur_larva.width_mean=(cur_larva.width_mean+Distances.WidthDist)/2;
				cur_larva.width_sum=cur_larva.width_sum+Distances.WidthDist;
				if (cur_larva.area_max < Distances.WidthDist)
				{
					cur_larva.area_max=Distances.WidthDist;
				}
				if (cur_larva.area_min > Distances.WidthDist)
				{
					cur_larva.area_min=Distances.WidthDist;
				}

				double greyVal=getGreyValue(larvaROI,blob,grey_frame);
				cur_larva.grey_value.push_back(greyVal);
				cur_larva.grey_value_mean=(cur_larva.grey_value_mean+greyVal)/2; 
				cur_larva.grey_value_sum=cur_larva.grey_value_sum+greyVal; 
				if (cur_larva.area_max < greyVal)
				{
					cur_larva.grey_value_max=greyVal;
				}
				if (cur_larva.area_min > greyVal)
				{
					cur_larva.grey_value_min=greyVal;
				}

				double perimeter=getPerimeter(blob);
				cur_larva.perimeter.push_back(perimeter);
				cur_larva.perimeter_mean=(cur_larva.perimeter_mean+perimeter)/2; 
				cur_larva.perimeter_sum=cur_larva.perimeter_sum+perimeter; 
				if (cur_larva.area_max < perimeter)
				{
					cur_larva.perimeter_max=perimeter;
				}
				if (cur_larva.area_min > perimeter)
				{
					cur_larva.perimeter_min=perimeter;
				}

				cur_larva.roundness.push_back((perimeter*perimeter)/(2*CV_PI*blob.area));

				// Construct a vector of points including both 
				// head and tail to decide which is which
				std::vector<cv::Point> startPoints;

				// Points must be corrected to match the ROI including the padding set
				startPoints.push_back(cv::Point(
							MAXPair.first.x-blob.minx+ROI_PADDING,
							MAXPair.first.y-blob.miny+ROI_PADDING)
						);
				startPoints.push_back(cv::Point(
							MAXPair.second.x-blob.minx+ROI_PADDING,
							MAXPair.second.y-blob.miny+ROI_PADDING)
						);
				//execute the function to find which is which and assign appropriately
				//to Head/Tail.
				//findHeadTail(startPoints,cur_larva,Head,Tail,true);
				findHeadTail(startPoints,cur_larva,Head,Tail);

				//Add head and tail to the history
				cur_larva.heads.push_back(Head);
				cur_larva.tails.push_back(Tail);

				//state that larva is not detected as part of a blob of larvae
				//NOTE: It is important to perform this value setting !AFTER!
				//      searching for the head tail because, the head tail 
				//      search behaves differently if the larva was previously
				//      part of a blob.
				cur_larva.inCluster.push_back(0);

				if(DEBUG_INFO!=0)
				{
					std::cout << CURRENT_FRAME <<
						" , " << cur_larva.larva_ID << 
						" , " << blob.area << 
						" , " << cur_larva.length.back() << 
						" , " << cur_larva.grey_value.back() << 
						" , " << cur_larva.perimeter.back() << 
						" , " << cur_larva.width.back() << 
						" , " << cur_larva.roundness.back() << 
						" , " << cur_larva.inCluster.back() << 
						std::endl;
				}

			}
			else
			{
				// Larva is a blob
				// IMPORTANT: When the larva IS a blob NOTHING is updated
				//            The only fields that get updated are:
				//             - isCluster
				//             - blob
				//             - centroid
				//             - centroid speeds
				//
				// IMPORTANT: When the larva is INSIDE a blob NOTHING is updated
				//            The only field that gets updated is:
				//             - inCluster
				//             - blob
				//
				cur_larva.isCluster=true;
				std::vector<unsigned int>::iterator cl;
				for ( cl=detected_clusters[ID].begin() ;cl!=detected_clusters[ID].end() ; cl++)
				{
					larvaObject &clusterLarva=detected_larvae[*cl];
					clusterLarva.inCluster.push_back(ID);
					clusterLarva.blobs.push_back(blob);
					if(DEBUG_INFO!=0)
					{
						std::cout << CURRENT_FRAME <<
							" , " << clusterLarva.larva_ID << 
							" , " << blob.area << 
							" , " << clusterLarva.length.back() << 
							" , " << clusterLarva.grey_value.back() << 
							" , " << clusterLarva.perimeter.back() << 
							" , " << clusterLarva.width.back() << 
							" , " << clusterLarva.roundness.back() << 
							" , " << clusterLarva.inCluster.back() << 
							std::endl;
					}
				}
			}
		}
		it++;
	}
}

inline double SQUARE(double n){
	return n*n;
}

//Currently only for the 4 sized vector
inline int vmin(std::vector<double>& vals, std::vector<int>& seq)
{
	for (int i=0 ; i<vals.size() ; i++ )
	{
		int sum=0;
		for (int j=0; j<vals.size() ; j++)
		{
			if (vals[i]>vals[j])
			{
				++sum;
			}
		}
		seq[sum]=i;
	}
}

namespace std{
	template <typename number >
		std::string printVector(std::vector<number> vec)
		{
			std::stringstream sstm;
			bool const is_number= std::is_arithmetic<number>::value;
			static_assert( is_number, "Provided type is not an arithmetic type");
			sstm << "[" ;
			typename std::vector<number>::const_iterator i=vec.begin();
			sstm << *i ;
			i++;
			for( ; i != vec.end(); ++i)
				sstm << ","<< *i;
			sstm << "]"; 
			return sstm.str();
		}
}

void diverge_match(unsigned int &candidate_larva_a, 
		unsigned int &candidate_larva_b,
		cvb::CvBlob  *newLarva1, 
		cvb::CvBlob *newLarva2)
{
	larvaObject &LarvaA=detected_larvae[candidate_larva_a];
	larvaObject &LarvaB=detected_larvae[candidate_larva_b];

	cv::Point mdpLarva1, mdpLarva2;
	mdpLarva1.x=newLarva1->centroid.x;
	mdpLarva1.y=newLarva1->centroid.y;
	mdpLarva2.x=newLarva2->centroid.x;
	mdpLarva2.y=newLarva2->centroid.y;

	double size_a=LarvaA.area_sum/LarvaA.area.size();
	double size_b=LarvaB.area_sum/LarvaB.area.size();;
	double size_1=newLarva1->area;
	double size_2=newLarva2->area;

	cv::Mat larvaROI1,larvaROI2;
	createLarvaContour(larvaROI1,*newLarva1);
	createLarvaContour(larvaROI2,*newLarva2);
	double grey_value_a=LarvaA.grey_value_sum/LarvaA.grey_value.size();
	double grey_value_b=LarvaB.grey_value_sum/LarvaB.grey_value.size();
	double grey_value_1=getGreyValue(larvaROI1,*newLarva1,grey_frame);
	double grey_value_2=getGreyValue(larvaROI2,*newLarva2,grey_frame);

	// Length
	std::vector<cv::Point> newLarva1Points;
	std::vector<cv::Point> newLarva2Points;
	blobToPointVector(*newLarva1,newLarva1Points);
	blobToPointVector(*newLarva2,newLarva2Points);
	larvaDistanceMap dstLarva1(newLarva1Points), dstLarva2(newLarva2Points);

	cv::Point centroid1;
	centroid1.x=newLarva1->centroid.x;
	centroid1.y=newLarva1->centroid.y;
	cv::Point centroid2;
	centroid2.x=newLarva2->centroid.x;
	centroid2.y=newLarva2->centroid.y;

	larvaSkel newLarvaSkel1(larvaROI1,centroid1);
	larvaSkel newLarvaSkel2(larvaROI2,centroid2);

	computeInnerDistances(*newLarva1,dstLarva1,newLarvaSkel1.MidPoint);
	computeInnerDistances(*newLarva2,dstLarva2,newLarvaSkel2.MidPoint);


	double length_a=LarvaA.length_sum/LarvaA.length.size();
	double length_b=LarvaB.length_sum/LarvaB.length.size();
	double length_1=dstLarva1.MaxDist;
	double length_2=dstLarva2.MaxDist;

	double perimeter_a=LarvaA.perimeter_sum/LarvaA.perimeter.size();
	double perimeter_b=LarvaB.perimeter_sum/LarvaB.perimeter.size();
	double perimeter_1=getPerimeter(*newLarva1);
	double perimeter_2=getPerimeter(*newLarva2);

	double width_a=LarvaA.width_sum/LarvaA.width.size();
	double width_b=LarvaB.width_sum/LarvaB.width.size();
	double width_1=dstLarva1.WidthDist;
	double width_2=dstLarva2.WidthDist;

	cv::Mat InputArrayA;
	cv::Mat InputArrayB;
	cv::hconcat(cv::Mat(LarvaA.area),cv::Mat(LarvaA.grey_value),InputArrayA);
	cv::hconcat(InputArrayA,cv::Mat(LarvaA.length),InputArrayA);
	cv::hconcat(InputArrayA,cv::Mat(LarvaA.perimeter),InputArrayA);
	cv::hconcat(InputArrayA,cv::Mat(LarvaA.width),InputArrayA);

	cv::hconcat(cv::Mat(LarvaB.area),cv::Mat(LarvaB.grey_value),InputArrayB);
	cv::hconcat(InputArrayB,cv::Mat(LarvaB.length),InputArrayB);
	cv::hconcat(InputArrayB,cv::Mat(LarvaB.perimeter),InputArrayB);
	cv::hconcat(InputArrayB,cv::Mat(LarvaB.width),InputArrayB);

	/*
		 std::cerr << InputArrayA << std::endl;
		 std::cerr << InputArrayB << std::endl;
		 std::cerr << "===========================================" << std::endl;
		 */
	std::vector<double> meanVecA={size_a , grey_value_a , length_a, perimeter_a, width_a };
	std::vector<double> meanVecB={size_b , grey_value_b , length_b, perimeter_b ,width_b };
	cv::Mat meanMatA(meanVecA);
	cv::Mat meanMatB(meanVecB);
	cv::Mat meanTMatA, meanTMatB;
	cv::transpose(meanMatA,meanTMatA);
	cv::transpose(meanMatB,meanTMatB);
	cv::Mat covarMatA;
	cv::Mat covarMatB;
	/* 
		 printVector(LarvaA.area);
		 printVector(LarvaA.grey_value);
		 printVector(LarvaA.perimeter);

		 printVector(LarvaB.area);
		 printVector(LarvaB.grey_value);
		 printVector(LarvaB.perimeter);
		 */
	cv::calcCovarMatrix(InputArrayA, covarMatA,meanTMatA,CV_COVAR_ROWS|CV_COVAR_NORMAL|CV_COVAR_USE_AVG);
	cv::calcCovarMatrix(InputArrayB, covarMatB,meanTMatB,CV_COVAR_ROWS|CV_COVAR_NORMAL|CV_COVAR_USE_AVG);

	/*
		 std::cerr << "CovarMatA" << covarMatA << std::endl;
		 std::cerr << "CovarMatB" << covarMatB << std::endl;
		 */
	cv::invert(covarMatA,covarMatA,cv::DECOMP_SVD);
	cv::invert(covarMatB,covarMatB,cv::DECOMP_SVD);

	std::vector<double> vec1 = { size_1, grey_value_1 , length_1, perimeter_1, width_1 };
	std::vector<double> vec2 = { size_2, grey_value_2 , length_2, perimeter_2, width_2 };

	cv::Mat mat1(vec1);
	cv::Mat mat2(vec2);

	double DistA1 = cv::Mahalanobis(mat1, meanMatA, covarMatA);
	double DistA2 = cv::Mahalanobis(mat2, meanMatA, covarMatA);
	double DistB1 = cv::Mahalanobis(mat1, meanMatB, covarMatB);
	double DistB2 = cv::Mahalanobis(mat2, meanMatB, covarMatB);

	double MHDiff = (DistA1+DistB2) - (DistA2+DistB1);
	std::vector<double> vals = { DistA1, DistB2, DistA2, DistB1 };
	std::vector<int> minv;
	minv.reserve(4);
	vmin(vals,minv);
	int mini=minv[0];
	if (DistA1+DistB2 > DistA2+DistB1)
		/*
			 double miniDiff = fabs(vals[minv[0]] - vals[minv[1]]);
			 if (mini >=2 && MHDiff < 0 )
			 {
			 if (fabs(MHDiff) > miniDiff)
			 {
			 mini=0;
			 }
			 }
			 else if (mini < 2 && MHDiff > 0 )
			 {
			 if (fabs(MHDiff) > miniDiff)
			 {
			 mini=2;
			 }
			 }

			 if (mini >=2 )
			 */
	{
		candidate_larva_a=newLarva2->label;
		candidate_larva_b=newLarva1->label;
		std::cerr << "MH<" << CURRENT_FRAME << "> : Assigning " << LarvaA.larva_ID << " -> 2 [" << newLarva2->centroid.x << ", " << newLarva2->centroid.y << "] MHDiff: " << MHDiff << " min i: " << mini <<  std::endl;
		std::cerr << "MH<" << CURRENT_FRAME << "> : Assigning " << LarvaB.larva_ID << " -> 1 [" << newLarva1->centroid.x << ", " << newLarva1->centroid.y << "] Dists: A1: " << DistA1 << " B2: " << DistB2 << ", A2: " << DistA2 << " B1: " << DistB1 << std::endl;
	}
			 else
			 {
				 candidate_larva_a=newLarva1->label;
				 candidate_larva_b=newLarva2->label;
				 std::cerr << "MH<" << CURRENT_FRAME << "> : Assigning " << LarvaA.larva_ID << " -> 1 [" << newLarva1->centroid.x << ", " << newLarva1->centroid.y << "] MHDiff: " << MHDiff << " min i: " << mini << std::endl;
				 std::cerr << "MH<" << CURRENT_FRAME << "> : Assigning " << LarvaB.larva_ID << " -> 2 [" << newLarva2->centroid.x << ", " << newLarva2->centroid.y << "] Dists: A1: " << DistA1 << " B2: " << DistB2 << ", A2: " << DistA2 << " B1: " << DistB1 << std::endl;
			 }
}



void larvae_track(cvb::CvBlobs &In,cvb::CvBlobs &Prev,cvb::CvBlobs &out)
{

	// Map to keep track of the larvae that were tracked succesfully so far
	std::map<unsigned int,std::vector<unsigned int> > used_map;
  // Resetting current_clusters and current_diverged to empty
  current_clusters.clear();
  current_diverged.clear();
  current_new.clear();
  
	cvb::CvBlobs::iterator prevIt=Prev.begin();
	while (prevIt!=Prev.end())
	{
		// preBlob is the larva we are trying to match at this point.
		// preBlob belongs to the previous frame and we are trying to find
		// which larva of the current frame is more likely to be the one
		cvb::CvBlob *preBlob=(*prevIt).second;

		// Used to store the label of the larvae which is closest to the preBlob
		unsigned int minLabel;
		// Helper variables to store and compute distances
		double XVAL,YVAL,cur=0;
		// Stores the minimal distance found (actually the manhattan distance).
		double min=200;

		// Here we check if the *prevIt was already assigned a match.
		// This can occur when diverging where the two diverging worms are
		// both assigned at the same time.

		if(out.find((*prevIt).first)!=out.end())
		{
			prevIt++;
			continue;
		}
		// Here we start looking for the blob in "In" (i.e. the latest larvae set)
		// that is closest to the preBlob. Typically (in 15fps and up ) distances
		// of matching larvae should be quite small (0.0X for the centroids).
		cvb::CvBlobs::iterator it=In.begin();
		while (it!=In.end())
		{
			cvb::CvBlob *blob;
			blob=((*it).second);
			// MIN_DISCARD_DISTANCE is used to quickly exclude larvae that are
			// too far to be considered
			if (((XVAL=fabs(blob->centroid.x - preBlob->centroid.x)) < MIN_DISCARD_DISTANCE ) &&
					((YVAL=fabs(blob->centroid.y - preBlob->centroid.y)) < MIN_DISCARD_DISTANCE ))
			{
				// we only use the manhattan distance which in this case should be sufficient
				cur=XVAL+YVAL;
				if (cur < min)
				{
					min=cur;
					minLabel=(*it).first;
				}
			}
			it++;
		}

		// We went through all the current batch of larvae and found minimum distance min
		// for the larva with label minLabel
		if (min<=5)
		{
			// Short distance. Indicates:
			// 1) Most likely normal motion
			// 2) Joining of larvae without much change in the centroid (TODO)
			// 3) Spliting of larvae without much change in the centoid (TODO)
			out[(*prevIt).first]=In[minLabel];
			out[(*prevIt).first]->label=(*prevIt).first;
			used_map[minLabel].push_back((*prevIt).first);
			//inactive[(*prevIt).first]=0;
		}
		else if (min<150)
		{
			// Rather large jump. Indicates:
			// 1) larvae converged
			// 2) larvae diverged
			// 3) TODO: Larva dissapeared and matched wrongly with another...

			// Check for converging
			// We check if more than one larvae match the new centroid and
			//  whether the middle of their centroids matches the centroid we have.

			// This is used only in the diverging case to point to the matching cluster.
			std::map< unsigned int, std::vector<unsigned int> >::iterator cluster;

			// If the larva was assigned previously to another one it may be that
			// we have a cluster! (used_map for the minLabel larvae has been assigned to
			// at least one larvae).
			if (used_map[minLabel].size()>0)
			{
				cvb::CvBlob *blob;
				blob=In[minLabel];
				// One or more larvae have matched to the new blob with ID minLabel
				// We will now test whether the centroids of the other matches give a
				//   good aproximation of the new centroid.
				//   TODO: Take care of case of more than 2 merging at one point.
				double newx = ((*prevIt).second->centroid.x + 
						Prev[used_map[minLabel][0]]->centroid.x);
				double newy = ((*prevIt).second->centroid.y + 
						Prev[used_map[minLabel][0]]->centroid.y);
				double XDIF,YDIF;
				// TODO: Use correct centroid calculation based on area not just middle
				if (((XDIF=fabs(2*blob->centroid.x - newx)) < 10 ) &&
						((YDIF=fabs(2*blob->centroid.y - newy)) < 10))
				{
					// We're in luck! The center's match :) We have a cluster.
					used_map[minLabel].push_back((*prevIt).first);
					// Add both larvae to the detected clusters
					detected_clusters[blob->label].push_back(used_map[minLabel][0]);
					detected_clusters[blob->label].push_back(used_map[minLabel][1]);
					current_clusters[blob->label].push_back(used_map[minLabel][0]);
					current_clusters[blob->label].push_back(used_map[minLabel][1]);
					// Create a second "instance" of the blob to represent both larvae
					//out[(*prevIt).first]=In[minLabel];
				}
				else
				{
					// May be that one larva dissapeared and matched wrongly with another nearby
					// In this case we should only keep the one that is closest!
					std::cerr << "PROBLEM!!! : Check how to resolve these cases where centers after do not match" << std::endl;
				}

			}
			// Check for diverging:
			// If the detected_clusters map contains a cluster for the 
			// previous blob with the one we are examining we might be 
			// looking into larvae that have separated.
			// TODO: nested clusters would also be nice...
			else if ((cluster=detected_clusters.find((*prevIt).first))!=detected_clusters.end())
			{
				it=In.begin();
				unsigned int secondDivergent=0;
				double XVAL,YVAL,cur=0;
				double min=40;
				while (it!=In.end())
				{
					if ((*it).first!=minLabel)
					{
						cvb::CvBlob *blob;
						blob=(*it).second;
						if (((XVAL=fabs(blob->centroid.x - preBlob->centroid.x)) < 90) &&
								((YVAL=fabs(blob->centroid.y - preBlob->centroid.y)) < 90))
						{
							cur=XVAL+YVAL;
							if (cur < min)
							{
								min=cur;
								secondDivergent=(*it).first;
							}
						}
					}
					it++;
				}
				if(secondDivergent==0)
				{
					std::cerr << "PROBLEM!!! : Did not find second divergent." << std::endl;
				}
				else
				{
					double newx = (In[minLabel]->centroid.x + In[secondDivergent]->centroid.x);
					double newy = (In[minLabel]->centroid.y + In[secondDivergent]->centroid.y);
					// TODO: Use correct centroid calculation based on area not just middle
					if (((XVAL=fabs(2*(*prevIt).second->centroid.x - newx)) < 20 ) &&
							((YVAL=fabs(2*(*prevIt).second->centroid.y - newy)) < 20))
					{
						// We're in luck! The center's match :) We have a diverging cluster.
						used_map[minLabel].push_back((*prevIt).first);
						used_map[secondDivergent].push_back((*prevIt).first);
						// Matching function here:
						unsigned int Larva1=(*cluster).second[0];
						unsigned int Larva2=(*cluster).second[1];
						diverge_match(Larva1,Larva2,In[minLabel],In[secondDivergent]);
						out[(*cluster).second[0]]=In[Larva1];
						out[(*cluster).second[1]]=In[Larva2];
						out[(*cluster).second[0]]->label=Larva1;
						out[(*cluster).second[1]]->label=Larva2;
            current_diverged[cluster->first].push_back(minLabel);
            current_diverged[cluster->first].push_back(secondDivergent);
						detected_clusters.erase(cluster);
					}
					else
					{
						std::cerr << "PROBLEM!!! : Check how to resolve these cases where centers do not match" << std::endl;
					}
				}
			}
			else
			{
				// Here we have the following cases:
				// 1) First Larva Spotted belonging to a too large cluster 
				//    -- Not the case since we have no unseen clusters
				// 1) Divergence of larvae that were clustered from the start <TODO>
				// 2) A small jump. Perhaps framerate related... <TODO>
				//used_map[minLabel].push_back((*prevIt).first);
				//out[(*prevIt).first]=In[minLabel];
				//out[(*prevIt).first]->label=(*prevIt).first;

				used_map[minLabel].push_back((*prevIt).first);
				out[++LARVAE_COUNT]=In[minLabel];
				out[LARVAE_COUNT]->label=LARVAE_COUNT;
        current_new.push_back(LARVAE_COUNT);
			}
			//out[(*prevIt).first]=In[minLabel];
			//out[(*prevIt).first]->label=(*prevIt).first;
			//used_map[minLabel]++;
			//detected_clusters[minLabel].push_back((*prevIt).first);
		}
		prevIt++;
	}

	cvb::CvBlobs::iterator it=In.begin();
	while (it!=In.end())
	{
		int ID=(*it).first;
		if (used_map.find(ID)==used_map.end())
		{
			out[++LARVAE_COUNT]=(*it).second;
			out[LARVAE_COUNT]->label=LARVAE_COUNT;
		}
		it++;
	}
	updateLarvae(out,Prev);

}


int main(int argv, char* argc[])
{
	bool SHOWTAGS=true;
	bool TRACK=false;
	bool STEP=true;
	bool showimg=true;
	int START_FRAME=0;
	cv::VideoCapture capture;
	char cam[3]="-c";
	std::ofstream summary;
	summary.open("lalala.summary");
	if(!strncmp(argc[1],cam,2))
	{
		capture.open(CV_CAP_DC1394);
		capture.set(CV_CAP_PROP_FPS,24);
		capture.set(CV_CAP_PROP_FRAME_WIDTH,1280);
		capture.set(CV_CAP_PROP_FRAME_HEIGHT,1024);
	}
	else
	{
		capture.open(argc[1]);
	}

	if (capture.get(CV_CAP_PROP_FPS)==0)
	{
		VIDEO_FPS=24.1;
	}
	else
	{
		if ( capture.get(CV_CAP_PROP_FPS) > 30 )
		{
			VIDEO_FPS=capture.get(CV_CAP_PROP_FPS)/2;
			//std::cout << VIDEO_FPS << std::endl;
		}
		else
		{
			VIDEO_FPS=capture.get(CV_CAP_PROP_FPS);
			//std::cout << VIDEO_FPS << std::endl;
		}
	}
	cv::Mat grey_bgFrame;
	unsigned int fg_threshold=10;
	unsigned int thresholdlow=90;
	unsigned int thresholdhigh=255;
	LARVAE_COUNT=0;

	if (!capture.isOpened())
		return 1;

	bool stop(false);
	cv::namedWindow("Extracted Frame");

	capture.read(bgFrame);
	cvtColor(bgFrame,bgFrame,CV_BGR2GRAY);
	if(bgFrame.channels()>1)
	{
		cvtColor(bgFrame,grey_bgFrame,CV_BGR2GRAY);
	}
	else
	{
		bgFrame.copyTo(grey_bgFrame);
	}

	//std::cout << bgFrame.rows << "," << bgFrame.cols << std::endl;
	std::vector<cv::Vec3f> circles;
	cv::HoughCircles(grey_bgFrame, circles, CV_HOUGH_GRADIENT,
			1,   // accumulator resolution (size of the image / 2)
			850,  // minimum distance between two circles
			50, // Canny high threshold
			180, // minimum number of votes
			bgFrame.rows/3, bgFrame.rows/2); // min and max radiusV

	std::vector<cv::Vec3f>::
		const_iterator itc= circles.begin();
	//std::cout << circles.size() << std::endl;
	if(circles.size()>0)
	{
		while (itc!=circles.end())
		{
			cv::circle(grey_bgFrame,
					cv::Point((*itc)[0], (*itc)[1]), // circle centre
					(*itc)[2]/2,       // circle radius
					cv::Scalar(), // color
					-1);              // thickness
			++itc;
		}
	}

	cvb::CvBlobs preBlobs;
	while (!stop)
	{
		// read next frame if any
		if (!capture.read(frame))
			break;

		if(TRACK)
		{
			if(START_FRAME==0)
			{
				gettimeofday(&tS,0);
				gettimeofday(&tP,0);
				START_FRAME=CURRENT_FRAME;
			}
			cvtColor(frame,frame,CV_BGR2GRAY);
			CURRENT_FRAME=capture.get(CV_CAP_PROP_POS_FRAMES);
			cv::Mat fg_frame;
			cv::Mat image;
			cv::Mat processed_frame;
			cv::Mat working_frame;
			cv::Mat fgROI;
			cv::Mat fg_image;
			cv::Mat mask;

			// Captured frame to BW (necessary for various filters and speed)
			if(frame.channels()>1)
			{
				cvtColor(frame,grey_frame,CV_BGR2GRAY);
			}
			else
			{
				frame.copyTo(grey_frame);
				cvtColor(grey_frame,frame,CV_GRAY2BGR);
			}

			//Background subtraction is done by eliminating areas outside of the
			//detected petri-dish
			// Here we subtract from the BW captured frame the background frame.
			//cv::addWeighted(grey_frame, 1.0, grey_bgFrame, -3.0, 0.0, fg_frame);
			cv::absdiff(grey_frame,grey_bgFrame,fg_frame);
			//cv::imshow("normWin",fg_frame);
			//fg_frame=grey_frame;
			fgROI=cv::Mat::zeros(fg_frame.rows , fg_frame.cols,fg_frame.depth());

			if(circles.size()>0)
			{
				cv::circle(fgROI, cv::Point(circles[0][0],circles[0][1]),int(circles[0][2]/1.0),cv::Scalar(255),-1);
				cv::circle(frame, cv::Point(circles[0][0],circles[0][1]),int(circles[0][2]/1.0),cv::Scalar(0,255,0),1);
			}

			fg_image=fg_frame&fgROI;

			cv::Mat fg_image_norm;

			cv::normalize(fg_image,fg_image_norm,0,255,cv::NORM_MINMAX);
			cv::threshold(fg_image,
					thresholded_frame,
					thresholdlow,
					thresholdhigh,
					cv::THRESH_BINARY|cv::THRESH_OTSU);
			//cv::blur(thresholded_frame,thresholded_frame,cv::Size(2,2));
			//cv::imshow("normWin",thresholded_frame);
			cvb::CvBlobs blobs;
			IplImage ipl_thresholded = thresholded_frame;
			labelImg=cvCreateImage(cvGetSize(&ipl_thresholded), IPL_DEPTH_LABEL, 1);

			unsigned int result=cvLabel(&ipl_thresholded, labelImg, blobs);
			cvb::cvFilterByArea(blobs, 45, 900);

			cvb::CvBlobs tracked_blobs;
			cvb::CvBlobs blob1;
			cvb::CvBlobs::iterator it=blobs.begin();

			if(preBlobs.size()>0)
			{
			  gettimeofday(&tC,0);
        FrameEllapsedTime = (tC.tv_sec - tP.tv_sec) + 
                            ((tC.tv_usec - tP.tv_usec)/1000000.0);
				larvae_track(blobs,preBlobs,tracked_blobs);
        gettimeofday(&tP,0);
				preBlobs=tracked_blobs;
				// output summary
				// Collect info
				gettimeofday(&tC,0);
				double elapsed=(tC.tv_sec - tS.tv_sec) + ((tC.tv_usec - tS.tv_usec)/1000000.0);
				double lifespanSUM=0;
				double speedSUM=0;
				double angularSpeedSUM=0;
				double lengthSUM=0;
				double lengthRelSUM=0;
				double widthSUM=0;
				double widthRelSUM=0;
				double aspectRatioSUM=0;
				double relAspectRatioSUM=0;
				double avgSizeSUM=0;
				double avgWiggleSUM=0;
				int larvaeToConsider=0;
				std::map<unsigned int,larvaObject>::iterator i=detected_larvae.begin();
				for (; i!=detected_larvae.end();i++)
				{
					larvaObject &cl=(*i).second;
					// avg duration
					if (cl.isCluster==true)
					 break;

					larvaeToConsider++;
					lifespanSUM+=(CURRENT_FRAME-cl.start_frame)/VIDEO_FPS;
					// avg speed
					float xvel,yvel;
					if(cl.centroid_speed_x.size()>0)
						xvel=cl.centroid_speed_x.back();
					else
						xvel=0.0;
					if(cl.centroid_speed_y.size()>0)
						yvel=cl.centroid_speed_y.back();
					else
						yvel=0.0;
					float sqvel=xvel*xvel+yvel*yvel;
					float vel[1];
					if(sqvel==0)
						vel[0]=0;
					else
					{
						ltsqrt(vel,&sqvel);
					}
					speedSUM+=vel[0];

					//avg angular speed
					angularSpeedSUM+=cl.angular_speed.back();
					
					//average length
					float length=cl.length.back();
					if(length!=0)
					{
						ltsqrt(vel,&length);
					}
          double realLength=vel[0];
					lengthSUM+=vel[0];
					lengthRelSUM+=cl.length.back()*cl.length.size()/cl.length_sum;
					
					//average width
					widthSUM+=cl.width.back();
					widthRelSUM+=cl.width.back()*cl.width.size()/cl.width_sum;
					
					//average aspect ratio
					aspectRatioSUM+=cl.width.back()/realLength;
					
					//average relative aspect ratio
					relAspectRatioSUM+=(cl.width.back()*cl.length_mean) /
					                (cl.width_mean*cl.length.back());
					//average wiggle
          double a1,a2;
          a1=angle(cl.heads.back(),
                   cl.lrvskels.back().Point20,
                   cl.tails.back());
          a2=angle(cl.heads.back(),
                   cl.lrvskels.back().Point80,
                   cl.tails.back());

					avgWiggleSUM += MAX(a1,a2);
					//average size
					avgSizeSUM+=cl.area.back();
				}

				std::setiosflags(std::ios::fixed);

				summary << CURRENT_FRAME-START_FRAME << " " ;
				summary << std::left 
                << std::fixed 
                << std::setfill('0') 
                << std::setprecision(3) 
                << elapsed 
                << "  ";

				summary << detected_larvae.size() << " ";

				summary << detected_larvae.size() << " ";

				summary << std::left 
                << std::fixed 
                << std::setfill('0') 
                << std::setprecision(1) 
                << lifespanSUM/larvaeToConsider 
                << "  ";

				//TODO angular speed
				summary << std::left 
                << std::fixed 
                << std::setfill('0') 
                << std::setprecision(2) 
                << speedSUM/larvaeToConsider 
                << " ";

				summary << std::left 
                << std::fixed 
                << std::setfill('0') 
                << std::setprecision(3) 
                << angularSpeedSUM/larvaeToConsider 
                << "  ";
				summary << std::left 
                << std::fixed 
                << std::setfill('0') 
                << std::setprecision(1) 
                << lengthSUM/larvaeToConsider 
                << " ";

				summary << std::left 
                << std::fixed 
                << std::setfill('0') 
                << std::setprecision(3) 
                << lengthRelSUM/larvaeToConsider 
                << "  ";

				summary << std::left 
                << std::fixed 
                << std::setfill('0') 
                << std::setprecision(1) 
                << widthSUM/larvaeToConsider 
                << " ";

				summary << std::left 
                << std::fixed 
                << std::setfill('0') 
                << std::setprecision(3) 
                << widthRelSUM/larvaeToConsider 
                << "  ";

				summary << std::left 
                << std::fixed 
                << std::setfill('0') 
                << std::setprecision(3) 
                << aspectRatioSUM/larvaeToConsider 
                << " ";

				summary << std::left 
                << std::fixed 
                << std::setfill('0') 
                << std::setprecision(3) 
                << relAspectRatioSUM/larvaeToConsider 
                << "  ";
        
				//TODO wiggle
				summary << std::left 
                << std::fixed 
                << std::setfill('0') 
                << std::setprecision(3) 
                << avgWiggleSUM/larvaeToConsider 
                << " ";

				summary << std::left 
                << std::fixed 
                << std::setfill('0') 
                << std::setprecision(3) 
                << avgSizeSUM/larvaeToConsider 
                << " ";
        if(current_new.size()>0 || 
           current_clusters.size()>0 ||
           current_diverged.size()>0)
        {
          std::cout << " %% ";
        }
        if(current_new.size()>0)
        {
          for (int i=0; i<current_new.size();i++)
          {
            std::cout << "0 " << current_new[i] << " ";
          }
        }
        if(current_clusters.size()>0)
        {
          std::map<unsigned int,std::vector<unsigned int> >::iterator ci;
          for (ci=current_clusters.begin(); ci!=current_clusters.end();ci++)
          {
            std::cout << ci->second[0] << " " << ci->first << " " ;
            std::cout << ci->second[1] << " " << ci->first << " " ;
          }
        }
				summary << std::endl;
			}
			else
			{
				preBlobs.clear();
				//normalize blobs for next use
				cvb::CvBlobs::iterator it=blobs.begin();
				int i=1;
				while (it!=blobs.end())
				{
					preBlobs[i]=(*it).second;
					preBlobs[i]->label=i;
					it++;
					i++;
				}
				LARVAE_COUNT=preBlobs.size();
			}

			if (SHOWTAGS!=0)
			{
				it=tracked_blobs.begin();
				while (it!=tracked_blobs.end())
				{
					std::stringstream sstm;
					cvb::CvBlob *blob=(*it).second;
					try{
						std::vector<unsigned int> &clusterMBs=detected_clusters.at((*it).first);
						sstm << printVector(clusterMBs);
						cv::putText(frame,
								sstm.str(),
								cv::Point2d(blob->centroid.x+12,blob->centroid.y+12),
								cv::FONT_HERSHEY_PLAIN,
								0.7,
								cv::Scalar(255,255,255),
								1,
								CV_AA);
					}
					catch (const std::out_of_range& oor)
					{
						//out of range. i.e. it's not there :)
						sstm << (*it).first;
						cv::putText(frame,
								sstm.str(),
								cv::Point2d(blob->centroid.x+12,blob->centroid.y+12),
								cv::FONT_HERSHEY_PLAIN,
								0.8,
								cv::Scalar(255,255,255),
								1,
								CV_AA);
					}

					// Here we define an extra PADDING
					int PAD=2;
					cv::Mat larvaROI(frame, 
							cv::Rect(blob->minx-ROI_PADDING-PAD,
								blob->miny-ROI_PADDING-PAD,
								blob->maxx-blob->minx+2*(ROI_PADDING+PAD),
								blob->maxy-blob->miny+2*(ROI_PADDING+PAD))
							);

					cv::circle(frame,
							cv::Point2d(blob->centroid.x,blob->centroid.y),
							1,
							cv::Scalar(255,0,0),
							-1);

					if(detected_larvae[it->first].isCluster==false)
					{
						cv::circle(larvaROI,
								cv::Point2d(
									detected_larvae[it->first].lrvskels.back().Point20.x+PAD,
									detected_larvae[it->first].lrvskels.back().Point20.y+PAD),
								1,
								cv::Scalar(255,255,0),
								-1);

						cv::circle(larvaROI,
								cv::Point2d(
									detected_larvae[it->first].heads.back().x+PAD,
									detected_larvae[it->first].heads.back().y+PAD),
								1,
								cv::Scalar(0,255,0),
								-1);

						cv::circle(larvaROI,
								cv::Point2d(
									detected_larvae[it->first].tails.back().x+PAD,
									detected_larvae[it->first].tails.back().y+PAD),
								1,
								cv::Scalar(0,0,255),
								-1);

						plotAngle(blob,larvaROI,PAD);
					}
					it++;
				}
			}

			cvReleaseImage(&labelImg);

		}
		//if (showimg)
		//{
		//cv::Mat flipframe;
		//cv::flip(frame,flipframe,0);
		//cv::imshow("Extracted Frame",flipframe);
		cv::imshow("Extracted Frame",frame);
		//cv::displayStatusBar("Extracted Frame","FOO",0);
		// showimg=false;
		//}
		//else
		//  showimg=true;

		int k;

		if (STEP==true)
			k=46;
		else
			k=cv::waitKey(1);
		if (k>=0)
		{
			if (k==32)
			{
				while (cv::waitKey(1)!=32){
					//No OP
				}
			}
			if (k==46)
			{
				STEP=true;
				while ((k=cv::waitKey(1))>=127 || k<0){
				}
				if (k!=46)
					STEP=false;
			}
			if (k==63234)
			{
				thresholdlow-=10;
				std::cerr << (int) thresholdlow << std::endl;
			}
			if (k==63235)
			{
				thresholdlow+=10;
				std::cerr << (int) thresholdlow << std::endl;
			}
			if (k==63232)
			{
				thresholdhigh+=10;
				std::cerr << (int) thresholdhigh << std::endl;
			}
			if (k==63233)
			{
				thresholdhigh-=10;
				std::cerr << (int) thresholdhigh << std::endl;
			}
			if (k==115)
			{
				SHOWTAGS=!SHOWTAGS;
			}
			if (k==116)
			{
				TRACK=!TRACK;
			}
			if (k==120)
			{
				exit(0);
			}
		}
	}

				summary.close();
}
