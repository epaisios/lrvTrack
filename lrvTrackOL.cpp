#include "lrvTrackOL.hpp"

/*
 * Function to draw spine points of larvaObject lrv on frame img
 */
void drawSpinePoints(Mat &img, larvaObject &lrv)
{
  if(lrv.lrvDistances.size()==0)
    return;
  vector<Point2f> const &spine=lrv.lrvDistances.back().Spine;
  if(spine.size()==0)
    return;
  vector<Point2f>::const_iterator it=spine.begin();
  for (;it!=spine.end();++it)
  {
    circle(img,
        *it,
        0,
        Scalar(255,0,255),
        -1);
  }
}

/*
 * Function to print brief contents of detected_larvae
 */
void dumpDetectedLarvae()
{
  cerr << "Size: " << detected_larvae.size() << endl;
  std::map<unsigned int,larvaObject>::const_iterator dlit;
  dlit=detected_larvae.begin();
  cerr << "Contents: " ;
  for(;dlit!=detected_larvae.end();++dlit)
  {
    cerr << dlit->first;
    cerr << " ";
    (dlit->second).dump();
  }
  cerr << endl;
}

/*
 * Quick function for returning the average of a vector
 */
double avgVec(vector<double> &vec)
{
  return (accumulate(vec.begin(),vec.end(),0)/vec.size());
}

/*
 * Quick function for returning the average of the last N
 * values of a vector.
 */
namespace std{
  template<typename data>
    double avgNVec(vector<data> &vec,unsigned int N=HISTORY_SIZE)
    {
      double SUM=0;
      unsigned int range;
      if (vec.size()>=N)
        range=N;
      else
        range=vec.size();

      return (accumulate(vec.rbegin(),vec.rbegin()+range,0)/range);
    }
}

/*
 * Simple factorial (we just need small numbers so this should be quick enough)
 */
unsigned int factorial(unsigned int num)
{
  unsigned int val=1;
  while(num>1)
    val=(num--)*val;

  return val;
}

/*
 * Number of combinations k of n
 */
unsigned int kofn(unsigned int k, unsigned int n)
{
  return factorial(n)/(factorial(k)*factorial(n-k));
}

/* This is not exactly the powerset. It is used for the detection of clustering
 * or divergence based on the combination that matches the data best.
 * We do not need the sets of 1 element and the complete SET in there
 * Input:
 *  IN: Vector of int (representing larvae IDs) to create the powersets from.
 *  OUT: Vector of vectors with the powersets (it does not contain powersets
 *       of 1 nor the complete vector as a set).
 */
void powersets(vector<unsigned int> &IN, vector<vector<unsigned int> > &OUT){
  for (unsigned int i=2 ; i<IN.size();i++)
  {
    vector<unsigned int> pointers;
    for(unsigned int k=0;k<i;k++)
    {
      pointers.push_back(k);
    }
    for (unsigned int j=0 ; j<kofn(i,IN.size());j++)
    {
      vector<unsigned int> cvec;
      for(unsigned int idx=0;idx<i;idx++)
      {
        cvec.push_back(IN[pointers[idx]]);
      }
      OUT.push_back(cvec);
      for(unsigned int inc=i;inc>0;inc--)
      {
        if(pointers[inc-1]<IN.size()-1-(i-inc))
        {
          pointers[inc-1]++;
          unsigned int add=0;
          for(unsigned int res=inc;res<i;res++)
          {
            add++;
            pointers[res]=pointers[inc-1]+add;
          }
          break;
        }
      }
    }
  }
}

/*
 * Function to see if the centre of blob1 matches the centre of blob2.
 * blob1 and blob2 are blobs coming from different frames and the 
 * function replies true if the blobs are close enough.
 * Input:
 *  blob1, blob2: the blobs to match
 *  val: the Manhattan distance between centres
 *  factor: defines the threshold for the match (factor*LARVA_OBJECT_LENGTH)
 */
bool centresMatch(
    cvb::CvBlob *blob1,
    cvb::CvBlob *blob2,
    double &val,
    double factor=LARVA_CENTRE_COMPARISON_FACTOR-1.0)
{
  stringstream DEBUG;
  double objectLength=
      min(max(blob1->maxx-blob1->minx,blob1->maxy-blob1->miny),
               max(blob2->maxx-blob2->minx,blob2->maxy-blob2->miny));

  DEBUG << "CentresMatchS ["<< blob1->label << ", " << blob2->label << "]: Length: " << objectLength << " Difx: " << fabs(blob1->centroid.x - blob2->centroid.x) << " Dify: " << fabs(blob1->centroid.y - blob2->centroid.y) << " Threshold: " << factor*LARVA_OBJECT_LENGTH;
  verbosePrint(DEBUG);

    val=fabs(blob1->centroid.x - blob2->centroid.x) + fabs(blob1->centroid.y - blob2->centroid.y) ;

  if (fabs(blob1->centroid.x - blob2->centroid.x) < factor*LARVA_OBJECT_LENGTH &&
      fabs(blob1->centroid.y - blob2->centroid.y)< factor*LARVA_OBJECT_LENGTH )
    {
      return true;
    }
  else
    {
      return false;
    }
}

/*
 * Function to see if the centres of blobs in "larvae" matches the centre of blob.
 * larvae and blob2 are derived from different frames and the 
 * function replies true if the barycentre if blobs in In and blob2
 * are close enough.
 * Input:
 *  In: Complete list of Blobs referenced by the larvae vector
 *  larvae: vector containing indices of larvae to consider for the comparison
 *          with blob. The indices refer to blobs in "In"
 *  blob2: the blob to match with larvae
 *  val: the Manhattan distance between barycentres
 *  factor: defines the threshold for the match (factor*LARVA_OBJECT_LENGTH)
 */
bool centresMatch(
    cvb::CvBlobs &In, 
    cvb::CvBlob *blob,
    vector<unsigned int> &larvae, 
    double &val,
    double factor=LARVA_CENTRE_COMPARISON_FACTOR-1)
{
  stringstream DEBUG;
  double xcomb=0, ycomb=0;
  double objectLength=
      max(blob->maxx-blob->minx,blob->maxy-blob->miny);
  double lrvAreaSum=0;
  if(larvae.size()==1)
  {
    xcomb=In[larvae[0]]->centroid.x;
    ycomb=In[larvae[0]]->centroid.y;
  }
  else
  {
    vector<unsigned int>::iterator it=larvae.begin();
    while(it!=larvae.end())
    {
      xcomb+=In[*it]->centroid.x*In[*it]->area;
      ycomb+=In[*it]->centroid.y*In[*it]->area;
      lrvAreaSum+=In[*it]->area;
      ++it;
    }
    xcomb=xcomb/lrvAreaSum;
    ycomb=ycomb/lrvAreaSum;
  }
  DEBUG << "CentresMatchP [" << blob->label << ", " << printVector(larvae) << "]: Length: " << objectLength << " Difx: " << fabs(blob->centroid.x - xcomb) << " Dify: " << fabs(blob->centroid.y - ycomb) << " Threshold: " << factor*LARVA_OBJECT_LENGTH;
  verbosePrint(DEBUG);
  val=fabs(blob->centroid.x - xcomb)+fabs(blob->centroid.y - ycomb);
  if (fabs(blob->centroid.x - xcomb)< factor* LARVA_OBJECT_LENGTH &&
      fabs(blob->centroid.y - ycomb)< factor* LARVA_OBJECT_LENGTH )
    {
      return true;
    }
  else
    {
      return false;
    }
}


/*
 * The function returns a guess on whether a blob is indeed a larva or not. The formula used
 * was optimized using optimization methods in R and real measured larvae from the Barcelona
 * data.
 *
 */
double is_larva(cvb::CvBlob *blob,bool out=false)
{
  vector<Point2f> newLarvaPoints;
  blobToPointVector(*blob,newLarvaPoints);
  larvaDistanceMap dstLarva(newLarvaPoints);
  fixContour(*blob,dstLarva,LRVTRACK_CONTOUR_RESOLUTION,colorFrame,previousFrame);
  Mat larvaROI;
  createLarvaContour(larvaROI,*blob);
  double greyVal=getGreyValue(larvaROI,*blob,greyFrame);
  double perimeter=getPerimeter(*blob);
  double a,l,w,p;
  a=blob->area;
  l=dstLarva.MaxDist;
  w=dstLarva.WidthDist;
  p=getPerimeter(*blob);
  cout << CURRENT_FRAME << ", " << blob->label << ", " << blob->area << " ," << l << ", " << w << ", " << greyVal << ", " << perimeter << endl;
  return 1.09*pow(a,0.4)
         -4.3*pow(l,0.8)
         +12*w
         +3*p; 
}


/*
 * Quick function to judge if the sizes of BLOB1 and BLOB2 are comparable.
 */
bool blobSizeIsRelevant(
    cvb::CvBlob *BLOB1,
    cvb::CvBlob *BLOB2,
    double ratio=LARVA_SIZE_COMPARISON_FACTOR)
{
  return (ratio*BLOB1->area > BLOB2->area &&
      ((2-ratio)*BLOB1->area < BLOB2->area));
}

/*
 * Quick function to judge if the sizes of blobs in larvae and BLOB are 
 * comparable.
 */
bool blobSizeIsRelevant(
    cvb::CvBlobs &In,
    cvb::CvBlob *BLOB,
    vector<unsigned int> &larvae,
    double ratio=LARVA_SIZE_COMPARISON_FACTOR)
{
  vector<unsigned int>::iterator IT=larvae.begin();
  double areaSUM=0;
  while (IT!=larvae.end())
  {
    areaSUM+=In[*IT]->area;
    ++IT;
  }
  return (ratio*BLOB->area > areaSUM &&
      ((2-ratio)*BLOB->area < areaSUM ));
}

//Returns all the larvae in the set Blobs within an area around "Blob".
// The area is determined by the size of the Blob (the longest side 
// of the bounding box multiplied by PADRatio).
//
// The vector nearbyLarvae is filled by those found sorted from closest
// to furthest.
void getNearbyLarvae(cvb::CvBlobs &Blobs, cvb::CvBlob *Blob, 
		            vector<unsigned int> &nearbyLarvae,bool pre=true,
                double PADRatio=1.5)
{
  vector<double> distances;
	double MaxDist = max(Blob->maxx-Blob->minx,Blob->maxy-Blob->miny);
	MaxDist=PADRatio*MaxDist;
  cvb::CvBlobs::iterator It=Blobs.begin();
  while (It!=Blobs.end())
	{
		cvb::CvBlob *cBlob=It->second;
    if (cBlob->centroid.x < (Blob->maxx + MaxDist/2) &&
        cBlob->centroid.x > (Blob->minx - MaxDist/2) &&
        cBlob->centroid.y < (Blob->maxy + MaxDist/2) &&
        cBlob->centroid.y > (Blob->miny - MaxDist/2) &&
        ( ((pre==false) && (assignedNew[It->first].size()<=0)) || 
          ((pre==true)  && (assignedPrevious[It->first].size()<=0)))
        )
    {
      double DIST=fast_abs(Blob->centroid.x - cBlob->centroid.x) +
        fast_abs(Blob->centroid.y - cBlob->centroid.y);
      vector<double>::iterator dIt=distances.begin();
      vector<unsigned int>::iterator nIt=nearbyLarvae.begin();
      if(nearbyLarvae.size()>0)
      {
        while(dIt!=distances.end())
        {
          if(DIST < *dIt)
          {
            distances.insert(dIt,DIST);
            nearbyLarvae.insert(nIt,It->first);
            break;
          }
          ++dIt;
          ++nIt;
        }
        if(dIt==distances.end())
        {
          nearbyLarvae.push_back(It->first);
          distances.push_back(DIST);
        }
      }
      else
      {
        distances.push_back(DIST);
        nearbyLarvae.push_back(It->first);
      }
    }
		++It;
	}
}

/***** Assigning Functions ********/

/*
 * These functions are just used as a quick way to set the correct IDs for
 * the assignedNew and assignedPrevious vectors used by the main tracking
 * algorithm. The assume the the IDs have already been correctly identified
 * and they do not change the global larvae structures (e.g. detected_larvae).
 */

/*
 * Function to assign an ID to a larva.
 * preID: is the ID of the larva from the previous frame to which preID
 *         matches
 * postID: is the ID of the blob of the current frame
 */
void assign_one(unsigned int preID,unsigned int postID)
{
  stringstream DEBUG;
  assignedPrevious[preID].push_back(postID);
  assignedNew[postID].push_back(preID);
  DEBUG << "Assigning: " << postID << " -> " << preID;
  verbosePrint(DEBUG);
  assignedPreMap[preID]=1;
}

/*
 * Function to match an ID to several IDs
 * preID: is the ID of the blob of the previous frame
 * postID: is the vector with the IDs of the larva in the new frame which
 *         diverged from preID
 */
void assign_one_OL(unsigned int preID,
                vector<unsigned int> postID)
{
  stringstream DEBUG;
  vector<unsigned int>::iterator postIT=postID.begin();
  assignedPreMap[preID]=postID.size();
  while(postIT!=postID.end())
  {
    unsigned int NewID=++LARVAE_COUNT;
    assignedNew[*postIT].push_back(NewID);
    DEBUG << "Assigning: " << *postIT << " -> " << NewID;
    verbosePrint(DEBUG);
    assignedPrevious[preID].push_back(NewID);
    detected_larvae[preID].childrenIDs.push_back(NewID);
    ++postIT;
  }
  DEBUG << "Cluster " << preID << " diverged into new larvae: " << printVector(assignedPrevious[preID]);
  verbosePrint(DEBUG);
  assignedPreMap[preID]=postID.size();
}

/*
 * Function to match an ID to several IDs
 * preID: is the ID of the blob of the previous frame
 * postID: is the vector with the IDs of the larva in the new frame which
 *         diverged from preID
 */
void assign_one(unsigned int preID,
                vector<unsigned int> postID)
{
  stringstream DEBUG;
  vector<unsigned int>::iterator postIT=postID.begin();
  assignedPreMap[preID]=postID.size();
  while(postIT!=postID.end())
  {
    unsigned int NewID=++LARVAE_COUNT;
    assignedNew[*postIT].push_back(NewID);
    DEBUG << "Assigning: " << *postIT << " -> " << NewID;
    verbosePrint(DEBUG);
    assignedPrevious[preID].push_back(NewID);
    ++postIT;
  }
  DEBUG << "Cluster " << preID << " diverged into new larvae: " << printVector(assignedPrevious[preID]);
  verbosePrint(DEBUG);
  assignedPreMap[preID]=postID.size();
}

/*
 * Function to match several IDs to an ID
 * preID: is the vector with the IDs of the larva in the old frame which
 *         collided to form blob with postID
 * postID: is the ID of the blob of the new frame
 */
void assign_one(vector<unsigned int> preID,
                unsigned int postID,unsigned int newID)
{
  stringstream DEBUG;
  vector<unsigned int>::iterator preIT=preID.begin();
  assignedNew[postID].push_back(newID);
  DEBUG << "New cluster " << postID << " from larvae: " << printVector(preID);
  verbosePrint(DEBUG);
  while(preIT!=preID.end())
  {
    assignedNew[postID].push_back(*preIT);
    assignedPrevious[*preIT].push_back(postID);
    assignedPreMap[*preIT]=1;
    ++preIT;
  }
}

/*
 * Function called whenever clustering is detected. Calls the
 * relevant assign one function and sets the correct details
 * in the detected_clusters vector. Assumes that correct IDs
 * have already been provided.
 * Input:
 *  POST_ID: the ID of the cluster
 *  IDs: the IDs of the blobs before the cluster
 */
void assign_clustering(
                      unsigned int POST_ID,
                      vector<unsigned int> &IDs
                      )
{
  unsigned int CLUSTER_ID=++LARVAE_COUNT;
  newClusters.push_back(CLUSTER_ID);
  vector<unsigned int> contents;
  contents.push_back(IDs.size());
  detected_clusters[CLUSTER_ID]=contents;
  map<unsigned int,vector<unsigned int> >::iterator dcIT;
  assign_one(IDs,POST_ID,CLUSTER_ID);
  vector<unsigned int>::iterator IT=IDs.begin();
  while (IT!=IDs.end())
  {
    if((dcIT=detected_clusters.find(*IT))!=detected_clusters.end())
    {
      detected_clusters[CLUSTER_ID].insert(
          detected_clusters[CLUSTER_ID].end(),
          dcIT->second.begin()+1,
          dcIT->second.end());
      // -2 because the cluster ID should be removed and replaced
      // by it's elements (-1) and another -1 because the 
      // detected_clusters contains as first element the size 
      // estimate of the cluster (the candidates can be more than that)
      detected_clusters[CLUSTER_ID][0]+=dcIT->second.size()-2;
      detected_clusters.erase(dcIT);
    }
    else
    {
      detected_clusters[CLUSTER_ID].push_back(*IT);
    }
    ++IT;
  }
}

/*
 * Function to calculate the Mahalanobis distance between
 * blob N and larva C.
 * Input:
 *  N: The ID of the Blob in the NEW CvBlobs structure which contains 
 *     the latest detected blobs.
 *  C: The ID of the larvae in the detected_larvae structure.
 * Output:
 *  the distance
 */
 double mh_dist(unsigned int N,unsigned int C)
{
  stringstream DEBUG;
  Mat candidateCovarMat;
  Mat candidateMeanMat;

  Mat newMeanMat;

  Mat Responses;
  Mat TrainArray;
  /*float size_avg=detected_larvae[C].area_sum/
    detected_larvae[C].area.size();
  float grey_value_avg=detected_larvae[C].grey_value_sum/
    detected_larvae[C].grey_value.size();
  float length_avg=detected_larvae[C].length_sum/
    detected_larvae[C].length.size();
  float perimeter_avg=detected_larvae[C].perimeter_sum/
    detected_larvae[C].perimeter.size();
  float width_avg=detected_larvae[C].width_sum/
    detected_larvae[C].width.size();*/

  float size_avg=avgNVec(detected_larvae[C].area);
  float grey_value_avg=avgNVec(detected_larvae[C].grey_value);
  float length_avg=avgNVec(detected_larvae[C].length);
  float perimeter_avg=avgNVec(detected_larvae[C].perimeter);
  float width_avg=avgNVec(detected_larvae[C].width);

  //float speed_x=avgNVec(detected_larvae[C].midpoint_speed_x);
  //float speed_y=avgNVec(detected_larvae[C].midpoint_speed_y);


  Mat InputArray;
  hconcat(Mat(detected_larvae[C].area),
      Mat(detected_larvae[C].grey_value),
      InputArray);

  hconcat(InputArray,
      Mat(detected_larvae[C].length),
      InputArray);

  hconcat(InputArray,
      Mat(detected_larvae[C].perimeter),
      InputArray);

  hconcat(InputArray,
      Mat(detected_larvae[C].width),
      InputArray);

  //hconcat(InputArray,
  //    Mat(detected_larvae[C].midpoint_speed_x),
  //    InputArray);

  //hconcat(InputArray,
  //    Mat(detected_larvae[C].midpoint_speed_y),
  //    InputArray);

  vector<float> meanVec;
  meanVec.push_back(size_avg);
  meanVec.push_back(grey_value_avg);
  meanVec.push_back(length_avg);
  meanVec.push_back(perimeter_avg);
  meanVec.push_back(width_avg);

  Mat(meanVec).copyTo(candidateMeanMat);
  Mat meanTMat;
  transpose(candidateMeanMat,meanTMat);


  calcCovarMatrix(InputArray, 
                  candidateCovarMat,
                  meanTMat,
                  CV_COVAR_ROWS|CV_COVAR_NORMAL|CV_COVAR_USE_AVG);
  candidateCovarMat.convertTo(candidateCovarMat,CV_32F);
  invert(candidateCovarMat,candidateCovarMat,DECOMP_SVD);

  Mat newSamplesMat;

  //Setup of new larva
  Mat larvaROI;
  vector<Point2f> newLarvaPoints;
  blobToPointVector(*NEW[N],newLarvaPoints);
  createLarvaContour(larvaROI,(*NEW[N]));
  larvaDistanceMap dstLarva(newLarvaPoints);
  Point2f centroid;
  centroid.x=NEW[N]->centroid.x;
  centroid.y=NEW[N]->centroid.y;
  //computeSpine(*NEW[N],dstLarva,grey_frame);
  fixContour(*NEW[N],
             dstLarva,
             LRVTRACK_CONTOUR_RESOLUTION,
             colorFrame,
             previousFrame);

  float newSize=NEW[N]->area;
  float newGreyval=getGreyValue(larvaROI,*NEW[N],greyFrame);
  float newLength=dstLarva.MaxDist;
  float newPerimeter=getPerimeter(*NEW[N]);
  float newWidth=dstLarva.WidthDist;

  meanVec.clear();
  meanVec.push_back(newSize);
  meanVec.push_back(newGreyval);
  meanVec.push_back(newLength);
  meanVec.push_back(newPerimeter);
  meanVec.push_back(newWidth);

  newMeanMat=Mat(meanVec);

 // cerr << candidateMeanMat << endl;
 // cerr << "==============================" << endl;
 // cerr << newMeanMat << endl;
 // cerr << "==============================" << endl;
 // cerr << candidateCovarMat << endl;

  double ret=Mahalanobis(candidateMeanMat,newMeanMat,candidateCovarMat);
  return ret;
}

/*
 * Function to ask whether assuming that blobP is assigned to blobN
 * is reasonable based on the speed calculated by this assignment.
 */
bool speedMatch(cvb::CvBlob &blobP, 
                cvb::CvBlob &blobN,
                double duration,
                double max_speed,
                double pre_speed_x=0.0,
                double pre_speed_y=0.0)
{
  double frames=(CURRENT_FRAME-detected_larvae[blobP.label].lastFrameWithStats);
  double mduration=frames/VIDEO_FPS;
  double uduration=max_speed*(1.0-(0.25-1/(frames+3)));
  double speedx = (blobP.centroid.x - blobN.centroid.x)/mduration;
  double speedy = (blobP.centroid.y - blobN.centroid.y)/mduration;
  
  speedx=speedx-pre_speed_x;
  speedy=speedy-pre_speed_y;

  double speed = sqrt(speedx*speedx + speedy*speedy);
  //cerr << "SpeedMatch: P:" << blobP.label << " N:" << blobN.label << 
  //  " M/S: " << uduration << "," << speed << endl;
  if(speed<uduration)
    return true;
  else
    return false;
}

/* Function to check if a Mapping makes sense. Currently uses speed
 * and a quick size check to figure this out. It may use other
 * options in the future.
 */
bool isMappingReasonable(lrvMapping &p,double duration)
{
      unsigned int lastIdx=detected_larvae[p.plrv].lastBlobWithStats;
      //dumpDetectedLarvae();
      cvb::CvBlob &blobP=detected_larvae[p.plrv].blobs[lastIdx];
      cvb::CvBlob &blobN=*NEW[p.nlrv];
      int csx=detected_larvae[p.plrv].centroid_speed_x.size();
      int csy=detected_larvae[p.plrv].centroid_speed_y.size();
      double pre_speed_x=avgNVec(detected_larvae[p.plrv].centroid_speed_x,3);
      double pre_speed_y=avgNVec(detected_larvae[p.plrv].centroid_speed_y,3);
      if(!speedMatch(
            blobP,
            blobN,
            duration,
            detected_larvae[p.plrv].max_centroid_speed,
            pre_speed_x,
            pre_speed_y) || 
          !blobSizeIsRelevant(&blobP,&blobN))
       return false;
      else
        return true;
}

/* Function to check if an assignment (several mappings) makes sense.
 * As in the above function we're using speed
 * and a quick size check to figure this out. It may use other
 * options in the future.
 */
bool isMappingReasonable(ltAssignments &m,double duration)
{
  for(unsigned int i=0;i<m.size();i++)
  {
      lrvMapping &p=m[i];
      unsigned int lastIdx=detected_larvae[p.plrv].lastBlobWithStats;
      cvb::CvBlob &blobP=detected_larvae[p.plrv].blobs[lastIdx];
      cvb::CvBlob &blobN=*NEW[p.nlrv];
      int csx=detected_larvae[p.plrv].centroid_speed_x.size();
      int csy=detected_larvae[p.plrv].centroid_speed_y.size();
      double pre_speed_x=avgNVec(detected_larvae[p.plrv].centroid_speed_x,3);
      double pre_speed_y=avgNVec(detected_larvae[p.plrv].centroid_speed_y,3);
      if(!speedMatch(
            blobP,
            blobN,
            duration,
            detected_larvae[p.plrv].max_centroid_speed) ||
          !blobSizeIsRelevant(&blobP,&blobN))
        return false;
  }
  return true;
}

/* Check if assignment contains a mapping of o->n
 */
bool assignmentContainsNewOld(ltAssignments &m,unsigned int n,unsigned int o)
{
  for(unsigned int i=0;i<m.size();i++)
  {
    if(m[i].nlrv==n || m[i].plrv==o)
      return true;
  }
  return false;
}

/* Given a vector of all possible mappings
 * returns the powerset of all possible assignments */
void pair_powersets(vector<lrvMapping> &IN, 
    vector<ltAssignments > &OUT)
{
  for (unsigned int i=1 ; i<=IN.size();i++)
  {
    vector<unsigned int> pointers;
    for(unsigned int k=0;k<i;k++)
    {
      pointers.push_back(k);
    }
    for (unsigned int j=0 ; j<kofn(i,IN.size());j++)
    {
      ltAssignments cvec;
      for(unsigned int idx=0;idx<i;idx++)
      {
        if(!assignmentContainsNewOld(cvec,
                                  IN[pointers[idx]].nlrv,
                                  IN[pointers[idx]].plrv))
          cvec.push_back(IN[pointers[idx]]);
      }
      if(cvec.size()==i)
        OUT.push_back(cvec);
      for(unsigned int inc=i;inc>0;inc--)
      {
        if(pointers[inc-1]<IN.size()-1-(i-inc))
        {
          pointers[inc-1]++;
          unsigned int add=0;
          for(unsigned int res=inc;res<i;res++)
          {
            add++;
            pointers[res]=pointers[inc-1]+add;
          }
          break;
        }
      }
    }
  }
}

//Adds a set of assignmets following a diverging event 
//to the possible assignments map.
void assignMapping(ltAssignments &m, 
                   map<unsigned int,unsigned int> &newAssignments)
{
  for(unsigned int i=0;i<m.size();i++)
  {
    newAssignments[m[i].nlrv]=m[i].plrv;
  }
}

//Return the accuracy of a set of mappings
double mappingAccuracy(ltAssignments &a)
{
  double SUM=0;
  for(unsigned int i=0;i<a.size();i++)
  {
    SUM+=a[i].getDistance();
  }
  return SUM/(1*(1+(a.size()-1)*1.5));
}

//Generate all the available reasonable mappings
void createReasonableMappings(
  vector<unsigned int> &candidateLarvae,
  vector<unsigned int> &newLarvae,
  double duration,
  vector<ltAssignments> &accepted_mappings)
{
    ltAssignments initial_pairs;
    vector<lrvMapping> valid_pairs;
    map<unsigned int,unsigned int> mappable_new;
    for(unsigned int i=0;i<candidateLarvae.size();i++)
    {
      for(unsigned int j=0;j<newLarvae.size();j++)
      {
        lrvMapping n(newLarvae[j],candidateLarvae[i]);
        //n.candidates=candidateLarvae;
        if(isMappingReasonable(n,duration))
        {
          n.setDistance();
          valid_pairs.push_back(n);
        }
      }
    }
  pair_powersets(valid_pairs,accepted_mappings);
}

/* This is not going to be used in the Off Line version.
 * Main function Handling diverging larvae
 * Input:
 *  candidateLarvae: vector with the IDs of the Larvae that
 *    were known to be involved in the collision 
 *  newLarvae: vector with the IDs of the unmapped blobs.
 * Output:
 *  newAssignments: mappints of the type [NEW]=OLD;
 */
void diverge_match(
  vector<unsigned int> &candidateLarvae,
  vector<unsigned int> &newLarvae,
  map<unsigned int, unsigned int> &newAssignments,
  double duration
  )
{
  vector<ltAssignments> valid_mappings;
  createReasonableMappings(
    candidateLarvae,
    newLarvae,
    duration,
    valid_mappings);

  map<double,unsigned int> mapping_accuracy;
  for(unsigned int i=0;i<valid_mappings.size();i++)
  {
    double acc=mappingAccuracy(valid_mappings[i]);
    //cerr << "=========================" << endl;
    //cerr << "Mapping " << i << ": " << endl;
    //for(unsigned int j=0;j<valid_mappings[i].size();j++)
    //{
    //  valid_mappings[i][j].print();
    //}
    //cerr << "Total: " << acc << endl;
    //cerr << "=========================" << endl;
    mapping_accuracy[acc]=i;
  }
  
  if(valid_mappings.size()!=0)
    assignMapping(valid_mappings[mapping_accuracy.begin()->second],
                newAssignments);
}


/*
 * For the Offline Version we will be treating every diverging case
 * as if the cluster were always new. The only difference would be
 * the registration of the fact that the new larvae derived from the
 * old larvae.
 */
void assign_divergingOL(cvb::CvBlobs &New,
                      unsigned int CLUSTER_ID,
                      vector<unsigned int> &IDs
                      )
{
  stringstream DEBUG;
  // We have the following cases here:
  //  1) CLUSTER_ID Corresponds to no cluster: This means
  //     the cluster is newly found (it started as a cluster)
  //  2) CLUSTER_ID Corresponds to an existing cluster:
  //     we assign all the new blobs as decendants of the previous blob
  map<unsigned int,vector<unsigned int> >::iterator dcIT;
  if((dcIT=detected_clusters.find(CLUSTER_ID))==detected_clusters.end())
  {
    //Not found new cluster NEW IDs to be given to the vector
    //elements
    DEBUG << "Cluster " << CLUSTER_ID << " is new. Assigning new IDs for diverged larvae";
    verbosePrint(DEBUG);
    assign_one(CLUSTER_ID,IDs);
    detected_larvae[CLUSTER_ID].isCluster=true;
  }
  else
  {
    assign_one_OL(CLUSTER_ID,IDs);
    detected_larvae[CLUSTER_ID].isCluster=true;
  }
}

void assign_diverging(cvb::CvBlobs &New,
                      unsigned int CLUSTER_ID,
                      vector<unsigned int> &IDs
                      )
{
  stringstream DEBUG;
  // We have the following cases here:
  //  1) CLUSTER_ID Corresponds to no cluster: This means
  //     the cluster is newly found (it started as a cluster)
  //  2) CLUSTER_ID Corresponds to an existing cluster:
  //      - We need to understand how it split.
  //        M -> M?
  //        M -> N<M?
  //        This will give us a clue about which ones are left.
  map<unsigned int,vector<unsigned int> >::iterator dcIT;
  if((dcIT=detected_clusters.find(CLUSTER_ID))==detected_clusters.end())
  {
    //Not found new cluster NEW IDs to be given to the vector
    //elements
    DEBUG << "Cluster " << CLUSTER_ID << 
      " is new. Assigning new IDs for diverged larvae";
    verbosePrint(DEBUG);
    assign_one(CLUSTER_ID,IDs);
    detected_larvae[CLUSTER_ID].isCluster=true;
  }
  else
  {
    vector<unsigned int> candidateLarvae(dcIT->second.begin()+1,
        dcIT->second.end());
    // New assignments are returned by the diverge_match_new function
    // [ NEW_ID -> OLD_ID ]
    map<unsigned int, unsigned int> newAssignments;
    /*
    cout << "Calling diverge match: CF: " << CURRENT_FRAME <<
      " BLOB SF: " << detected_larvae[dcIT->first].start_frame << endl;
      */
    unsigned int FRAMEDIFF = CURRENT_FRAME - 
                             detected_larvae[dcIT->first].start_frame;
    diverge_match(candidateLarvae,
        IDs,
        newAssignments,
        FRAMEDIFF/VIDEO_FPS
        );

    DEBUG << "New assignments returned by diverge matching: " << 
             endl << printUIMap(newAssignments);
    verbosePrint(DEBUG);
    // Perform the assignments

    //newCluster to register the assignments and what is left in the cluster
    vector<unsigned int> 
      newCluster(dcIT->second.begin()+1,dcIT->second.end()); 

    //Initiate a vector the the IDs to asign to keep track of what is assigned.
    vector<unsigned int> newIDs=IDs;

    for (map<unsigned int, unsigned int>::iterator naIT=
        newAssignments.begin();
        naIT!=newAssignments.end();
        ++naIT)
    {
      if(naIT->second!=0)
      {
        assign_one(naIT->second,naIT->first);
        assignedPreMap[CLUSTER_ID]=assignedPreMap[CLUSTER_ID]+1;
        assignedPrevious[CLUSTER_ID].push_back(naIT->second);

        for(vector<unsigned int>::iterator erIT=newCluster.begin();
            erIT!=newCluster.end();++erIT)
        {
          if(*erIT==naIT->second)
          {
            newCluster.erase(erIT);
            break;
          }
        }
        for(vector<unsigned int>::iterator erIT=newIDs.begin();
            erIT!=newIDs.end();++erIT)
        {
          if(*erIT==naIT->first)
          {
            newIDs.erase(erIT);
            break;
          }
        }
      }
    }

    // We have performed the assignents and what is left is contained in
    // the newIDs and the newCluster vectors
    if(newIDs.size()==newCluster.size() && newCluster.size()==0)
    {
      //perfect assignment. return from function
      detected_clusters.erase(dcIT);
      return;
    }
    if(newIDs.size()==newCluster.size() && newCluster.size()==1)
    {
      // Object diverged and left only one unassigned and 
      // only one is left in the cluster
      // We force the assignment
      assign_one(newCluster[0],newIDs[0]);
      assignedPreMap[CLUSTER_ID]=assignedPreMap[CLUSTER_ID]+1;
      assignedPrevious[CLUSTER_ID].push_back(newCluster[0]);
      detected_clusters.erase(dcIT);
      return;
    }
    if(newIDs.size()==1 && newCluster.size()>1)
    {
      //only one object remains to be assigned but our cluster was bigger
      unsigned int CLUSTER_ID=++LARVAE_COUNT;

      detected_clusters[CLUSTER_ID].push_back(newCluster.size());
      detected_clusters[CLUSTER_ID].insert(
          detected_clusters[CLUSTER_ID].end(),
          newCluster.begin(),
          newCluster.end());

      assign_one(newCluster,newIDs[0],CLUSTER_ID);
      detected_clusters.erase(dcIT);
    }
    //TODO: Figure out the case for newIDs.size()>1
    //  a) give new numbers to those that look like larvae
    //  b) create new clusters with dummy? larvae (if necessary i.e. if )
    //     plus those that were there
    //
    // IMPORTANT!! LIMIT THE SIZE
    if(newIDs.size()==2 && newCluster.size()==2)
    {
      if(is_larva(New[newIDs[0]])>1000 && is_larva(New[newIDs[1]])>1000)
      {
        //unsigned int CLUSTER_ID=++LARVAE_COUNT;
        //detected_clusters[CLUSTER_ID]=newCluster;
        //assign_one(newCluster,newIDs[0],CLUSTER_ID);
        detected_clusters.erase(dcIT);
      }
      if(is_larva(New[newIDs[0]])<1000 && is_larva(New[newIDs[1]])>1000)
      {
        // Unhandled
      }
    }
    // Unhandled should produce an message so that we debug it
    cerr << "Unhandled case were more than 1 IDs were unassigned"
      << endl;
  }
}


int detect_diverging(vector<unsigned int> &preLarvaeNearby,
                       vector<unsigned int> &newLarvaeNearby,
                       cvb::CvBlobs &Pre,
                       cvb::CvBlobs &New)
{
  verbosePrint("Trying to detect diverging clusters");

  stringstream DEBUG;
  DEBUG<< "Size of newLarvaeNearby: "  << newLarvaeNearby.size();
  verbosePrint(DEBUG);

  if(newLarvaeNearby.size()<=1)
    return -1; // No diverging is possible
  // Check the complete set first
  vector<unsigned int>::iterator pIT=preLarvaeNearby.begin();
  while(pIT!=preLarvaeNearby.end())
  {
    DEBUG << "Checking if nodes " << printVector(newLarvaeNearby) << " diverged from: " << *pIT;
    verbosePrint(DEBUG);
    double v;
    if(centresMatch(New,Pre[*pIT],newLarvaeNearby,v))
    {
      // Centres of all candidates match with new blob
      // cluster contains all. We can return
      DEBUG << "Node " << *pIT << " matches " << printVector(newLarvaeNearby);
      verbosePrint(DEBUG);
      assign_divergingOL(New,*pIT,newLarvaeNearby);
      break;
    }
    if(newLarvaeNearby.size()>2)
    {
      DEBUG << "Checking powersets of " << printVector(newLarvaeNearby) << " that diverged from: " << *pIT;
      verbosePrint(DEBUG);
      vector<vector<unsigned int> > pSETS;
      powersets(newLarvaeNearby,pSETS);
      vector<vector<unsigned int> >::iterator pSIT=pSETS.begin();
      while(pSIT!=pSETS.end())
      {
        if(centresMatch(New,Pre[*pIT],*pSIT,v))
        {
          // Centres of all candidates match with new blob
          // cluster contains all. We can return
          assign_divergingOL(New,*pIT,*pSIT);
          break;
        }
        ++pSIT;
      }
    }
    ++pIT;
  }
  return 0;
}

int detect_clustering(vector<unsigned int> &preLarvaeNearby,
                       vector<unsigned int> &newLarvaeNearby,
                       cvb::CvBlobs &Pre,
                       cvb::CvBlobs &New)
{
  verbosePrint("Trying to detect clusters");

  stringstream DEBUG;
  if(preLarvaeNearby.size()<=1)
    return -1; // No clustering is possible
  
  // Check the complete set first
  vector<unsigned int>::iterator nIT=newLarvaeNearby.begin();
  while(nIT!=newLarvaeNearby.end())
  {
    double v;
    if(centresMatch(Pre,New[*nIT],preLarvaeNearby,v) &&
       blobSizeIsRelevant(Pre,New[*nIT],preLarvaeNearby))
    {
      // Centres of all candidates match with new blob
      // cluster contains all. We can return
      DEBUG << "Centres of new larva: " << *nIT << " match centres of " 
       << printVector(preLarvaeNearby) << ". Assigning clustering";
      verbosePrint(DEBUG);

      assign_clustering(*nIT,preLarvaeNearby);
      //continue;
      break;
    }
    vector<vector<unsigned int> > pSETS;
    powersets(preLarvaeNearby,pSETS);
    vector<vector<unsigned int> >::iterator pSIT=pSETS.begin();
    while(pSIT!=pSETS.end())
    {
      verbosePrint("Trying to detect subsets for clustering");
      if(centresMatch(Pre,New[*nIT],*pSIT,v) &&
          blobSizeIsRelevant(Pre,New[*nIT],*pSIT))
      {
        // Centres of all candidates match with new blob
        // cluster contains subset pSIT. We can return
        DEBUG << "Centres of new larva: " << *nIT << " match centres of subset " 
          << printVector(*pSIT) << ". Assigning clustering";
        verbosePrint(DEBUG);
        assign_clustering(*nIT,*pSIT);
        break;
      }
      ++pSIT;
    }
    ++nIT;
  }
  return 0;
}

void findHeadTail(larvaObject &lrv,
                  Point2f &Head,
                  Point2f &Tail,
                  bool force_SurroundingValSearch=false)
{
  int max=0,min=65535;
  unsigned int i=0;
  cvb::CvBlob blob=lrv.blobs.back();
  vector<Point2f> &spine=lrv.lrvDistances.back().Spine;
  Point2f bp(lrv.blobs.back().minx,lrv.blobs.back().miny);
  Point2f sp_front=spine[0]-bp;
  Point2f sp_back=spine.back()-bp;
  if( lrv.start_frame==CURRENT_FRAME ||
      lrv.inCluster.back()>0 ||
      force_SurroundingValSearch ||
      lrv.roundness.back()<=2.0 ||
      (lrv.heads.back().x == 0 && lrv.heads.back().y==0 &&
       lrv.tails.back().x == 0 && lrv.tails.back().y==0)
    )
    {
      int BIAS=0;
      double curvatureBias=BIAS*lrv.lrvDistances.back().curvatureBias;
      double fsz=getSurroundingSize(sp_front,blob,colorFrame,previousFrame);
      double bsz=getSurroundingSize(sp_back,blob,colorFrame,previousFrame);

      if(curvatureBias==BIAS)
      {
        if(fsz-curvatureBias < bsz)
        {
          Tail=sp_back;
          Head=sp_front;
        }
        else
        {
          Tail=sp_front;
          Head=sp_back;
        }
      }
      else if(curvatureBias==0)
      {
        if(fsz > bsz-BIAS)
        {
          Tail=sp_front;
          Head=sp_back;
        }
        else
        {
          Tail=sp_back;
          Head=sp_front;
        }
      }
      else
        cerr << "PROBLEM AREAS AROUND HEAD AND TAIL ARE THE SAME" << endl;
    }
  else
    {
      double hfdiff=diff(lrv.heads.back(),sp_front);
      double hbdiff=diff(lrv.heads.back(),sp_back);
      double tfdiff=diff(lrv.tails.back(),sp_front);
      double tbdiff=diff(lrv.tails.back(),sp_back);
      if (hfdiff+tbdiff<hbdiff+tfdiff)
      {
        Head=sp_front;
        Tail=sp_back;
      }
      else
      {
        Head=sp_back;
        Tail=sp_front;
      }

    }
  double hwavg,twavg;
  double curvBias;
  if(Head==lrv.lrvDistances.back().Spine[0])
  {
    hwavg=lrv.lrvDistances.back().firstHalfWidthsSum;
    twavg=lrv.lrvDistances.back().secondHalfWidthsSum;
    curvBias=lrv.lrvDistances.back().curvatureBias;
  }
  else
  {
    hwavg=lrv.lrvDistances.back().secondHalfWidthsSum;
    twavg=lrv.lrvDistances.back().firstHalfWidthsSum;
    curvBias=1-lrv.lrvDistances.back().curvatureBias;
  }

  /*cerr << CURRENT_FRAME << " , " << 
               lrv.larva_ID <<  " , " <<
               getSurroundingSize(Tail,lrv.blobs.back(),origFrame,previousFrame)- 
               getSurroundingSize(Head,lrv.blobs.back(),origFrame,previousFrame) << " , " <<
               2*curvBias-1.0 << " , " <<
               twavg-hwavg << " , " <<
               lrv.roundness.back() <<
               endl;*/

}


void updateOneLarva(cvb::CvBlobs &In,
                    cvb::CvBlobs &Prev,
                    cvb::CvBlobs::iterator it,
                    tbb::concurrent_hash_map<unsigned int, larvaObject> &NEW_LARVA)
{
  unsigned int ID=(*it).first;
  cvb::CvBlob blob=*((*it).second);
  Mat larvaROI,cntPoints;
  createLarvaContour(larvaROI,blob);
  createLarvaContourPoints(cntPoints,blob);

  map<unsigned int,larvaObject>::iterator curLarva;
  // NEW LARVA OBJECT!
  if ((curLarva=detected_larvae.find(ID))==detected_larvae.end()
      ||
      Prev==In)
  {
    // Create and allocate the new object
    larvaObject newLarva;

    // Set the frame of it's existence
    newLarva.start_frame=CURRENT_FRAME;

    // Give the larva the necessary ID
    newLarva.larva_ID=ID;
    newLarva.old_ID=blob.n20;

    // Add the blob of the larva to its blob history
    newLarva.blobs.push_back(blob);

    // State that the larva is not in a blob
    newLarva.inCluster.push_back(false);

    Point2f centroid=Point2f(
        (blob.centroid.x-blob.minx+ROI_PADDING),
        (blob.centroid.y-blob.miny+ROI_PADDING)
        );

    Point2f centroidf=Point2f(
        (blob.centroid.x),
        (blob.centroid.y)
        );

    newLarva.capture_times.push_back(CURRENT_FRAME/VIDEO_FPS);

    // Initialize the speed to 0
    newLarva.midpoint_speed_x.push_back(0);
    newLarva.midpoint_speed_y.push_back(0);
    newLarva.max_midpoint_speed=0;

    newLarva.centroid_speed_x.push_back(0);
    newLarva.centroid_speed_y.push_back(0);
    newLarva.max_centroid_speed=0;

    newLarva.centroids.push_back(centroid);
    newLarva.centroids_full.push_back(centroidf);

    if (detected_clusters.find(ID)==detected_clusters.end())
    {

      ++newLarva.lifetimeWithStats;
      newLarva.lastBlobWithStats=0;

      newLarva.isCluster=false;

      //Initialize the area values
      newLarva.area.push_back(blob.area);
      newLarva.area_mean=blob.area;
      newLarva.area_sum=blob.area;
      newLarva.area_max=newLarva.area_min=blob.area;
      newLarva.area_min=newLarva.area_min=blob.area;

      double greyVal=getGreyValue(larvaROI,blob,greyFrame);
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
      newLarva.centroid_distance_x.push_back(0);
      newLarva.centroid_distance_y.push_back(0);
      newLarva.centroid_distance_x_sum=0;
      newLarva.centroid_distance_y_sum=0;

      newLarva.roundness.push_back((perimeter*perimeter)/(2*CV_PI*blob.area));
      cerr << newLarva.larva_ID << ": " << newLarva.roundness.back() << endl;
      // In this block we compute the inner spine for the larva
      vector<Point2f> cntPoints;
      blobToPointVector(blob,cntPoints);
      larvaDistanceMap Distances(cntPoints);
      //computeSpine(blob,Distances,frame);
      fixContour(blob,Distances,LRVTRACK_CONTOUR_RESOLUTION,colorFrame,previousFrame);
      newLarva.lrvDistances.push_back(Distances);
      newLarva.length.push_back(Distances.MaxDist);
      newLarva.length_mean = Distances.MaxDist;
      newLarva.length_sum= Distances.MaxDist;
      newLarva.length_max = Distances.MaxDist;
      newLarva.length_min = Distances.MaxDist;
      PointPair MAXPair=Distances.MaxDistPoints;

      Point2f Head,Tail;
      vector<Point2f> startPoints;
      /*startPoints.push_back(Point2f(
            MAXPair.first.x-blob.minx+ROI_PADDING,
            MAXPair.first.y-blob.miny+ROI_PADDING)
          );
      startPoints.push_back(Point2f(
            MAXPair.second.x-blob.minx+ROI_PADDING,
            MAXPair.second.y-blob.miny+ROI_PADDING)
          );*/
      newLarva.angular_speed.push_back(0);

      //findHeadTail(newLarva,Head,Tail);
      //findHeadTail(newLarva,Distances,Head,Tail);
      newLarva.heads.push_back(Head);
      newLarva.tails.push_back(Tail);

      Point2f MP;
      MP.x=Distances.MidPoint.x-newLarva.blobs.back().minx;
      MP.y=Distances.MidPoint.y-newLarva.blobs.back().miny;
      Point2f AxS(MP.x,Tail.y);
      newLarva.headBodyAngle.push_back(angleD(Tail,MP,Head));
      newLarva.orientationAngle.push_back(cvb::cvAngle(&blob));

      newLarva.width.push_back(Distances.WidthDist);
      newLarva.width_mean = Distances.WidthDist;
      newLarva.width_sum= Distances.WidthDist;
      newLarva.width_max = Distances.WidthDist;
      newLarva.width_min = Distances.WidthDist;


      if(LRVTRACK_CSV_OUTPUT)
      {
        csvfile << CURRENT_FRAME/VIDEO_FPS <<
          " , " << newLarva.larva_ID;
          vector<Point2f> &spine=newLarva.lrvDistances.back().Spine;
          if(newLarva.heads.back() == spine.back())
          {
            Point2f &t=spine[0];
            csvfile << " , " << (t.x-circles[bestCircle][0])*PIXEL_SIZE_IN_MM << 
                       " , " << (t.y-circles[bestCircle][1])*PIXEL_SIZE_IN_MM;
            vector<Point2f>::iterator it=spine.begin();
            it+=2;
            for(;it!=spine.end();it++)
            {
              //Point2f cp=*it-t;
              Point2f cp=*it;
              csvfile << " , " << (cp.x-circles[bestCircle][0])*PIXEL_SIZE_IN_MM << " , " 
                               << (cp.y-circles[bestCircle][1])*PIXEL_SIZE_IN_MM;
            }

          }
          else
          {
            Point2f &t=spine.back();
            csvfile << " , " << (t.x-circles[bestCircle][0])*PIXEL_SIZE_IN_MM << 
                       " , " << (t.y-circles[bestCircle][1])*PIXEL_SIZE_IN_MM;
            vector<Point2f>::reverse_iterator it=spine.rbegin();
            it+=2;
            for(;it!=spine.rend();it++)
            {
              //Point2f cp=*it-t;
              Point2f cp=*it;
              csvfile << " , " << (cp.x-circles[bestCircle][0])*PIXEL_SIZE_IN_MM << " , " 
                               << (cp.y-circles[bestCircle][1])*PIXEL_SIZE_IN_MM;
            }

          }
          //csvfile << " , " << newLarva.centroids.back().x + blob.minx <<
          //" , " << newLarva.centroids.back().y + blob.miny <<
          //" , " << newLarva.inCluster.back() <<
          csvfile << endl;
      }
      newLarva.lastFrameWithStats=CURRENT_FRAME;
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
      vector<unsigned int>::iterator cl;
      for ( cl=detected_clusters[ID].begin()+1 ; cl!=detected_clusters[ID].end() ; ++cl)
      {
        larvaObject &clusterLarva=detected_larvae[*cl];
        clusterLarva.inCluster.push_back(ID);
        clusterLarva.blobs.push_back(blob);
        /*
        if(LRVTRACK_CSV_OUTPUT)
        {
          csvfile << CURRENT_FRAME <<
            " , " << clusterLarva.larva_ID <<
            " , " << clusterLarva.area.back() <<
            " , " << clusterLarva.length.back() <<
            " , " << clusterLarva.width.back() <<
            " , " << clusterLarva.headBodyAngle.back() <<
            " , " << clusterLarva.orientationAngle.back() <<
            " , " << clusterLarva.midpoint_speed_x.back() <<
            " , " << clusterLarva.midpoint_speed_y.back() <<
            " , " << clusterLarva.centroids.size()-1 <<
            " , " << clusterLarva.centroids.back().x + blob.minx <<
            " , " << clusterLarva.centroids.back().y + blob.miny <<
            " , " << 0 <<
            " , " << 0 <<
            " , " << 0 <<
            " , " << 0 <<
            " , " << clusterLarva.inCluster.back() <<
            " , " << clusterLarva.isCluster <<
            endl;
        }
        */
      }

    }
    //detected_larvae[ID]=newLarva;
    //NEW[ID]=newLarva;
    tbb::concurrent_hash_map<unsigned int,larvaObject>::accessor a;
    NEW_LARVA.insert(a,ID);
    a->second=newLarva;
  }
  // UPDATED LARVA OBJECT!
  else
  {
    //Reference to current larva
    larvaObject &cur_larva=(*curLarva).second;
    //Pointer for the previous blob
    cvb::CvBlob &preBlob = cur_larva.blobs.back();

    // Set the ID of the larvaObject to the ID found TODO:Probably unnecessary
    cur_larva.larva_ID=ID;
    cur_larva.old_ID=blob.n20;

    // Add the current blob to the blobs history of the larva
    cur_larva.blobs.push_back(blob);

    // Create the skeleton of the larva and add it to the skeletons history
    Point2f centroid=Point2f(
        (blob.centroid.x-blob.minx+ROI_PADDING),
        (blob.centroid.y-blob.miny+ROI_PADDING)
        );

    Point2f centroidf=Point2f(
        (blob.centroid.x-blob.minx+ROI_PADDING),
        (blob.centroid.y-blob.miny+ROI_PADDING)
        );

    cur_larva.centroids_full.push_back(centroidf);
    cur_larva.centroids.push_back(centroid);

    cur_larva.capture_times.push_back(CURRENT_FRAME/VIDEO_FPS);

    // Look if larva is a blob
    if ( (!cur_larva.isCluster) &&
        (detected_clusters.find(ID)==detected_clusters.end()))
    {

      // If not then:
      //  Update area values for larva.

      ++cur_larva.lifetimeWithStats;
      cur_larva.lastBlobWithStats=cur_larva.blobs.size()-1;
      cur_larva.isCluster=false;

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
      //larvaSkel newLarvaSkel(larvaROI,centroid);
      //cur_larva.lrvskels.push_back(newLarvaSkel);

      cur_larva.centroid_distance_x.push_back(fabs(blob.centroid.x - preBlob.centroid.x));
      cur_larva.centroid_distance_y.push_back(fabs(blob.centroid.y - preBlob.centroid.y));
      cur_larva.centroid_distance_x_sum+=fabs(blob.centroid.x - preBlob.centroid.x);
      cur_larva.centroid_distance_x_sum+=fabs(blob.centroid.y - preBlob.centroid.y);


      // Point coordinates for head/tail
      Point2f Head,Tail;

      // Map to keep the distances of each point to all the others
      // Pair of points to keep the points with the Maximum distance (i.e. head and tail :) )
      vector<Point2f> cntPoints;
      blobToPointVector(blob,cntPoints);
      larvaDistanceMap Distances(cntPoints);
      // Compute all the inner distances for the larva
      //computeInnerDistances(blob,Distances,newLarvaSkel.MidPoint);
      //computeSpine(blob,Distances,frame);
      fixContour(blob,
                 Distances,
                 LRVTRACK_CONTOUR_RESOLUTION,
                 colorFrame,
                 previousFrame,
                 &cur_larva.heads,
                 &cur_larva.tails,
                 &cur_larva.blobs);
      double TIMEFRAME;
      if(cur_larva.inCluster.back()==0)
      {
        TIMEFRAME=1/VIDEO_FPS;
      }
      else
      {
        TIMEFRAME=(CURRENT_FRAME-cur_larva.lastFrameWithStats)/VIDEO_FPS;
      }
      cur_larva.midpoint_speed_x.push_back(
          (Distances.MidPoint.x - cur_larva.lrvDistances.back().MidPoint.x)
          /TIMEFRAME);
      cur_larva.midpoint_speed_y.push_back(
          (Distances.MidPoint.y - cur_larva.lrvDistances.back().MidPoint.y)
          /TIMEFRAME);

      cur_larva.centroid_speed_x.push_back(
          (blob.centroid.x - cur_larva.blobs[cur_larva.blobs.size()-2].centroid.x)
          /TIMEFRAME);
      cur_larva.centroid_speed_y.push_back(
          (blob.centroid.y - cur_larva.blobs[cur_larva.blobs.size()-2].centroid.y)
          /TIMEFRAME);
      
      double mspeed=sqrt(cur_larva.midpoint_speed_x.back()*
                      cur_larva.midpoint_speed_x.back() +
                      cur_larva.midpoint_speed_y.back()*
                      cur_larva.midpoint_speed_y.back());
      if(fabs(cur_larva.max_midpoint_speed) < fabs(mspeed))
      {
        cur_larva.max_midpoint_speed=mspeed;
      }

      double cspeed=sqrt(cur_larva.centroid_speed_x.back()*
                      cur_larva.centroid_speed_x.back() +
                      cur_larva.centroid_speed_y.back()*
                      cur_larva.centroid_speed_y.back());
      if(fabs(cur_larva.max_centroid_speed) < fabs(cspeed))
      {
        cur_larva.max_centroid_speed=cspeed;
      }

      cur_larva.lrvDistances.push_back(Distances);
      cur_larva.length.push_back(Distances.MaxDist);
      cur_larva.length_mean=(cur_larva.length_mean+Distances.MaxDist)/2;
      cur_larva.length_sum=cur_larva.length_sum+Distances.MaxDist;
      if (cur_larva.length_max < Distances.MaxDist)
      {
        cur_larva.length_max=Distances.MaxDist;
      }
      if (cur_larva.length_min > Distances.MaxDist)
      {
        cur_larva.length_min=Distances.MaxDist;
      }
      PointPair MAXPair=Distances.MaxDistPoints;
      
      cur_larva.width.push_back(Distances.WidthDist);
      cur_larva.width_mean=(cur_larva.width_mean+Distances.WidthDist)/2;
      cur_larva.width_sum=cur_larva.width_sum+Distances.WidthDist;
      if (cur_larva.width_max < Distances.WidthDist)
      {
        cur_larva.width_max=Distances.WidthDist;
      }
      if (cur_larva.width_min > Distances.WidthDist)
      {
        cur_larva.width_min=Distances.WidthDist;
      }

      double greyVal=getGreyValue(larvaROI,blob,greyFrame);
      cur_larva.grey_value.push_back(greyVal);
      cur_larva.grey_value_mean=(cur_larva.grey_value_mean+greyVal)/2;
      cur_larva.grey_value_sum=cur_larva.grey_value_sum+greyVal;
      if (cur_larva.grey_value_max < greyVal)
      {
        cur_larva.grey_value_max=greyVal;
      }
      if (cur_larva.grey_value_min > greyVal)
      {
        cur_larva.grey_value_min=greyVal;
      }

      double perimeter=getPerimeter(blob);
      cur_larva.perimeter.push_back(perimeter);
      cur_larva.perimeter_mean=(cur_larva.perimeter_mean+perimeter)/2;
      cur_larva.perimeter_sum=cur_larva.perimeter_sum+perimeter;
      if (cur_larva.perimeter_max < perimeter)
      {
        cur_larva.perimeter_max=perimeter;
      }
      if (cur_larva.perimeter_min > perimeter)
      {
        cur_larva.perimeter_min=perimeter;
      }

      cur_larva.roundness.push_back((perimeter*perimeter)/(2*CV_PI*blob.area));
      /*if(cur_larva.roundness.back() >= 2.7 || 
         cur_larva.length.back() > cur_larva.length_max/1.8)
        cerr << cur_larva.larva_ID << ": " << cur_larva.roundness.back() << endl;
      else
      {
        cerr<< "\033[1;31m";
        cerr <<"ROUND "<<cur_larva.larva_ID<<": "<< cur_larva.roundness.back();
        cerr << "\033[0m" << endl;
      }*/
      // Construct a vector of points including both
      // head and tail to decide which is which
      vector<Point2f> startPoints;

      // Points must be corrected to match the ROI including the padding set
      /*startPoints.push_back(Point2f(
            MAXPair.first.x-blob.minx+ROI_PADDING,
            MAXPair.first.y-blob.miny+ROI_PADDING)
          );
      startPoints.push_back(Point2f(
            MAXPair.second.x-blob.minx+ROI_PADDING,
            MAXPair.second.y-blob.miny+ROI_PADDING)
          );*/
      //execute the function to find which is which and assign appropriately
      //to Head/Tail.
      findHeadTail(cur_larva,Head,Tail);
      //findHeadTail(cur_larva,Distances,Head,Tail);

      //Add head and tail to the history
      cur_larva.heads.push_back(Head);
      cur_larva.tails.push_back(Tail);

      Point2f MP;
      MP.x=Distances.MidPoint.x-cur_larva.blobs.back().minx;
      MP.y=Distances.MidPoint.y-cur_larva.blobs.back().miny;
      cur_larva.headBodyAngle.push_back(angleD(Tail,MP,Head));
      Point2f AxS(MP.x,Tail.y);
      cur_larva.orientationAngle.push_back(cvb::cvAngle(&blob));

      double curAngle=cvb::cvAngle(&blob);
      double preAngle=cvb::cvAngle(&preBlob);

      cur_larva.angular_speed.push_back(fast_abs(curAngle-preAngle)*VIDEO_FPS);

      //state that larva is not detected as part of a blob of larvae
      //NOTE: It is important to perform this value setting !AFTER!
      //      searching for the head tail because, the head tail
      //      search behaves differently if the larva was previously
      //      part of a blob.
      cur_larva.inCluster.push_back(0);

      if(LRVTRACK_CSV_OUTPUT)
      {
        csvfile << CURRENT_FRAME/VIDEO_FPS <<
          " , " << cur_larva.larva_ID;
          vector<Point2f> &spine=cur_larva.lrvDistances.back().Spine;
          if(cur_larva.heads.back() == spine.back())
          {
            Point2f &t=spine[0];
            csvfile << " , " << (t.x-circles[bestCircle][0])*PIXEL_SIZE_IN_MM << 
                       " , " << (t.y-circles[bestCircle][1])*PIXEL_SIZE_IN_MM;
            vector<Point2f>::iterator it=spine.begin();
            it+=2;
            for(;it!=spine.end();it++)
            {
              //Point2f cp=*it-t;
              Point2f cp=*it;
              csvfile << " , " << (cp.x-circles[bestCircle][0])*PIXEL_SIZE_IN_MM << " , " 
                               << (cp.y-circles[bestCircle][1])*PIXEL_SIZE_IN_MM;
            }

          }
          else
          {
            Point2f &t=spine.back();
            csvfile << " , " << (t.x-circles[bestCircle][0])*PIXEL_SIZE_IN_MM << 
                       " , " << (t.y-circles[bestCircle][1])*PIXEL_SIZE_IN_MM;
            vector<Point2f>::reverse_iterator it=spine.rbegin();
            it+=2;
            for(;it!=spine.rend();it++)
            {
              //Point2f cp=*it-t;
              Point2f cp=*it;
              csvfile << " , " << (cp.x-circles[bestCircle][0])*PIXEL_SIZE_IN_MM << " , " 
                               << (cp.y-circles[bestCircle][1])*PIXEL_SIZE_IN_MM;
            }

          }
          //csvfile << " , " << cur_larva.centroids.back().x + blob.minx <<
          //" , " << cur_larva.centroids.back().y + blob.miny <<
          //" , " << cur_larva.inCluster.back() <<
          csvfile << endl;
      }
      cur_larva.lastFrameWithStats=CURRENT_FRAME;
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
      vector<unsigned int>::iterator cl;
      if(detected_clusters[ID].size()<2)
      {
        cerr << "We are working with a cluster but the size of registered nodes is too small" << endl;
        return;
      }
      
      for ( cl=detected_clusters[ID].begin()+1 ; cl!=detected_clusters[ID].end() ; ++cl)
      {
        larvaObject &clusterLarva=detected_larvae[*cl];
        clusterLarva.inCluster.push_back(ID);
        clusterLarva.blobs.push_back(blob);
        /*
        if(LRVTRACK_CSV_OUTPUT)
        {
          csvfile << CURRENT_FRAME <<
            " , " << clusterLarva.larva_ID <<
            " , " << clusterLarva.headBodyAngle.back() <<
            " , " << clusterLarva.orientationAngle.back() <<
            " , " << clusterLarva.midpoint_speed_x.back() <<
            " , " << clusterLarva.midpoint_speed_y.back() <<
            " , " << clusterLarva.centroids.size()-1 <<
            " , " << clusterLarva.centroids.back().x + blob.minx <<
            " , " << clusterLarva.centroids.back().y + blob.miny <<
            " , " << 0 <<
            " , " << 0 <<
            " , " << 0 <<
            " , " << 0 <<
            " , " << clusterLarva.inCluster.back() <<
            " , " << clusterLarva.isCluster <<
            endl;
        }
        */
      }
    }
  }
}

void bg_without_larvae(Mat &fgImg)
{
  IplImage *fflabelImg;

  cvb::CvBlobs blobs;
  IplImage ipl_thresholded = fgImg;

  fflabelImg=cvCreateImage(
      cvGetSize(&ipl_thresholded), IPL_DEPTH_LABEL, 1);

  cvLabel(&ipl_thresholded, fflabelImg, blobs);
  cvb::cvFilterByArea(blobs, 40, 1000);
  cvb::CvBlobs::iterator it=blobs.begin();
  double c,max=0.0;
  while (it!=blobs.end())
  {
    Mat ROI;
    cvb::CvBlob &blob=*(it->second);
    createLarvaROI(fgImg, ROI, blob);
    createLarvaContour(ROI,
                        blob,
                        CV_8UC1,
                        0,
                        true,
                        Scalar(0),
                        8);
    
    ++it;
  }
}

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
void process_frame(Mat &newFrame,
                  Mat &greyBgFrame,
                  Mat &showFrame,
                  Mat &processedFrame)
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

      if(!showFrame.empty())
      {
        // Add the circle as a marker to the showFrame
        circle(showFrame,
            Point2f(circles[bestCircle][0],circles[bestCircle][1]),
            int(circles[bestCircle][2]),
            Scalar(0,255,0),1);
      }
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

//Match lost larvae
unsigned int matchLostLarvae(unsigned int newLabel)
{
  for(vector<unsigned int>::iterator lIT=lost_larvae.begin();
      lIT!=lost_larvae.end();lIT++)
  {
      cvb::CvBlob &blobN=*NEW[newLabel];
      if(detected_larvae[*lIT].blobs.size()==0)
        return 0;
      cvb::CvBlob &blobP=detected_larvae[*lIT].blobs.back();
      double duration=(CURRENT_FRAME 
             - detected_larvae[*lIT].lastFrameWithStats)/VIDEO_FPS;
      if(speedMatch(
            blobP,
            blobN,
            duration,
            detected_larvae[*lIT].max_centroid_speed) && 
          blobSizeIsRelevant(&blobP,&blobN))
      {
        unsigned int ret=*lIT;
        lost_larvae.erase(lIT);
        return ret;
      }
  }
  return 0;
}

void updateLarvae(cvb::CvBlobs &In, cvb::CvBlobs &Prev)
{
  cvb::CvBlobs::iterator it=In.begin();
  tbb::concurrent_hash_map<unsigned int, larvaObject> NEW_LARVA;
  larvaeUpdateBody body(In,Prev,NEW_LARVA);
  //parallel_for_(Range(0, In.size()), body);
  parallel_for(BlockedRange(0, In.size()), body);
  
  tbb::concurrent_hash_map<unsigned int, larvaObject>::iterator nit=
    NEW_LARVA.begin();
  while (nit!=NEW_LARVA.end())
  {
    detected_larvae[nit->first]=nit->second;
    ++nit;
    }
}



//
// TODO:
// IMPORTANT NOTE!!
//  We have essentially three mappings:
//   1) The previous frame blobs -> new frame blobs
//   2) The new frame blobs -> previous frame blobs
//   3) The new frame blobs -> updated numbers
//
//   The updated numbers are:
//    a) the numbers that remained from the previous frame.
//    b) the numbers that diverged from frames before the previous
//    c) new numbers (increasing LARVAE_COUNT)
//
void newLarvaeTrack(cvb::CvBlobs &In, cvb::CvBlobs &Prev, cvb::CvBlobs &out)
{

  stringstream DEBUG;

  assignedNew.clear();
  assignedPrevious.clear();
  assignedPreMap.clear();
  newInFrame.clear();
  newClusters.clear();

  cvb::CvBlobs::iterator prevIt=Prev.begin();
  while (prevIt!=Prev.end())
  {
    assignedPreMap[prevIt->first]=0;
    ++prevIt;
  }

  verbosePrint("Starting Tracking loop");

	prevIt=Prev.begin();
  while (prevIt!=Prev.end())
  {
    //If the larvae from the previous frame has not been assigned
    //  -- It may be that it is assigned during a previous larva
    //     inspection.
    if(assignedPreMap[prevIt->first]==0)
    {
      cvb::CvBlob *preBlob=prevIt->second;
      unsigned int preID=prevIt->first;

      // Get larva around it before and after
      vector<unsigned int> preLarvaeNearby;
      vector<unsigned int> postLarvaeNearby;
      getNearbyLarvae(Prev,preBlob,preLarvaeNearby,true);
      getNearbyLarvae(In,preBlob,postLarvaeNearby,false);

      DEBUG << "Around larva " << prevIt->first << " we found " << preLarvaeNearby.size() << " larvae from the previous frame and " << postLarvaeNearby.size() << " from the new frame.";
      verbosePrint(DEBUG);

      DEBUG << "Around larva " << prevIt->first << " we found pre: " << printVector(preLarvaeNearby) << " and post: " << printVector(postLarvaeNearby);
      verbosePrint(DEBUG);
      // Great case collection now:
      
      // CASE 1 -> 1
      if(preLarvaeNearby.size()==1 && postLarvaeNearby.size()==1)
      {
        if(blobSizeIsRelevant(In[postLarvaeNearby[0]],preBlob))
        {
          double v;
          if(centresMatch(In,Prev[preLarvaeNearby[0]],postLarvaeNearby,v,0.3))
          {
            DEBUG << "Larva " << prevIt->first << " from previous frame matches " << postLarvaeNearby[0] << " from new frame. Assigning.";
            verbosePrint(DEBUG);
            assign_one(preID,postLarvaeNearby[0]);
          }
          else
          {
            // Centres do NOT match! Maybe it disappeared.
            // TODO: Check
            verbosePrint("1-1: with size similar\
                and centres not matching. !UNHANDLED!");
          }
        }
        else
        {
          // Sizes are too different:
          //  Cases:
          //    1) Object dissapeared and another came into sight.
          //       Case will be handled when the other one is handled
          //    2) Merge with a large object.
          //       Case will be handled when the large object is handled.
          verbosePrint("1-1: with size too different. No assignment.");
        }
      }
      // END OF 1-1 CASE
      // GENERAL CASE:
      else
      {
        // Try to assign obvious ones:
        vector<unsigned int>::iterator preNearIt=preLarvaeNearby.begin();
        while(preNearIt!=preLarvaeNearby.end())
        {
          bool erased=false;
          vector<unsigned int>::iterator postNearIt=postLarvaeNearby.begin();
          while(postNearIt!=postLarvaeNearby.end())
          {
            double val1,valb,vala;
            bool cm1,cmb,cma;

            cm1=centresMatch(Prev[*preNearIt],In[*postNearIt],val1,0.2);
            cmb=centresMatch(Prev,In[*postNearIt],preLarvaeNearby,valb);
            cma=centresMatch(In,Prev[*preNearIt],postLarvaeNearby,vala);

            if ( blobSizeIsRelevant(Prev[*preNearIt],In[*postNearIt]) &&
                 min(val1,valb)==val1 && 
                 min(val1,vala)==val1 && 
                 cm1)
            {
              //1-1 and one extra in sight
              //New Larva matches preID larva therefore assign and ignore the
              //other.
              DEBUG << "Larva " << *preNearIt << " from previous frame matches " << *postNearIt << " from new frame. Assigning";
              verbosePrint(DEBUG);
              assign_one(*preNearIt,*postNearIt);
              // erase the obvious assignments
              try{
                preLarvaeNearby.erase(preNearIt);
                postLarvaeNearby.erase(postNearIt);
              }
              catch(...)
              {
                cerr << "Problem erasing: " << *preNearIt << " or ";
                cerr << *postNearIt << endl;
              }
              erased=true;
            }
            else{
              DEBUG << "Larva " << *preNearIt << " from previous frame doesn't match " << *postNearIt << " from new frame.";
              verbosePrint(DEBUG);
              erased=false;
              ++postNearIt;
            }
          }
          if(!erased)
          {
            ++preNearIt;
            break;
          }
        }
        // The rest are either appearing/disappearing/clustering/diverging
        // BUT if the main one was indeed assigned then we are in luck 
        // and we move on
        DEBUG << "After initial assignment assignedPreMap: " << printUIMap(assignedPreMap) << " for larva with Id: " << preID;
        verbosePrint(DEBUG);
        if(assignedPreMap[preID]>0)
        {
          continue;
        }

        detect_clustering(preLarvaeNearby,postLarvaeNearby,Prev, In);

        DEBUG << "After clustering check assignedPreMap: " << printUIMap(assignedPreMap) << " for larva with Id: " << preID;
        verbosePrint(DEBUG);
        if(assignedPreMap[preID]>0)
        {
          continue;
        }

        detect_diverging(preLarvaeNearby,postLarvaeNearby,Prev, In);
        DEBUG << "After diverging check assignedPreMap: " << printUIMap(assignedPreMap) << " for larva with Id: " << preID;
        verbosePrint(DEBUG);
        if(assignedPreMap[preID]>0)
        {
          continue;
        }
        else
          verbosePrint("FOUND TESTCASE FOR MIXED DIVERGING/CONVERGING");
      }
    }
    ++prevIt;
    verbosePrint("assignedPrevious: ");
    map<unsigned int, vector<unsigned int> >::iterator prI=assignedPrevious.begin();
    while(prI!=assignedPrevious.end())
    {
      if(prI->second.size()>0)
      {
        DEBUG << prI->first << " -> " << printVector(prI->second);
        verbosePrint(DEBUG);
      }
      else
      {
        DEBUG << prI->first << " -> []" ;
        lost_larvae.push_back(prI->first);
        verbosePrint(DEBUG);
      }
      prI++;
    }

    verbosePrint("assignedNew: ");
    prI=assignedNew.begin();
    while(prI!=assignedNew.end())
    {
      if(prI->second.size()>0)
      {
        DEBUG << prI->first << " -> " << printVector(prI->second);
        verbosePrint(DEBUG);
      }
      else
      {
        DEBUG << prI->first << " -> []" ;
        
        verbosePrint(DEBUG);
      }
      prI++;
    }
  }


  cvb::CvBlobs::iterator bIT=In.begin();
  while(bIT!=In.end())
  {
    if (assignedNew[bIT->first].size()>0)
    {
      bIT->second->n20=bIT->second->label;
      bIT->second->label=assignedNew[bIT->first][0];
      out[assignedNew[bIT->first][0]]=bIT->second;
    }
    else
    {
      unsigned int Match=matchLostLarvae(bIT->first);
      if(Match==0)
      {
        unsigned int NEWID=++LARVAE_COUNT;
        verbosePrint("Extra Assignment:");
        DEBUG << bIT->first << " -> " << NEWID;
        newInFrame.push_back(NEWID);
        verbosePrint(DEBUG);
        bIT->second->n20=bIT->second->label;
        bIT->second->label=NEWID;
        out[NEWID]=bIT->second;
      }
      else
      {
        bIT->second->n20=bIT->second->label;
        bIT->second->label=Match;
        out[Match]=bIT->second;
      }
    }
    bIT++;
  }

  updateLarvae(out,Prev);
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
  previousFrame=greyFrame;
#ifdef LRVTRACK_WITH_OPENCL
  bool retval=capture.read(output);
  ocl::oclMat xchange;
  ocl::oclMat ocl_output(output);
  ocl::oclMat ocl_colorOutput;
  if(ocl_output.channels()==3)
  {
    ocl_colorOutput=ocl_output;
    ocl_output=xchange;
    ocl::cvtColor(ocl_colorOutput,ocl_output,CV_BGR2GRAY);
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
    cout << "HERE 1" << endl;
    ocl_output=ocl_output*MF+PF;
    cout << "HERE 2" << endl;
    ocl::addWeighted(ocl_output,1.0,ocl_output,-2.0,THRESHOLD,ocl_output);
    cout << "HERE 3" << endl;
    ocl_output.convertTo(ocl_output,CV_32F);
    ocl::pow(ocl_output,cpow,ocl_output);
    ocl_output.convertTo(ocl_output,CV_8UC1);
    ocl::cvtColor(ocl_output,ocl_colorOutput,CV_GRAY2BGR);
  }
  else
  {
    if(!xchange.empty())
      ocl::cvtColor(ocl_output,ocl_colorOutput,CV_GRAY2BGR);
  }
  output=ocl_output;
  colorOutput=ocl_colorOutput;
  return retval;
#else
  bool retval=capture.read(output);
  Mat xchange;
  if(output.channels()==3)
  {
    colorOutput=output;
    output=xchange;
    cvtColor(colorOutput,output,CV_BGR2GRAY);
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
    cvtColor(output,colorOutput,CV_GRAY2BGR);
  }
  else
  {
    if(!xchange.empty())
      cvtColor(output,colorOutput,CV_GRAY2BGR);
  }
  return retval;
#endif
}

//Function to extract the background
// The function retrieves a frame from the capture device and 
// returns the background into the bgFrame image.
//   Input:
//     * capture: the video capture device
//   Output:
//     * greyBgFrame: the Suggested background
void extract_background(VideoCapture &capture,
                       Mat &greyBgFrame)
{
#ifdef LRVTRACK_WITH_OPENCL
#else
  Mat origFrame;
  // Grab a frame from the stream
  if(!get_next_frame(capture,greyBgFrame,origFrame))
  {
    //TODO: Error handling
    exit(1);
  }
  greyBgFrame.copyTo(origFrame);
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
      double thresholdlow,thresholdhigh=0;
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

  //Initialize a first background frame. Everything that is in the circle
  //is black (i.e. will be excluded from the background). 
  //Outside the circle is white.
  if(circles.size()>0)
  {
    circle(greyBgFrame,
        Point2f(circles[bestCircle][0], circles[bestCircle][1]),
        circles[bestCircle][2]*0.9,
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
  Mat showFrame;
  process_frame(origFrame,
               greyBgFrame,
               showFrame,
               tempFrame);
  // We remove the larvae from the background
  bg_without_larvae(tempFrame);

  // We add the greyBgFrame with the tempFrame to get the full
  // background image
  add(greyBgFrame,tempFrame,greyBgFrame);
#endif

 return;
  
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
  bool STOP=false;

#ifdef LRVTRACK_WITH_CUDA
  gpu::printShortCudaDeviceInfo(gpu::getDevice());
#elif defined(LRVTRACK_WITH_OPENCL)
  ocl::DevicesInfo oclinfo;
  ocl::getOpenCLDevices(oclinfo);
  cout << "Found " << oclinfo.size() << " OpenCL devices. Using: ";
  cout << oclinfo[0]->deviceName << endl;
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
  extract_background(capture,bgFrame);
  cvb::CvBlobs preBlobs;
  while (!STOP)
  {
    if (!get_next_frame(capture,greyFrame,colorFrame))
      break;
    process_frame(greyFrame,bgFrame,colorFrame,processedFrame);
    IplImage ipl_thresholded = processedFrame;
    labelImg=cvCreateImage(
        cvGetSize(&ipl_thresholded), IPL_DEPTH_LABEL, 1);
    NEW.clear();
    cvLabel(&ipl_thresholded, labelImg, NEW);
    cvb::cvFilterByArea(NEW, 40, 1000);
    cvb::CvBlobs::iterator blIT=NEW.begin();
    cvb::CvBlobs tracked_blobs;
    cvb::CvBlobs blob1;
    cvb::CvBlobs::iterator it=NEW.begin();

    if(preBlobs.size()>0)
    {
      //larvae_track(blobs,preBlobs,tracked_blobs);
      newLarvaeTrack(NEW,preBlobs,tracked_blobs);

      if(LRVTRACK_CHOREOGRAPHY_OUTPUT)
      {
        //printSummary(preBlobs,NEW,false);
      }

      preBlobs=tracked_blobs;
    }
    else
    {
      preBlobs.clear();
      //normalize blobs for next use
      cvb::CvBlobs::iterator it=NEW.begin();
      int i=1;
      while (it!=NEW.end())
      {
        preBlobs[i]=(*it).second;
        preBlobs[i]->label=i;
        it++;
        i++;
      }
      updateLarvae(preBlobs,preBlobs);
      LARVAE_COUNT=preBlobs.size();
    }

    //dumpDetectedLarvae();
    if(LRVTRACK_CHOREOGRAPHY_OUTPUT)
    {
      //printSummary(preBlobs,NEW,true);
    }


    if (SHOWTAGS!=0)
    {
      it=tracked_blobs.begin();
      while (it!=tracked_blobs.end())
      {
        stringstream sstm;
        cvb::CvBlob *blob=it->second;
        if(LRVTRACK_CHOREOGRAPHY_OUTPUT)
        {
          //printBlobFile(detected_larvae[it->first]);
        }
        try
        {
          vector<unsigned int> &clusterMBs=
            detected_clusters.at(it->first);
          sstm << it->first << printVector(clusterMBs,1);
          /*sstm << "(" << detected_larvae[it->first].old_ID <<")";*/
          putText(colorFrame,
              sstm.str(),
              Point2f(blob->centroid.x+12,blob->centroid.y+12),
              FONT_HERSHEY_PLAIN,
              0.7,
              Scalar(255,255,255),
              1,
              CV_AA);
        }
        catch (const out_of_range& oor)
        {
          //out of range. i.e. it's not there :)
          sstm << (*it).first;
          /*sstm << "(" << detected_larvae[it->first].old_ID <<")";*/
          putText(colorFrame,
              sstm.str(),
              Point2f(blob->centroid.x+12,blob->centroid.y+12),
              FONT_HERSHEY_PLAIN,
              0.8,
              Scalar(255,255,255),
              1,
              CV_AA);
        }

        /*
           putText(frame,
           "*",
           Point2f(blob->centroid.x+12,blob->centroid.y-12),
           FONT_HERSHEY_PLAIN,
           0.8,
           Scalar(100,100,255),
           1,
           CV_AA);
           */
        // Here we define an extra PADDING
        int PAD=2;
        Mat larvaROI;
        try{
          larvaROI=Mat(colorFrame,
              Rect(blob->minx-ROI_PADDING-PAD,
                blob->miny-ROI_PADDING-PAD,
                blob->maxx-blob->minx+1+2*(ROI_PADDING+PAD),
                blob->maxy-blob->miny+1+2*(ROI_PADDING+PAD))
              );
        }
        catch(...)
        {
          cerr << "larvaROI failed. continuing" << endl;
          break;
        }
        Point2f cr(blob->centroid.x-blob->minx,blob->centroid.y-blob->miny);
        //larvaSkel testSkel(larvaROI,cr);
        //testSkel.drawSkeleton(larvaROI,Scalar(0,0,255));
        //imshow("ROI",larvaROI);
        //waitKey(1);
        //if(is_larva(blob)<IS_LARVA_THRESHOLD)
        //{
        //createLarvaContour(larvaROI,*blob,CV_8UC3,0,false,
        //    Scalar(0,255,0),8);
        drawSpinePoints(colorFrame,detected_larvae[it->first]);
        /*
           circle(frame,
           Point2f(blob->centroid.x,blob->centroid.y),
           0,
           Scalar(255,0,0),
           -1);
           */
        //}

        if(detected_larvae[it->first].isCluster==false)
        {
          circle(larvaROI,
              Point2f(
                detected_larvae[it->first].heads.back().x+PAD,
                detected_larvae[it->first].heads.back().y+PAD),
              1,
              Scalar(0,255,0),
              -1);

          circle(larvaROI,
              Point2f(
                detected_larvae[it->first].tails.back().x+PAD,
                detected_larvae[it->first].tails.back().y+PAD),
              1,
              Scalar(0,0,255),
              -1);

          //plotAngle(blob,larvaROI,PAD);
        }
        it++;
      }
    }

    cvReleaseImage(&labelImg);
    imshow("Extracted Frame",colorFrame);
    int k;

    if (STEP==true)
      k=46;
    else
      k=cv::waitKey(1);
    if (k>=0)
    {
      if (k==' ')
      {
        if(!TRACK)
          while (cv::waitKey(1)!=' ')
          {
            //No OP
          }
      }
      if (k=='.')
      {
        STEP=true;
        while ((k=cv::waitKey(1))>=127 || k<0)
        {
        }
        if (k!='.')
          STEP=false;
      }
      if (k=='s')
      {
        SHOWTAGS=!SHOWTAGS;
      }
      if (k=='t')
      {
        TRACK=!TRACK;
      }
      if (k=='x')
      {
        exit(0);
      }
    }
  }

}
