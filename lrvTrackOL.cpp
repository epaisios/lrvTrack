#include "lrvTrackOL.hpp"
#include "lrvTrackPartitionGenerator.hpp"

/*
 * Function to draw spine points of larvaObject lrv on frame img
 */
void drawSpinePoints(Mat &img, larvaObject &lrv,size_t idx)
{
  if(lrv.lrvDistances.size()==0)
    return;
  if(idx>lrv.lrvDistances.size())
    return;
  vector<Point2f> const &spine=lrv.lrvDistances[idx].Spine;
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

void drawSpinePoints(Mat &img, larvaObject &lrv)
{
  drawSpinePoints(img,lrv,lrv.lrvDistances.size()-1);
}

/*
 * Function to print brief contents of detected_larvae
 */
void dumpDetectedLarvae()
{
  BOOST_LOG_TRIVIAL(trace) << "Size: " << detected_larvae.size();
  std::map<size_t,larvaObject>::const_iterator dlit;
  dlit=detected_larvae.begin();
  stringstream d;
  d << "Contents: ";
  for(;dlit!=detected_larvae.end();++dlit)
  {
    d << dlit->first;
    d << " ";
    (dlit->second).dump();
  }
  BOOST_LOG_TRIVIAL(trace) << d.str() << endl;
}

void detectHeadTail()
{
  // Loop
}

/*void detectHeadTailLarva(larvaObject &l)
{
  double distance;
  for(size_t i=0;i<l.centroid_speed_x.size();i++)
  {
    Point2f distz=l.lrvDistances[i].Spine[0] - l.lrvDistances[i].MidPoint;
    Point2f distb=l.lrvDistances[i].Spine.back() - l.lrvDistances[i].MidPoint;
    Point2f mspeed(l.midpoint_speed_x,l.midpoint_speed_y);
  }
}*/

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
    double avgNVec(vector<data> &vec,size_t N=HISTORY_SIZE)
    {
      double SUM=0;
      size_t range;
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
uint64_t factorial(size_t n)
{
      return factorial_vec[n];
}

/*
 * Number of combinations k of n
 */
size_t kofn(size_t k, size_t n)
{
  return factorial(n)/(factorial(k)*factorial(n-k));
}

/*
 * Stirling number of the second kind
 * Set of N elements partitioned into k-blocks
 */

size_t stirling_2(size_t n, size_t k)
{
  int sum=0;
  for(size_t j=0;j<=k;++j)
    sum+=(int)(2*(((int) j+1)%2)-1) * (int)(kofn(j,k)*pow(k-j,n));
  return sum/factorial(k);
}

/* This is not exactly the powerset. It is used for the detection of clustering
 * or divergence based on the combination that matches the data best.
 * We do not need the sets of 1 element and the complete SET in there
 * Input:
 *  IN: Vector of int (representing larvae IDs) to create the powersets from.
 *  OUT: Vector of vectors with the powersets (it does not contain powersets
 *       of 1 nor the complete vector as a set).
 */
void powersets(vector<size_t> &IN, vector<vector<size_t> > &OUT){
  for (size_t i=2 ; i<IN.size();i++)
  {
    vector<size_t> pointers;
    for(size_t k=0;k<i;k++)
    {
      pointers.push_back(k);
    }
    for (size_t j=0 ; j<kofn(i,IN.size());j++)
    {
      vector<size_t> cvec;
      for(size_t idx=0;idx<i;idx++)
      {
        cvec.push_back(IN[pointers[idx]]);
      }
      OUT.push_back(cvec);
      for(size_t inc=i;inc>0;inc--)
      {
        if(pointers[inc-1]<IN.size()-1-(i-inc))
        {
          pointers[inc-1]++;
          size_t add=0;
          for(size_t res=inc;res<i;res++)
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
 * Function to handle various keystrokes for each frame
 */
void handleKeys(bool &STEP,bool &STOP, bool &TRACK, bool &SHOWTAGS)
{
    int k;
    if (STEP==true)
      k=46;
    else
      k=waitKey(1);
    if (k>=0)
    {
      if (k==' ')
      {
        if(!TRACK)
          while (waitKey(1)!=' ')
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

bool check_roundness(size_t lID, 
                     double spDist, 
                     double curRnd,
                     double RndRat,
                     double prespDist,
                     double variance)
{
      /*BOOST_LOG_TRIVIAL(debug) << "RND, S, " << CURRENT_FRAME 
        << ", " << lID
        << ", " << curRnd
        << ", " << RndRat
        << ", " << spDist
        << ", " << prespDist
        << ", " << variance;*/
    if(curRnd<=2.75)
    {
      return true;
    }
    else
      return false;
}

double is_larva(larvaObject &f)
{
  double avg_a=0.0;
  double avg_l=0.0;
  double avg_w=0.0;
  double avg_p=0.0;
  double avg_wavg=0.0;
  double avg_wsd=0.0;
  double total=0.0;
  size_t blen=f.blobs.size();
  for(size_t i=0;i<f.blobs.size();i++)
  {
    avg_a+=f.blobs[i].area;
    avg_l+=f.lrvDistances[i].MaxDist;
    avg_w+=f.lrvDistances[i].WidthDist;
    avg_p+=f.perimeter[i];
    avg_wavg+=f.lrvDistances[i].WidthAvg;
    avg_wsd+=f.lrvDistances[i].WidthSD;
  }
    avg_a=avg_a/blen;
    avg_l=avg_l/blen;
    avg_w=avg_w/blen;
    avg_p=avg_p/blen;
    avg_wavg=avg_wavg/blen;
    avg_wsd=avg_wsd/blen;
    total+=1.09*pow(avg_a,0.4)
         -4.3*pow(avg_l,0.8)
         +12*avg_w
         +3*avg_p;
    BOOST_LOG_TRIVIAL(debug) << "ISLARVA, " << f.larva_ID << ", "
                      << avg_a << ", " << avg_l << ", "
                      << avg_w << ", " << avg_p << ", " 
                      << avg_wavg << ", " << avg_wsd << ", " << total;
  return total;
}

void apply_model()
{
      for( auto &OOIp: detected_larvae) //OOI : Object of Interest
      {
        larvaObject &OOI=OOIp.second;
        if(OOI.start_frame==CURRENT_FRAME && 
           parent_blobs[OOI.larva_ID].size()>1 )
        {
          //Find predecessor's (if any and if collision)
          //and then use the model :)
          vector<cvb::CvBlob> vec;
          for(auto &parents: parent_blobs[OOI.larva_ID])
          {
            cvb::CvBlob &nblob=detected_larvae[parents].blobs.back();
            vec.push_back(nblob);
          }
          for(auto &parent: parent_blobs[OOI.larva_ID])
          {
            larvaFit lf;
            larvaObject &o = detected_larvae[parent];
            lf.setup(o);
            lf.setloops(LRVTRACK_DSTEP,LRVTRACK_WSTEP);
            lf.optimize(vec, OOI.blobs[0]);
          }
        }
      }

}
void showTags2()
{
      map<size_t,larvaObject>::iterator it=detected_larvae.begin();
      while (it!=detected_larvae.end())
      {
        if(it->second.start_frame>CURRENT_FRAME || 
           it->second.end_frame<CURRENT_FRAME)
        {
          it++;
          continue;
        }
        stringstream sstm;
        size_t c_index=CURRENT_FRAME-it->second.start_frame;
        cvb::CvBlob* blob=
                &(it->second.blobs[c_index]);
        if(LRVTRACK_CHOREOGRAPHY_OUTPUT)
        {
          //printBlobFile(detected_larvae[it->first]);
        }
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
          BOOST_LOG_TRIVIAL(warning) << "larvaROI failed. continuing";
          break;
        }
        Point2f cr(blob->centroid.x-blob->minx,blob->centroid.y-blob->miny);
        //larvaSkel testSkel(larvaROI,cr);
        //testSkel.drawSkeleton(larvaROI,Scalar(0,0,255));
        //imshow("ROI",larvaROI);
        //waitKey(1);
        //if(is_larva(blob)<IS_LARVA_THRESHOLD)
        //{
        /*
           circle(frame,
           Point2f(blob->centroid.x,blob->centroid.y),
           0,
           Scalar(255,0,0),
           -1);
           */
        //}

        if(it->second.isCluster==false)
        {
          circle(larvaROI,
              Point2f(
                it->second.heads[c_index].x+PAD,
                it->second.heads[c_index].y+PAD),
              2,
              Scalar(255,255,0),
              -1);

          circle(larvaROI,
              Point2f(
                it->second.tails[c_index].x+PAD,
                it->second.tails[c_index].y+PAD),
              2,
              Scalar(0,0,255),
              -1);

          createLarvaContour(larvaROI,*blob,CV_8UC3,PAD,false,
              Scalar(0,255,0),8);
          drawSpinePoints(colorFrame,it->second,c_index);
          BOOST_LOG_TRIVIAL(debug) << "LRVINFO:" << 
            it->second.larva_ID << ", " <<
            it->second.area[c_index] << ", " <<
            it->second.length[c_index] << ", " <<
            it->second.width[c_index] << ", " <<
            it->second.grey_value[c_index] << ", " <<
            it->second.perimeter[c_index];
          //plotAngle(blob,larvaROI,PAD);
        }
        it++;
      }
}

void showTags(const cvb::CvBlobs &tracked_blobs)
{
      cvb::CvBlobs::const_iterator it=tracked_blobs.begin();
      while (it!=tracked_blobs.end())
      {
        stringstream sstm;
        cvb::CvBlob *blob=it->second;
        if(LRVTRACK_CHOREOGRAPHY_OUTPUT)
        {
          //printBlobFile(detected_larvae[it->first]);
        }
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
          BOOST_LOG_TRIVIAL(warning) << "larvaROI failed. continuing";
          break;
        }
        Point2f cr(blob->centroid.x-blob->minx,blob->centroid.y-blob->miny);
        //larvaSkel testSkel(larvaROI,cr);
        //testSkel.drawSkeleton(larvaROI,Scalar(0,0,255));
        //imshow("ROI",larvaROI);
        //waitKey(1);
        //if(is_larva(blob)<IS_LARVA_THRESHOLD)
        //{
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
              Scalar(255,0,0),
              -1);

          circle(larvaROI,
              Point2f(
                detected_larvae[it->first].tails.back().x+PAD,
                detected_larvae[it->first].tails.back().y+PAD),
              1,
              Scalar(0,0,255),
              -1);

          createLarvaContour(larvaROI,*blob,CV_8UC3,PAD,false,
              Scalar(0,255,0),8);
          drawSpinePoints(colorFrame,detected_larvae[it->first]);

          //plotAngle(blob,larvaROI,PAD);
        }
        it++;
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
  double objectLength=
      min(max(blob1->maxx-blob1->minx,blob1->maxy-blob1->miny),
               max(blob2->maxx-blob2->minx,blob2->maxy-blob2->miny));

  BOOST_LOG_TRIVIAL(trace) << "CentresMatchS ["<< blob1->label << ", "
                           << blob2->label << "]: Length: " << objectLength 
                           << " Difx: " 
                           << fabs(blob1->centroid.x - blob2->centroid.x) 
                           << " Dify: " 
                           << fabs(blob1->centroid.y - blob2->centroid.y) 
                           << " Threshold: " 
                           << factor*LARVA_OBJECT_LENGTH;

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
    vector<size_t> &larvae, 
    double &val,
    double factor=LARVA_CENTRE_COMPARISON_FACTOR-1)
{
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
    vector<size_t>::iterator it=larvae.begin();
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
  BOOST_LOG_TRIVIAL(trace) << "CentresMatchP [" << blob->label 
                           << ", " << printVector(larvae) << "]: Length: " 
                           << objectLength << " Difx: " 
                           << fabs(blob->centroid.x - xcomb) 
                           << " Dify: " << fabs(blob->centroid.y - ycomb) 
                           << " Threshold: " << factor*LARVA_OBJECT_LENGTH;
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


bool larvaToRing(cvb::CvBlob &blob)
{
  cv::Point2f cntr(circles[bestCircle][0],circles[bestCircle][1]);
  if(p2fdist(blob.centroid.x,
             blob.centroid.y,
             circles[bestCircle][0],
             circles[bestCircle][1])
      >
      circles[bestCircle][2]*0.95)
  {
    std::vector<cv::Point2f> p;
    blobToPointVector(blob,p);
    for(size_t i=0;i<p.size();i++)
    {
      Point2f cp=p[i];
      double dst=p2fdist(p[i],cntr);
      if(dst>=circles[bestCircle][2]-2.0)
        return true;
    }
  }
  return false;
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
    vector<size_t> &larvae,
    double ratio=LARVA_SIZE_COMPARISON_FACTOR)
{
  vector<size_t>::iterator IT=larvae.begin();
  double areaSUM=0;
  while (IT!=larvae.end())
  {
    areaSUM+=In[*IT]->area;
    ++IT;
  }
  return (ratio*BLOB->area > areaSUM &&
      ((2-ratio)*BLOB->area < areaSUM ));
}

/* 
 * Returns all the larvae in the set Blobs within an area around "Blob".
 * The area is determined by the size of the Blob (the longest side 
 * of the bounding box multiplied by PADRatio).
 *
 * The vector nearbyLarvae is filled by those found sorted from closest
 * to furthest.
 */
void getNearbyLarvae(cvb::CvBlobs &Blobs, cvb::CvBlob *Blob, 
		            vector<size_t> &nearbyLarvae,bool pre=true,
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
      vector<size_t>::iterator nIt=nearbyLarvae.begin();
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

/***** Assigning Functions ********

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
void assign_one(size_t preID,size_t postID)
{
  assignedPrevious[preID].push_back(postID);
  assignedNew[postID].push_back(preID);
  BOOST_LOG_TRIVIAL(trace) << "Assigning: " << postID << " -> " << preID;
  assignedPreMap[preID]=1;
}

/*
 * Function to match an ID to several IDs
 * preID: is the ID of the blob of the previous frame
 * postID: is the vector with the IDs of the larva in the new frame which
 *         diverged from preID
 */
void assign_one(size_t preID,
                vector<size_t> postID)
{
  vector<size_t>::iterator postIT=postID.begin();
  assignedPreMap[preID]=postID.size();
  //detected_larvae[preID].isCluster=true;
  while(postIT!=postID.end())
  {
    size_t NewID=++LARVAE_COUNT;
    assignedNew[*postIT].push_back(NewID);
    BOOST_LOG_TRIVIAL(trace) << "Assigning: " << *postIT << " -> " << NewID;
    assignedPrevious[preID].push_back(NewID);
    detected_larvae[preID].diverged_to.push_back(NewID);
    parent_blobs[NewID].push_back(preID);
    children_blobs[preID].push_back(NewID);;
    ++postIT;
  }
  BOOST_LOG_TRIVIAL(debug) << CURRENT_FRAME << ", "
                           <<  preID 
                           << " -> "
                           << printVector(assignedPrevious[preID]);
  assignedPreMap[preID]=postID.size();
}

/*
 * Function to match several IDs to an ID
 * preID: is the vector with the IDs of the larva in the old frame which
 *         collided to form blob with postID
 * postID: is the ID of the blob of the new frame
 */
void assign_one(vector<size_t> preID,
                size_t postID,size_t newID)
{
  vector<size_t>::iterator preIT=preID.begin();
  assignedNew[postID].push_back(newID);
  BOOST_LOG_TRIVIAL(debug) << CURRENT_FRAME << ", "
                           << printVector(preID) << " -> "
                           << newID;
  while(preIT!=preID.end())
  {
    assignedNew[postID].push_back(*preIT);
    assignedPrevious[*preIT].push_back(postID);
    assignedPreMap[*preIT]=1;
    detected_larvae[*preIT].converged_to=newID;
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
                      size_t POST_ID,
                      vector<size_t> &IDs
                      )
{
  size_t CLUSTER_ID=++LARVAE_COUNT;
  //newClusters.push_back(CLUSTER_ID);
  //detected_converged[CLUSTER_ID]=IDs;
  parent_blobs[CLUSTER_ID]=IDs;
  for(auto &i:IDs)
    children_blobs[i].push_back(CLUSTER_ID);
  assign_one(IDs,POST_ID,CLUSTER_ID);
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
 double mh_dist(size_t N,size_t C)
{
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
      size_t lastIdx=detected_larvae[p.plrv].lastBlobWithStats;
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
  for(size_t i=0;i<m.size();i++)
  {
      lrvMapping &p=m[i];
      size_t lastIdx=detected_larvae[p.plrv].lastBlobWithStats;
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
bool assignmentContainsNewOld(ltAssignments &m,size_t n,size_t o)
{
  for(size_t i=0;i<m.size();i++)
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
  for (size_t i=1 ; i<=IN.size();i++)
  {
    vector<size_t> pointers;
    for(size_t k=0;k<i;k++)
    {
      pointers.push_back(k);
    }
    for (size_t j=0 ; j<kofn(i,IN.size());j++)
    {
      ltAssignments cvec;
      for(size_t idx=0;idx<i;idx++)
      {
        if(!assignmentContainsNewOld(cvec,
                                  IN[pointers[idx]].nlrv,
                                  IN[pointers[idx]].plrv))
          cvec.push_back(IN[pointers[idx]]);
      }
      if(cvec.size()==i)
        OUT.push_back(cvec);
      for(size_t inc=i;inc>0;inc--)
      {
        if(pointers[inc-1]<IN.size()-1-(i-inc))
        {
          pointers[inc-1]++;
          size_t add=0;
          for(size_t res=inc;res<i;res++)
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
                   map<size_t,size_t> &newAssignments)
{
  for(size_t i=0;i<m.size();i++)
  {
    newAssignments[m[i].nlrv]=m[i].plrv;
  }
}

//Return the accuracy of a set of mappings
double mappingAccuracy(ltAssignments &a)
{
  double SUM=0;
  for(size_t i=0;i<a.size();i++)
  {
    SUM+=a[i].getDistance();
  }
  return SUM/(1*(1+(a.size()-1)*1.5));
}

//Generate all the available reasonable mappings
void createReasonableMappings(
  vector<size_t> &candidateLarvae,
  vector<size_t> &newLarvae,
  double duration,
  vector<ltAssignments> &accepted_mappings)
{
    ltAssignments initial_pairs;
    vector<lrvMapping> valid_pairs;
    map<size_t,size_t> mappable_new;
    for(size_t i=0;i<candidateLarvae.size();i++)
    {
      for(size_t j=0;j<newLarvae.size();j++)
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
  vector<size_t> &candidateLarvae,
  vector<size_t> &newLarvae,
  map<size_t, size_t> &newAssignments,
  double duration
  )
{
  vector<ltAssignments> valid_mappings;
  createReasonableMappings(
    candidateLarvae,
    newLarvae,
    duration,
    valid_mappings);

  map<double,size_t> mapping_accuracy;
  for(size_t i=0;i<valid_mappings.size();i++)
  {
    double acc=mappingAccuracy(valid_mappings[i]);
    //cerr << "=========================" << endl;
    //cerr << "Mapping " << i << ": " << endl;
    //for(size_t j=0;j<valid_mappings[i].size();j++)
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
                      size_t CLUSTER_ID,
                      vector<size_t> &IDs
                      )
{
  // We have the following cases here:
  //  1) CLUSTER_ID Corresponds to no cluster: This means
  //     the cluster is newly found (it started as a cluster)
  //  2) CLUSTER_ID Corresponds to an existing cluster:
  //     we assign all the new blobs as decendants of the previous blob
  map<size_t,vector<size_t> >::iterator dcIT;
  //Not found new cluster NEW IDs to be given to the vector
  //elements
  BOOST_LOG_TRIVIAL(trace) << "Cluster " << CLUSTER_ID 
                  << " diverged. Assigning new IDs for diverged larvae";
  assign_one(CLUSTER_ID,IDs);
  //detected_larvae[CLUSTER_ID].isCluster=true;
}

int detect_diverging(vector<size_t> &preLarvaeNearby,
                       vector<size_t> &newLarvaeNearby,
                       cvb::CvBlobs &Pre,
                       cvb::CvBlobs &New)
{
  BOOST_LOG_TRIVIAL(trace) << "Trying to detect diverging clusters";

  BOOST_LOG_TRIVIAL(trace)<< "Size of newLarvaeNearby: "  
                          << newLarvaeNearby.size();

  if(newLarvaeNearby.size()<=1)
    return -1; // No diverging is possible
  // Check the complete set first
  vector<size_t>::iterator pIT=preLarvaeNearby.begin();
  while(pIT!=preLarvaeNearby.end())
  {
    //BOOST_LOG_TRIVIAL(trace) << "Checking if nodes " 
    //                         << printVector(newLarvaeNearby) 
    //                         << " diverged from: " << *pIT;
    double v;
    if(centresMatch(New,Pre[*pIT],newLarvaeNearby,v))
    {
      // Centres of all candidates match with new blob
      // cluster contains all. We can return
      //BOOST_LOG_TRIVIAL(trace) << "Node " << *pIT << " matches " 
      //                         << printVector(newLarvaeNearby);
      assign_divergingOL(New,*pIT,newLarvaeNearby);
      break;
    }
    if(newLarvaeNearby.size()>2)
    {
      //BOOST_LOG_TRIVIAL(trace) << "Checking powersets of " 
      //                         << printVector(newLarvaeNearby) 
      //                         << " that diverged from: " << *pIT;
      vector<vector<size_t> > pSETS;
      powersets(newLarvaeNearby,pSETS);
      vector<vector<size_t> >::iterator pSIT=pSETS.begin();
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

int detect_clustering(vector<size_t> &preLarvaeNearby,
                       vector<size_t> &newLarvaeNearby,
                       cvb::CvBlobs &Pre,
                       cvb::CvBlobs &New)
{
  BOOST_LOG_TRIVIAL(trace) << "Trying to detect clusters";

  if(preLarvaeNearby.size()<=1)
    return -1; // No clustering is possible
  
  // Check the complete set first
  vector<size_t>::iterator nIT=newLarvaeNearby.begin();
  while(nIT!=newLarvaeNearby.end())
  {
    double v;
    if(centresMatch(Pre,New[*nIT],preLarvaeNearby,v) &&
       blobSizeIsRelevant(Pre,New[*nIT],preLarvaeNearby))
    {
      // Centres of all candidates match with new blob
      // cluster contains all. We can return
      //BOOST_LOG_TRIVIAL(trace) << "Centres of new larva: "
      //                         << *nIT << " match centres of " 
      //                         << printVector(preLarvaeNearby) 
      //                         << ". Assigning clustering";

      assign_clustering(*nIT,preLarvaeNearby);
      /* LarvaFit Test starts here: */
      
      /* LarvaFit Test ends here: */

      //continue;
      break;
    }
    vector<vector<size_t> > pSETS;
    powersets(preLarvaeNearby,pSETS);
    vector<vector<size_t> >::iterator pSIT=pSETS.begin();
    while(pSIT!=pSETS.end())
    {
      BOOST_LOG_TRIVIAL(trace) 
        << "Trying to detect subsets for clustering";
      if(centresMatch(Pre,New[*nIT],*pSIT,v) &&
          blobSizeIsRelevant(Pre,New[*nIT],*pSIT))
      {
        // Centres of all candidates match with new blob
        // cluster contains subset pSIT. We can return
        //BOOST_LOG_TRIVIAL(trace) << "Centres of new larva: " 
        //                         << *nIT << " match centres of subset " 
        //                         << printVector(*pSIT) 
        //                         << ". Assigning clustering"
        //                         << endl;
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
  size_t i=0;
  cvb::CvBlob blob=lrv.blobs.back();
  vector<Point2f> &spine=lrv.lrvDistances.back().Spine;
  Point2f bp(lrv.blobs.back().minx,lrv.blobs.back().miny);
  Point2f sp_front=spine[0]-bp;
  Point2f sp_back=spine.back()-bp;
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
}


void updateOneLarva(cvb::CvBlobs &In,
                    cvb::CvBlobs &Prev,
                    cvb::CvBlobs::iterator it,
                    tbb::concurrent_hash_map<size_t, larvaObject> &NEW_LARVA)
{
  size_t ID=(*it).first;
  cvb::CvBlob blob=*((*it).second);
  Mat larvaROI,cntPoints;
  createLarvaContour(larvaROI,blob);
  createLarvaContourPoints(cntPoints,blob);

  map<size_t,larvaObject>::iterator curLarva;
  // NEW LARVA OBJECT!
  if ((curLarva=detected_larvae.find(ID))==detected_larvae.end()
      ||
      Prev==In)
  {
    // Create and allocate the new object
    larvaObject newLarva;

    // Set the frame of it's existence
    newLarva.start_frame=CURRENT_FRAME;
    newLarva.end_frame=CURRENT_FRAME;

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
    //cerr << newLarva.larva_ID << ": " << newLarva.roundness.back() << endl;
    // In this block we compute the inner spine for the larva
    vector<Point2f> cntPoints;
    blobToPointVector(blob,cntPoints);
    larvaDistanceMap Distances(cntPoints);
    //computeSpine(blob,Distances,frame);
    fixContour(blob,Distances,
        LRVTRACK_CONTOUR_RESOLUTION,
        colorFrame,previousFrame);
    newLarva.lrvDistances.push_back(Distances);
    newLarva.length.push_back(Distances.MaxDist);
    newLarva.length_mean = Distances.MaxDist;
    newLarva.length_sum= Distances.MaxDist;
    newLarva.length_max = Distances.MaxDist;
    newLarva.length_min = Distances.MaxDist;
    Point2f Head,Tail;
    Point2f bp(newLarva.blobs.back().minx,newLarva.blobs.back().miny);
    Head=Distances.Spine[0]-bp;
    Tail=Distances.Spine.back()-bp;
    newLarva.heads_brightness.push_back(getSurroundingSize(Head,blob,greyFrame));
    newLarva.tails_brightness.push_back(getSurroundingSize(Tail,blob,greyFrame));
    newLarva.angular_speed.push_back(0);

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

    newLarva.round_flag.push_back(false);


    /*if(LRVTRACK_CSV_OUTPUT)
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
    }*/
    newLarva.lastFrameWithStats=CURRENT_FRAME;

    //detected_larvae[ID]=newLarva;
    //NEW[ID]=newLarva;
    tbb::concurrent_hash_map<size_t,larvaObject>::accessor a;
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
    cur_larva.end_frame++;

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


    // If not then:
    //  Update area values for larva.

    ++cur_larva.lifetimeWithStats;
    cur_larva.lastBlobWithStats=cur_larva.blobs.size()-1;
    //cur_larva.isCluster=false;

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

    cur_larva.centroid_distance_x.push_back(
        fabs(blob.centroid.x - preBlob.centroid.x));
    cur_larva.centroid_distance_y.push_back(
        fabs(blob.centroid.y - preBlob.centroid.y));
    cur_larva.centroid_distance_x_sum+=
      fabs(blob.centroid.x - preBlob.centroid.x);
    cur_larva.centroid_distance_x_sum+=
      fabs(blob.centroid.y - preBlob.centroid.y);


    // Point coordinates for head/tail
    Point2f Head,Tail;

    // Map to keep the distances of each point to all the others
    // Pair of points to keep the points with the Maximum distance 
    // (i.e. head and tail :) )
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

    double prespDist=p2fdist(cur_larva.heads.back(),cur_larva.tails.back());
    //execute the function to find which is which and assign appropriately
    //to Head/Tail.
    Head=Distances.Spine[0];
    Tail=Distances.Spine.back();
    findHeadTail(cur_larva,Head,Tail);
    cur_larva.heads.push_back(Head);
    cur_larva.tails.push_back(Tail);
    cur_larva.heads_brightness.push_back(getSurroundingSize(Head,blob,greyFrame));
    cur_larva.tails_brightness.push_back(getSurroundingSize(Tail,blob,greyFrame));

    double preRnd=cur_larva.roundness.back();
    double preLen=cur_larva.length[cur_larva.length.size()-2];
    cur_larva.roundness.push_back((perimeter*perimeter)/(2*CV_PI*blob.area));
    double curRnd=cur_larva.roundness.back();
    double curLen=cur_larva.length.back();
    double RndRat=curRnd/preRnd;
    double spDist=p2fdist(cur_larva.heads.back(),cur_larva.tails.back());
    double var=Distances.curvatureVariance;
    //double spDist=distanceBetweenSpines(
    //             Distances,
    //             cur_larva.lrvDistances[cur_larva.lrvDistances.size()-2]);

    //Add head and tail to the history
    if(check_roundness(cur_larva.larva_ID,spDist,curRnd,RndRat,prespDist,var))
    {
      cur_larva.round_flag.push_back(true);
    }
    else if(cur_larva.round_flag.back()==true && 
           !check_roundness(cur_larva.larva_ID,spDist,curRnd,RndRat,prespDist,var))
    {
      /*BOOST_LOG_TRIVIAL(debug) << "RND, E, " << CURRENT_FRAME 
        << ", " << cur_larva.larva_ID
        << ", " << curRnd
        << ", " << RndRat
        << ", " << spDist;*/
      cur_larva.round_flag.push_back(false);
    }
    else if(cur_larva.round_flag.back()==true)
    {
      /*BOOST_LOG_TRIVIAL(debug) << "RND, -, " << CURRENT_FRAME 
        << ", " << cur_larva.larva_ID
        << ", " << curRnd
        << ", " << RndRat
        << ", " << spDist;*/
      cur_larva.round_flag.push_back(true);
    }
    else if(cur_larva.round_flag.back()==false)
    {
      /*BOOST_LOG_TRIVIAL(debug) << "RND, -, " << CURRENT_FRAME 
        << ", " << cur_larva.larva_ID
        << ", " << curRnd
        << ", " << RndRat
        << ", " << spDist;*/
      cur_larva.round_flag.push_back(false);
    }
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

    /*if(LRVTRACK_CSV_OUTPUT)
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
    }*/

    cur_larva.lastFrameWithStats=CURRENT_FRAME;
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
  cvReleaseImage(&labelImg);
}

/* The function removes the background from the new frame
   and applies thresholding to derive the image in which we look
   for blobs.
    Input:
       * newFrame: the new frame from the get_next_frame function
       * greyBgFrame: a grey image of the background
       * showFrame: the BGR image returned to the user
    Output:
       * processedFrame: the frame containing only the greyscale shapes
                         of the larvae
*/
void process_frame(Mat &newFrame,
                  Mat &greyBgFrame,
                  Mat &showFrame,
                  Mat &processedFrame)
{
    int THRESH=THRESH_BINARY;
    size_t CMAX=255;
    size_t NB=137;
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
    //addWeighted(newFrame, 1.0, greyBgFrame, -3.0, 0.0, fgFrame);

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
      for(size_t i=0; i<cups.size(); ++i)
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
    //equalizeHist(processedFrame, processedFrame);
}

//Match lost larvae
size_t matchLostLarvae(size_t newLabel)
{
  for(vector<size_t>::iterator lIT=lost_larvae.begin();
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
        size_t ret=*lIT;
        lost_larvae.erase(lIT);
        return ret;
      }
  }
  return 0;
}

void updateLarvae(cvb::CvBlobs &In, cvb::CvBlobs &Prev)
{
  //cvb::CvBlobs::iterator it=In.begin();
  tbb::concurrent_hash_map<size_t, larvaObject> NEW_LARVA;
  larvaeUpdateBody body(In,Prev,NEW_LARVA);
  //parallel_for_(Range(0, In.size()), body);
  parallel_for(BlockedRange(0, In.size()), body);
  
  tbb::concurrent_hash_map<size_t, larvaObject>::iterator nit=
    NEW_LARVA.begin();
  while (nit!=NEW_LARVA.end())
  {
    detected_larvae[nit->first]=nit->second;
    ++nit;
    }
}



/*
 * TODO:
 * IMPORTANT NOTE!!
 *  We have essentially three mappings:
 *   1) The previous frame blobs -> new frame blobs
 *   2) The new frame blobs -> previous frame blobs
 *   3) The new frame blobs -> updated numbers
 * 
 *  The updated numbers are:
 *   a) the numbers that remained from the previous frame.
 *   b) the numbers that diverged from frames before the previous
 *   c) new numbers (increasing LARVAE_COUNT)
 */
void newLarvaeTrack(cvb::CvBlobs &In, 
                    cvb::CvBlobs &Prev, 
                    cvb::CvBlobs &out)
{


  assignedNew.clear();
  assignedPrevious.clear();
  assignedPreMap.clear();
  newInFrame.clear();
  //newClusters.clear();

  cvb::CvBlobs::iterator prevIt=Prev.begin();
  while (prevIt!=Prev.end())
  {
    assignedPreMap[prevIt->first]=0;
    ++prevIt;
  }

  BOOST_LOG_TRIVIAL(trace) << "========  Starting Tracking loop ==========";

	prevIt=Prev.begin();
  while (prevIt!=Prev.end())
  {
    //If the larvae from the previous frame has not been assigned
    //  -- It may be that it is assigned during a previous larva
    //     inspection.
    if(assignedPreMap[prevIt->first]==0)
    {
      cvb::CvBlob *preBlob=prevIt->second;
      size_t preID=prevIt->first;

      // Get larva around it before and after
      vector<size_t> preLarvaeNearby;
      vector<size_t> postLarvaeNearby;
      getNearbyLarvae(Prev,preBlob,preLarvaeNearby,true);
      getNearbyLarvae(In,preBlob,postLarvaeNearby,false);

      //BOOST_LOG_TRIVIAL(trace) << "Around larva " << prevIt->first 
      //                         << " we found " << preLarvaeNearby.size() 
      //                         << " larvae from the previous frame and " 
      //                         << postLarvaeNearby.size() 
      //                         << " from the new frame.";

      //BOOST_LOG_TRIVIAL(trace) << "Around larva " << prevIt->first 
      //                         << " we found pre: " 
      //                         << printVector(preLarvaeNearby) 
      //                         << " and post: " 
      //                         << printVector(postLarvaeNearby);
      // Great case collection now:
      
      // CASE 1 -> 1
      if(preLarvaeNearby.size()==1 && postLarvaeNearby.size()==1)
      {
        if(blobSizeIsRelevant(In[postLarvaeNearby[0]],preBlob))
        {
          double v;
          if(centresMatch(In,Prev[preLarvaeNearby[0]],postLarvaeNearby,v,0.3))
          {
            //BOOST_LOG_TRIVIAL(trace) << "Larva " << prevIt->first 
            //                         << " from previous frame matches " 
            //                         << postLarvaeNearby[0] 
            //                         << " from new frame. Assigning.";
            assign_one(preID,postLarvaeNearby[0]);
          }
          else
          {
            // Centres do NOT match! Maybe it disappeared.
            // TODO: Check
            BOOST_LOG_TRIVIAL(trace) 
              << "1-1: with size similar and centres not matching. !UNHANDLED!";
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
          BOOST_LOG_TRIVIAL(trace) 
            << "1-1: with size too different. No assignment.";
        }
      }
      // END OF 1-1 CASE
      // GENERAL CASE:
      else
      {
        // Try to assign obvious ones:
        vector<size_t>::iterator preNearIt=preLarvaeNearby.begin();
        while(preNearIt!=preLarvaeNearby.end())
        {
          bool erased=false;
          vector<size_t>::iterator postNearIt=postLarvaeNearby.begin();
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
              //BOOST_LOG_TRIVIAL(trace) << "Larva " << *preNearIt 
              //                         << " from previous frame matches " 
              //                         << *postNearIt 
              //                         << " from new frame. Assigning";
              assign_one(*preNearIt,*postNearIt);
              // erase the obvious assignments
              try{
                preLarvaeNearby.erase(preNearIt);
                postLarvaeNearby.erase(postNearIt);
              }
              catch(...)
              {
                BOOST_LOG_TRIVIAL(warning) << "Problem erasing: " 
                  << *preNearIt << " or " << *postNearIt;
              }
              erased=true;
            }
            else{
              //BOOST_LOG_TRIVIAL(trace) << "Larva " << *preNearIt 
              //                         <<" from previous frame doesn't match "
              //                         << *postNearIt << " from new frame.";
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
        //BOOST_LOG_TRIVIAL(trace)<< "After initial assignment assignedPreMap: "
        //                        << printUIMap(assignedPreMap) 
        //                        << " for larva with Id: " 
        //                        << preID;
        if(assignedPreMap[preID]>0)
        {
          continue;
        }

        detect_clustering(preLarvaeNearby,postLarvaeNearby,Prev, In);

        //BOOST_LOG_TRIVIAL(trace) << "After clustering check assignedPreMap: "
        //                         << printUIMap(assignedPreMap) 
        //                         << " for larva with Id: " << preID;
        if(assignedPreMap[preID]>0)
        {
          continue;
        }

        detect_diverging(preLarvaeNearby,postLarvaeNearby,Prev, In);
        //BOOST_LOG_TRIVIAL(trace) << "After diverging check assignedPreMap: " 
        //                         << printUIMap(assignedPreMap) 
        //                         << " for larva with Id: " << preID;
        if(assignedPreMap[preID]>0)
        {
          continue;
        }
        else
          BOOST_LOG_TRIVIAL(trace) 
            << "FOUND TESTCASE FOR MIXED DIVERGING/CONVERGING";
      }
    }
    ++prevIt;
    map<size_t, vector<size_t> >::iterator prI=assignedPrevious.begin();
    while(prI!=assignedPrevious.end())
    {
      if(prI->second.size()>0)
      {
      }
      else
      {
        lost_larvae.push_back(prI->first);
      }
      prI++;
    }

    /*prI=assignedNew.begin();
    while(prI!=assignedNew.end())
    {
      if(prI->second.size()>0)
      {
        //BOOST_LOG_TRIVIAL(trace) << prI->first << " -> " 
        //                         << printVector(prI->second);
      }
      else
      {
        //BOOST_LOG_TRIVIAL(trace) << prI->first << " -> []" ;
      }
      prI++;
    }*/
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
      //size_t Match=matchLostLarvae(bIT->first);
      //if(Match==0)
      //{
        size_t NEWID=++LARVAE_COUNT;
        newInFrame.push_back(NEWID);
        bIT->second->n20=bIT->second->label;
        bIT->second->label=NEWID;
        out[NEWID]=bIT->second;
      /*}
      else
      {
        bIT->second->n20=bIT->second->label;
        bIT->second->label=Match;
        out[Match]=bIT->second;
      }*/
    }
    bIT++;
  }

  updateLarvae(out,Prev);
}

/* Function to extract and process each frame so that the background is dark
 * and the forground white.
 *   Input: 
 *       * capture: the capture device
 *   Output:
 *       * output: the processed frame
 *       * origFrame: the RGB version of the image
 */
bool get_next_frame(VideoCapture &capture, Mat &output, Mat &colorOutput)
{
  CURRENT_FRAME++;
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
    size_t PF=0;
    int THRESHOLD=255;
    double cpow=1.25;
    ocl_output=ocl_output*MF+PF;
    ocl::addWeighted(ocl_output,1.0,ocl_output,-2.0,THRESHOLD,ocl_output);
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
    size_t PF=0;
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

void extract_background_offline(VideoCapture &capture,
                       Mat &greyBgFrame)
{
  Mat tmpFrame;
  Mat origFrame;
  Mat tmpFrame32f;
  size_t width=capture.get(CV_CAP_PROP_FRAME_WIDTH);
  size_t height=capture.get(CV_CAP_PROP_FRAME_HEIGHT);
  Mat resultframe = Mat::zeros(height,width,CV_32FC1);

  size_t count=capture.get(CV_CAP_PROP_FRAME_COUNT);
  while(get_next_frame(capture,tmpFrame,origFrame)) 
  {
    tmpFrame.convertTo(tmpFrame32f,CV_32FC1);
    resultframe+=tmpFrame32f; 
  }
  capture.set(CV_CAP_PROP_POS_FRAMES ,0);
  CURRENT_FRAME=0;
  resultframe*=(1.0/count);
  normalize(resultframe,greyBgFrame,0,255,CV_MINMAX);
  greyBgFrame.convertTo(greyBgFrame,CV_8UC1);
  imshow("background",greyBgFrame);
  waitKey(10000);
}

/*
 * Function to extract the background
 * The function retrieves a frame from the capture device and 
 *  returns the background into the bgFrame image.
 *    Input:
 *     * capture: the video capture device
 *    Output:
 *     * greyBgFrame: the Suggested background
 *
 */
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
  size_t CMAX=255;
  size_t NB=275;
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
  size_t mysz=circles.size();
  for(size_t i=0; i<circles.size();i++)
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

/*
 * Function to take care of the various input possibilities. 
 * Set up the parameters in the capture device for direct camera 
 * input/file input and correct FPS.
 */
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
  TOTAL_FRAMES=capture.get(CV_CAP_PROP_FRAME_COUNT);
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
      ("dstep,d",
       po::value<int>(&LRVTRACK_DSTEP),
       "Dstep.")
      ("wstep,w",
       po::value<double>(&LRVTRACK_WSTEP),
       "Wstep.")
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
          BOOST_LOG_TRIVIAL(error) << "Camera listing not yet implemented." ;
          exit(1);
        }

      if (vm.count("file-input")<=0 && vm.count("camera-input")<=0)
        {
          BOOST_LOG_TRIVIAL(error) << "Error: No Input Specified. Exiting..." ;
          exit(2);
        }
      else if (vm.count("file-input")>0 && vm.count("camera-input")>0)
        {
          BOOST_LOG_TRIVIAL(error) 
            << "Error: Ambiguous Input Specified. Exiting...";
          exit(2);
        }

      if (vm.count("file-input")>0)
        {
          LRVTRACK_FILE_INPUT=vm["file-input"].as<string>();
          if (LRVTRACK_FILE_INPUT=="")
            {
              BOOST_LOG_TRIVIAL(error) 
              << "Error: Input file flag given but no file specified. Exiting...";
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
      BOOST_LOG_TRIVIAL(error) << "Problem parsing options: " << e.what();
      exit(1);
    }

  return 0;
}

/* Second Pass 
 *
 * The functions below belong to the second pass of data
 *  - Verify larvae and non-larvae objects
 *  - Detect heads/tails
 *
*/

size_t evalPartitionSize(std::vector<std::vector<size_t> > &v,
                          std::vector<size_t> &blobs)
{
  size_t vsz=v.size();
  size_t bsz=blobs.size();
  size_t totalError=0;
  size_t b=0;
  for(auto &i:v)
  {
    size_t expected_size=0;
    for(auto &l:i)
      expected_size+=detected_larvae[l].area.back();
    totalError+=abs(detected_larvae[blobs[b++]].area.back()-
                    expected_size);
  }
  return totalError;
}

void findEndObjectsFwHelper(larvaObject &t,
                            size_t c, 
                            size_t objnum,
                            std::map<size_t,size_t> &blobSizeMap,
                            std::vector<size_t> &e)
{
  if(e.size()==objnum)
    return;
  BOOST_LOG_TRIVIAL(debug) << "Verify: FindFW: Main Obj: " <<
      t.larva_ID << " CurObj: " << c << " Objnum: " << objnum;
  if(blobSizeMap[c]==1)
  {
    BOOST_LOG_TRIVIAL(debug) << 
      "Verify: FindFW: blobSize 1 and not Cluster so we push: " << c;
    e.push_back(c);
  }
  else
  {
    if(!children_blobs[c].empty())
    {
      if(e.size()==objnum)
        return;
      BOOST_LOG_TRIVIAL(debug) << 
        "Verify: FindFW: children not empty: " << c;
      BOOST_LOG_TRIVIAL(debug) << "Verify: FindFW: " 
        << printVector(children_blobs[c]);
      for(auto &n: children_blobs[c])
      {
        if(e.size()==c)
          break;
        findEndObjectsFwHelper(t,n,objnum,blobSizeMap,e);
      }
    }
  }
}

void findEndObjectsBwHelper(larvaObject &t,
                            size_t c, 
                            size_t objnum,
                            std::map<size_t,size_t> &blobSizeMap,
                            std::vector<size_t> &e)
{
  if(e.size()==objnum)
    return;
  BOOST_LOG_TRIVIAL(debug) << "Verify: FindBW: Main Obj: " <<
      t.larva_ID << " CurObj: " << c << " Objnum: " << objnum;
  if(blobSizeMap[c]==1)
  {
    BOOST_LOG_TRIVIAL(debug) << 
      "Verify: FindBW: blobSize 1 and not Cluster so we push: " << c;
    e.push_back(c);
  }
  else
  {
    if(!parent_blobs[c].empty())
    {
      BOOST_LOG_TRIVIAL(debug) << 
        "Verify: FindBW: parents not empty: " << c;
      BOOST_LOG_TRIVIAL(debug) << "Verify: FindBW: " 
        << printVector(parent_blobs[c]);
      for(auto &n: parent_blobs[c])
      {
        findEndObjectsBwHelper(t,n,objnum,blobSizeMap,e);
      }
    }
  }
}

void findEndObjects(larvaObject &t, 
                    std::vector<size_t> &e,
                    std::map<size_t,size_t> &blobSizeMap,
                    size_t obj_num,
                    bool forward)
{
  if(forward)
    findEndObjectsFwHelper(t,t.larva_ID,obj_num,blobSizeMap,e);
  else
    findEndObjectsBwHelper(t,t.larva_ID,obj_num,blobSizeMap,e);

  if(e.size()!=t.blobSize)
  {
    BOOST_LOG_TRIVIAL(debug) << "WARNING: Cannot find end Objects "
      << "Larva: " << t.larva_ID << " with size " << t.blobSize;
  }
}

void clarifyLarvae(larvaObject &l,
                   std::vector<size_t> &neighbors, 
                   size_t neighborsize,
                   std::vector<size_t> &allfound,
                   std::map<size_t,size_t> &blobSizeMap,
                   bool &change)
{
    if(neighborsize>allfound.size())
    {
      BOOST_LOG_TRIVIAL(debug) << "Clarify: For Larva: " <<
        l.larva_ID << " WARN Nsize=" <<
        neighbors.size() << ", ALL=" << allfound.size();
      return;
    }
    partition_generator<size_t> pgs(allfound,neighbors.size());
    pgs.m_partitions_of_n();
    vector<vector <size_t> > &BestPartition=pgs.RES[0];
    size_t minError=SIZE_MAX;
    for(auto &pg: pgs.RES)
    {
      size_t ne;
      ne=evalPartitionSize(pg,neighbors);
      if(ne<minError)
      {
        minError=ne;
        BestPartition=pg;
        change=1;
      }
    }
    size_t newsz=0;
    for(auto i=0;i<neighbors.size();i++)
    {
      size_t sz=0;
      for(auto &s:BestPartition[i])
        sz+=detected_larvae[s].blobSize;

      newsz+=sz;
      detected_larvae[neighbors[i]].blobSize=sz;
      blobSizeMap[neighbors[i]]=sz;
      detected_larvae[neighbors[i]].isCluster=(sz-1);
      BOOST_LOG_TRIVIAL(debug) << "Verify: FIX: OBJECT" 
        << neighbors[i]
        << " updated with blobSize: " << sz;
      change=true;
    }
    //l.blobSize=newsz;
    //if(newsz>1)
    //  l.isCluster=true;
    //BOOST_LOG_TRIVIAL(debug) << "Verify: FIX: OBJECT" 
    //  << l.larva_ID
    //  << " updated with blobSize: " << newsz;
}

/*void verifyLarva(larvaObject &f,bool &change)
{
  if(!f.diverged_to.empty())
  {
    size_t sum=0;
    for(auto &v: f.diverged_to)
      sum+=detected_larvae[v].blobSize;
    BOOST_LOG_TRIVIAL(debug) << "Verify: Checking: CLUSTER" 
      << f.larva_ID
      << " has blobSize: " << f.blobSize
      << " and the sum of those diverging from it is: " << sum;
    if(f.blobSize!=sum)
    {
      BOOST_LOG_TRIVIAL(debug) << "Verify: WARN: CLUSTER" 
        << f.larva_ID
        << " has blobSize: " << f.blobSize
        << " and the sum of those diverging from it is: " << sum;
      BOOST_LOG_TRIVIAL(debug) << "Verify: WARN: From: " 
        << printVector(f.diverged_to);
      std::vector<size_t> endObjs;
      f.blobSize=std::max(sum,(size_t) f.blobSize);

      if(f.blobSize>1)
        f.isCluster=true;
      findEndObjects(f,endObjs,1);
      if(endObjs.size()!=f.blobSize)
      {
        endObjs.clear();
        findEndObjects(f,endObjs,0);
      }
      if(endObjs.size()!=f.blobSize)
      {
        BOOST_LOG_TRIVIAL(debug) << "Verify: WARN: for " << f.larva_ID
          << " we cannot find enough end objects." ;
        return;
      }
        
      BOOST_LOG_TRIVIAL(debug) << "Verify: Clarifying: " << f.larva_ID 
        << " with end objects: "<< printVector(endObjs)
        << " for ObjectSize: " << f.blobSize;
      clarifyLarvae(f,f.diverged_to,sum,endObjs,change);
    }
  }
  if(parent_blobs[f.larva_ID].size()>1)
  {
    size_t sum=0;
    for(auto &v: parent_blobs[f.larva_ID])
      sum+=detected_larvae[v].blobSize;
    BOOST_LOG_TRIVIAL(debug) << "Verify: Checking: CLUSTER" 
      << f.larva_ID
      << " has blobSize: " << f.blobSize
      << " and the sum of those diverging from it is: " << sum;
    if(f.blobSize!=sum)
    {
      BOOST_LOG_TRIVIAL(debug) << "Verify: WARN: CLUSTER" 
        << f.larva_ID
        << " has blobSize: " << f.blobSize
        << " and the sum of those diverging from it is: " << sum;
      BOOST_LOG_TRIVIAL(debug) << "Verify: WARN: From: " 
        << printVector(parent_blobs[f.larva_ID]);
      std::vector<size_t> endObjs;
      f.blobSize=std::max(sum,(size_t) f.blobSize);

      if(f.blobSize>1)
        f.isCluster=true;
      findEndObjects(f,endObjs,0);
      if(endObjs.size()!=f.blobSize)
      {
        endObjs.clear();
        findEndObjects(f,endObjs,0);
      }
      if(endObjs.size()!=f.blobSize)
      {
        BOOST_LOG_TRIVIAL(debug) << "Verify: WARN: for " << f.larva_ID
          << " we cannot find enough end objects." ;
        return;
      }
        
      BOOST_LOG_TRIVIAL(debug) << "Verify: Clarifying: " << f.larva_ID 
        << " with end objects: "<< printVector(endObjs)
        << " for ObjectSize: " << f.blobSize;
      clarifyLarvae(f,parent_blobs[f.larva_ID],sum,endObjs,change);
    }
  }
}*/

void determineHeadTail(larvaObject &f)
{
  stringstream d;
  try{
  d << endl;
  d << "|=============== Determine head/tail of: " << f.larva_ID 
    << " ============= " << endl;
  if(f.isCluster==true)
  {
    d << "|=============== Blob " << f.larva_ID << " is a cluster ========"
      << endl ;
    d << "|========================================= " << endl;
    BOOST_LOG_TRIVIAL(debug) << d.str();
    return;
  }
  size_t range=f.blobs.size();
  std::vector<size_t> breaks;
  for(size_t i=0;i<range;i++)
  {
    if(f.round_flag[i]==true)
      breaks.push_back(i);
    while(f.round_flag[i]==true)
      i++;
  }

  breaks.push_back(range-1);
   d << "|  Range: " << range << endl;
   d << "|  Breaks: " << breaks.size() << " " << printVector(breaks) << endl;

  for(size_t i=0;i<breaks.size();i++)
  {
    d << "|    Break " << i << ":" << endl;
    size_t endBreak=breaks[i];
    size_t preBreak;
    if(i>0)
      preBreak=breaks[i-1]+1;
    else
      preBreak=0;

    Point2f startPos(f.blobs[preBreak].centroid.x,
                     f.blobs[preBreak].centroid.y);
    Point2f endPos(f.blobs[endBreak].centroid.x,
                   f.blobs[endBreak].centroid.y);

    Point2f startbp(f.blobs[preBreak].minx,f.blobs[preBreak].miny);
    Point2f endbp(f.blobs[endBreak].minx,f.blobs[endBreak].miny);

    d << "|      Start[ " << preBreak <<"]: Centroid" << startPos << " "
      << "IHead" << f.heads[preBreak]+startbp << " "
      << "ITail" << f.tails[preBreak]+startbp 
      << endl;
    d << "|      End[ " << endBreak <<"]: Centroid" << endPos << " "
      << "IHead" << f.heads[endBreak]+endbp << " "
      << "ITail" << f.tails[endBreak]+endbp
      << endl;

    double sumHeadDist=0.0; //Distance to endpoint
    double sumTailDist=0.0;
    double sumHeadGrey=0.0; //Grey value sum
    double sumTailGrey=0.0;
    double sumHeadDiff=0.0; //Difference of current centroid 
                            //from previous endpoints
    double sumTailDiff=0.0;
    double sumHeadCnt=0.0;  //Distance covered by endpoint
    double sumTailCnt=0.0;

    for(size_t j=preBreak;j<=endBreak;j++)
    {
      Point2f bpPoint(f.blobs[j].minx,f.blobs[j].miny);
      Point2f hPoint=f.heads[j]+bpPoint;
      Point2f tPoint=f.tails[j]+bpPoint;
      sumHeadDist+=p2fdist(hPoint,endPos);
      sumTailDist+=p2fdist(tPoint,endPos);
      sumHeadGrey+=f.heads_brightness[j];
      sumTailGrey+=f.tails_brightness[j];
      if(j>preBreak)
      {
        Point2f prebpPoint(f.blobs[j-1].minx,f.blobs[j-1].miny);
        Point2f prehPoint=f.heads[j-1]+prebpPoint;
        Point2f pretPoint=f.tails[j-1]+prebpPoint;
        Point2f centroid(f.blobs[j].centroid.x,f.blobs[j].centroid.x);
        Point2f midpoint=f.lrvDistances[j].MidPoint;
        sumHeadDiff+=p2fdist(midpoint,prehPoint);
        sumTailDiff+=p2fdist(midpoint,pretPoint);
        sumHeadCnt+=p2fdist(hPoint,prehPoint);
        sumTailCnt+=p2fdist(tPoint,pretPoint);
      }
    }
    double vote_for_nochange=0;

    if(sumTailDiff>sumHeadDiff)
      vote_for_nochange+=2;
    else
      vote_for_nochange-=2;

    if(sumTailGrey>sumHeadGrey)
      vote_for_nochange++;
    else
      vote_for_nochange--;

    if(sumHeadCnt>sumTailCnt)
      vote_for_nochange++;
    else
      vote_for_nochange--;

    /*vote_for_nochange+=2.0*(sumTailDiff-sumHeadDiff)/
                            std::max(sumHeadDiff,sumTailDiff);

    vote_for_nochange+=1.0*(sumTailGrey-sumHeadGrey)/
                            std::max(sumHeadGrey,sumTailGrey);

    vote_for_nochange+=1.0*(sumHeadCnt-sumTailCnt)/(
                           std::max(sumHeadCnt,sumTailCnt)*(endBreak-preBreak+1));*/
    size_t duration=endBreak-preBreak;
    d << "|" << endl;
    d << "|      LarvaID, Dhead, GHead, DiffH, CntH, Duration, Votes" << endl;
    d << "|      " <<f.larva_ID*1000 + breaks.size() << ", ";
    d << sumHeadDist/duration << " ,";
    d << sumHeadGrey/duration << ", ";
    d << sumHeadDiff/duration << ", ";
    d << sumHeadCnt/duration << ", " << endBreak-preBreak <<", ";
    d << (vote_for_nochange>=0) << endl;
    d << "|      LarvaID, Dtail, GTail, DiffT, CntT, Duration, Votes" << endl;
    d << "|      " <<f.larva_ID*1000 + breaks.size() << ", ";
    d << sumTailDist/duration << ", ";
    d << sumTailGrey/duration << ", ";
    d << sumTailDiff/duration << ", ";
    d << sumTailCnt/duration << ", " << endBreak-preBreak <<", ";
    d << (vote_for_nochange<0) << endl;
    d << "|" << endl;

    /*d << "|      LarvaID, Dtail, GTail, DiffT, CntT, Duration, Votes" << endl;
    d << "|      " <<f.larva_ID*1000 + breaks.size() << ", ";
    d << sumTailDist - sumHeadDist << ", ";
    d << sumTailGrey - sumHeadGrey << ", ";
    d << sumTailDiff - sumHeadDiff << ", ";
    d << sumTailCnt - sumHeadCnt << ", " << endBreak-preBreak <<", ";
    d << (vote_for_nochange>=0) << endl;
    d << "|" << endl;*/

    if(vote_for_nochange<0 && breaks.size()==1)
    {
      f.heads.swap(f.tails);
    }
    else if(vote_for_nochange<0)
    {
      for(size_t k=preBreak;k<=endBreak;k++)
      {
        cv::Point2f b=f.heads[k];
        f.heads[k]=f.tails[k];
        f.tails[k]=b;
      }
    }
    d << "|      Start[ " << preBreak <<"]: Centroid" << startPos << " "
      << "RHead" << f.heads[preBreak]+startbp << " "
      << "RTail" << f.tails[preBreak]+startbp 
      << endl;
    d << "|      End[ " << endBreak <<"]: Centroid" << endPos << " "
      << "RHead" << f.heads[endBreak]+endbp << " "
      << "RTail" << f.tails[endBreak]+endbp
      << endl;
    d << "|========================================= " << endl;
  }
  BOOST_LOG_TRIVIAL(debug) << d.str();
  }
  catch(...)
  {
    BOOST_LOG_TRIVIAL(debug) << d.str();
    exit(0);
  }
}

#define blobSum(l,g) accumulate(l.begin(), \
                              l.end(),0,\
                              [&g](size_t total, size_t cur)\
                              {return g[cur]+total;} );

void updateBlobSize(larvaObject &l,
                 map<size_t,size_t> &blobSizeMap,
                 bool &change)
{
  size_t lID=l.larva_ID;
  size_t parents=blobSum(parent_blobs[lID],blobSizeMap);
  size_t children=blobSum(children_blobs[lID],blobSizeMap);
  int size_as_child=0;
  int size_as_parent=0;

  if(parents==0 && children==0)
  {
    if(blobSizeMap[lID]==0)
    {
      change=true;
      blobSizeMap[lID]=1;
    }
    return;
  }

  //Node as a child
  if(parents!=0)
  {
    //One parent and being the onlychild
    // Weird case that perhaps doesn't occur.
    if(parent_blobs[lID].size()==1 && 
       children_blobs[parent_blobs[lID].front()].size()==1)
    {
      size_as_child=blobSizeMap[parent_blobs[lID].front()];
    }
    else if(parent_blobs[lID].size()==1 && 
       children_blobs[parent_blobs[lID].front()].size()>1)
    {
      /*size_as_child=blobSizeMap[parent_blobs[lID].front()];
      for(auto &pc: children_blobs[parent_blobs[lID].front()])
      {
        if(pc!=lID)
          size_as_child-=blobSizeMap[pc];
      }*/
      size_as_child=0;
    }
    else if(parent_blobs[lID].size()>1)
    {
      bool onlychild=true;
      for(auto &pp: parent_blobs[lID])
      {
        if(children_blobs[pp].size()>1)
        {
          onlychild=false;
          break;
        }
        if(children_blobs[pp].front()!=lID)
        {
          BOOST_LOG_TRIVIAL(debug) << "Verify: " <<
            "Warning one parent but different child";
          onlychild=false;
          break;
        }
      }
      if(onlychild)
      {
        size_as_child=blobSum(parent_blobs[lID],blobSizeMap);
      }
      else
      {
          BOOST_LOG_TRIVIAL(debug) << "Verify: " <<
            "Warning two parents with foster siblings. Is unhandled";
      }
    }
  }

  //Node as parent
  if(children!=0)
  {
    //One child and being the onlyparent
    // Weird case that perhaps doesn't occur.
    if(children_blobs[lID].size()==1 && 
       parent_blobs[children_blobs[lID].front()].size()==1)
    {
      size_as_parent=blobSizeMap[children_blobs[lID].front()];
    }
    else if(children_blobs[lID].size()==1 && 
       parent_blobs[children_blobs[lID].front()].size()>1)
    {
      /*size_as_parent=blobSizeMap[children_blobs[lID].front()];
      for(auto &pc: parent_blobs[children_blobs[lID].front()])
      {
        if(children_blobs[pc].size()>1)
        {
          size_as_parent=0;//A bit of a mess :)
          break;
        }
        //if(pc!=lID) //Also unhandled
         // size_as_parent-=blobSizeMap[pc];
      }*/
      size_as_parent=0;
    }
    else if(children_blobs[lID].size()>1)
    {
      bool singleparent=true;
      for(auto &pp: children_blobs[lID])
      {
        if(parent_blobs[pp].size()>1)
        {
          singleparent=false;
          break;
        }
        if(parent_blobs[pp].front()!=lID)
        {
          BOOST_LOG_TRIVIAL(debug) << "Verify: " <<
            "Warning one child but different parent";
          singleparent=false;
          break;
        }
      }
      if(singleparent)
      {
        size_as_parent=blobSum(children_blobs[lID],blobSizeMap);
      }
      else
      {
          BOOST_LOG_TRIVIAL(debug) << "Verify: " <<
            "Warning two children with other parents. Is unhandled";
      }
    }
  }

  BOOST_LOG_TRIVIAL(debug) << "Verify: L" 
    << lID
    << " with blobSize: " << blobSizeMap[lID]
    << " parents: " << printVector(parent_blobs[lID])
    << " as parent: " << size_as_parent
    << " children: " << printVector(children_blobs[lID])
    << " as child: " << size_as_child;

  if(size_as_parent>size_as_child && size_as_child>0)
  {
    BOOST_LOG_TRIVIAL(debug) << "Verify: L" 
      << lID
      << " parents>children clarifying";
    vector<size_t> allFwBlobs;
    findEndObjects(l,allFwBlobs,blobSizeMap,size_as_parent,1);
    BOOST_LOG_TRIVIAL(debug) << "Verify: L" 
      << lID
      << " found end Objects:" << printVector(allFwBlobs);
    clarifyLarvae(l,
        parent_blobs[lID],
        parents,
        allFwBlobs,
        blobSizeMap,
        change);
    if(blobSizeMap[lID]<size_as_parent)
    {
      blobSizeMap[lID]=size_as_parent;
      change=true;
    }
  }

  else if(size_as_parent<size_as_child && size_as_parent>0)
  {
    BOOST_LOG_TRIVIAL(debug) << "Verify: L" 
      << lID
      << " parents<children clarifying";
    vector<size_t> allBwBlobs;
    findEndObjects(l,allBwBlobs,blobSizeMap,size_as_child,0);
    BOOST_LOG_TRIVIAL(debug) << "Verify: L" 
      << lID
      << " found end Objects:" << printVector(allBwBlobs);
    clarifyLarvae(l,
        children_blobs[lID],
        children,
        allBwBlobs,
        blobSizeMap,
        change);
    if(blobSizeMap[lID]<size_as_child)
    {
      blobSizeMap[lID]=size_as_child;
      change=true;
    }
  }
  else
  {
    if(max(size_as_parent,size_as_child)>blobSizeMap[lID])
    {
      blobSizeMap[lID]=max(size_as_parent,size_as_child);
      change=true;
    }
  }
}

void secondPass()
{
  /* We loop over each larva object and try to identify head/tail and wether it's a larva or not */

  std::map<size_t,size_t> blobSizeMap;
  std::map<size_t, std::vector<size_t> > detailed_clusters;
  //Initialize
  for (auto &p: detected_larvae)
  {
    size_t lID=p.second.larva_ID;
    blobSizeMap[lID]=max(parent_blobs[lID].size(),
                         children_blobs[lID].size());
    if(blobSizeMap[lID]==0)
      blobSizeMap[lID]=1;

    BOOST_LOG_TRIVIAL(debug) << "Verify: Initialize: L" 
      << lID
      << " with blobSize: " << blobSizeMap[lID];
  }
  size_t iterations=1;
  bool change=true; 
  while(change && iterations<20)
  {
    change=false; 
    BOOST_LOG_TRIVIAL(debug) << "Verify: Entering iteration: " 
                             << iterations;
    change=0;
    BOOST_LOG_TRIVIAL(debug) << "Verify: Forward Scan";
    for (auto p=detected_larvae.rbegin();
              p!=detected_larvae.rend();p++)
    {
      updateBlobSize(p->second,blobSizeMap,change);
    }
    iterations++;
  }

  for (auto &p: detected_larvae)
  {
    if(blobSizeMap[p.second.larva_ID]>1)
    {
      p.second.isCluster=1;
      p.second.blobSize=blobSizeMap[p.second.larva_ID];
    }
  }
  
  std::map<size_t,larvaObject>::iterator dlit;
  dlit=detected_larvae.begin();
  for(;dlit!=detected_larvae.end();++dlit)
  {
    determineHeadTail(dlit->second);
  }
}

int main(int argc, char* argv[])
{
  //waitKey(10000);
  bool SHOWTAGS=true;
  bool TRACK=false;
  bool STEP=true;
  bool STOP=false;
  VideoWriter vidOut,vidPOut;

  log_init();

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
    BOOST_LOG_TRIVIAL(error) 
      << "Error setting up the capture device (camera or video file)";
    exit(1);
  }
  extract_background(capture,bgFrame);

  //Initialize the frame size
  FRAME_COLS=bgFrame.cols;
  FRAME_ROWS=bgFrame.rows;
  
  cvb::CvBlobs preBlobs;
  if (LRVTRACK_SAVE_PROCESSED_VIDEO!="" && !vidPOut.isOpened())
  {
    vidPOut.open("result.avi",
        VIDEO_CODEC,
        VIDEO_FPS,
        bgFrame.size());
    if (!vidPOut.isOpened())
    {
      std::cerr << "Error opening video output file. " <<
        "Problems with video codec. Exiting..." << std::endl;
      exit(5);
    }
  }
  //larvaFit lf;
  double totalFitError=0.0;
  while (!STOP)
  {
    if (!get_next_frame(capture,greyFrame,colorFrame))
      break;
    process_frame(greyFrame,bgFrame,colorFrame,processedFrame);
    imshow("PF",processedFrame);
    waitKey(1);
    IplImage ipl_thresholded = processedFrame;
    labelImg=cvCreateImage(
        cvGetSize(&ipl_thresholded), IPL_DEPTH_LABEL, 1);
    NEW.clear();
    cvLabel(&ipl_thresholded, labelImg, NEW);
    cvb::cvFilterByArea(NEW, 40, 1200);
    cvb::CvBlobs tracked_blobs;
    cvb::CvBlobs blob1;
    cvb::CvBlobs::iterator blIT=NEW.begin();
    for(;blIT!=NEW.end();)
    {
      if(larvaToRing(*blIT->second))
      {
        NEW.erase(blIT++);
      }
      else
        blIT++;
    }

    cvb::CvBlobs::iterator it=NEW.begin();
    if(preBlobs.size()>0)
    {
      newLarvaeTrack(NEW,preBlobs,tracked_blobs);

      if(LRVTRACK_CHOREOGRAPHY_OUTPUT)
      {
        //printSummary(preBlobs,NEW,false);
      }
      preBlobs=tracked_blobs;
      //stringstream sstr;
      //sstr<<lf;
      //BOOST_LOG_TRIVIAL(debug) << sstr;
      //unsigned long long a=getmsofday();
      //totalFitError+=lf.optimize(detected_larvae[LRVTRACK_VERBOSE_LEVEL].blobs.back());
      //BOOST_LOG_TRIVIAL(debug) << "FitTime: " << getmsofday()-a << "micros";
      //lf.showContour();
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
      //initialize(detected_larvae[3]);
      LARVAE_COUNT=preBlobs.size();
      //lf.setup(detected_larvae[LRVTRACK_VERBOSE_LEVEL]);
      //lf.setloops(LRVTRACK_DSTEP,LRVTRACK_WSTEP);
      //lf.optimize(detected_larvae[LRVTRACK_VERBOSE_LEVEL].blobs.back());
    }

    //dumpDetectedLarvae();
    if(LRVTRACK_CHOREOGRAPHY_OUTPUT)
    {
      //printSummary(preBlobs,NEW,true);
    }


    if (SHOWTAGS!=0)
    {
      showTags(tracked_blobs);
    }

    cvReleaseImage(&labelImg);
    //imshow("Extracted Frame",colorFrame);
    //vidPOut << colorFrame;
    //handleKeys(STEP,STOP,TRACK,SHOWTAGS);
    std::cerr << "F: " << CURRENT_FRAME << "/" << TOTAL_FRAMES << "\r";
  }
  std::cerr << endl;
  BOOST_LOG_TRIVIAL(debug) << "Total Fit Error:" << totalFitError ;
  
  //Second Pass
  secondPass();

  if(setup_capture_input(capture) == -1)
  {
    BOOST_LOG_TRIVIAL(error) 
      << "Error setting up the capture device (camera or video file)";
    exit(1);
  }
  CURRENT_FRAME=0;
  while (!STOP)
  {
    if (!get_next_frame(capture,greyFrame,colorFrame))
      break;
    apply_model(); 
  }

  // Third pass
  if(setup_capture_input(capture) == -1)
  {
    BOOST_LOG_TRIVIAL(error) 
      << "Error setting up the capture device (camera or video file)";
    exit(1);
  }
  CURRENT_FRAME=0;
  while (!STOP)
  {
    if (!get_next_frame(capture,greyFrame,colorFrame))
      break;
    /*cv::circle(colorFrame,
               Point(circles[bestCircle][0],circles[bestCircle][1]),
               circles[bestCircle][2],
               Scalar(255,100,100));*/
    if (SHOWTAGS!=0)
    {
      showTags2();
    }
    
    stringstream d;
    d<<CURRENT_FRAME;
    putText(colorFrame,
        d.str(),
        Point2f(20,50),
        FONT_HERSHEY_PLAIN,
        2.0,
        Scalar(0,0,0),
        2,
        CV_AA);

    stringstream n;
    vector<int> compression_params;
    compression_params.push_back(CV_IMWRITE_PNG_COMPRESSION);
    compression_params.push_back(0);
    n << "img/ff" << CURRENT_FRAME << ".png";
    try {
      imwrite(n.str(), colorFrame, compression_params);
    }
    catch (runtime_error& ex) {
      fprintf(stderr, 
          "Exception converting image to PNG format: %s\n", 
          ex.what());
      return 1;
    }
    //vidPOut << colorFrame;
  }
}
