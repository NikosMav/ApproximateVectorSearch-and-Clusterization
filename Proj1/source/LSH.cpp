#include "../headers/LSH.h"
#include "../headers/Hashtable.h"
#include "../headers/Aux.h"
#include "../headers/distFunc.h"

#include <set>
#include <algorithm>

using namespace std;

LSH::LSH(unsigned int L, int numberOfBuckets){
    //keep L so we know it just in case
    this->L = L;

    //we must make a vector that keeps L different hashtables
    this->hashtables = new vector<Hashtable *>(this->L, NULL);
    
    int currHashtable;

    for(currHashtable = 0; currHashtable < this->L; currHashtable++){
        this->hashtables->at(currHashtable) = new Hashtable(numberOfBuckets);
    }
}

LSH::~LSH(){
    for(Hashtable *hashtable : *this->hashtables){
        delete hashtable;
    }
    delete this->hashtables;
}


//insert image to every hash table
void LSH::insert(Image *image){
    int currHashtable;

    for(currHashtable = 0; currHashtable < this->L; currHashtable++){
        this->hashtables->at(currHashtable)->insert(image);
    }
}

Image *LSH::NearestNeighbor(Image *queryImage){
    //signifies no neighbor found yet
    Image *currBestImage = NULL;
    vector<Image *> *currNeighbors = NULL;
    Hashtable *currHashtable = NULL;

    unsigned int currBestDistance = 0;
    unsigned int newDistance;

    int currHashtableIndex;

    //for every hashtable
    for(currHashtableIndex = 0; currHashtableIndex < this->L; currHashtableIndex ++){
        currHashtable = this->hashtables->at(currHashtableIndex);
        //get neighbors of queryImage in such hashtable
        currNeighbors = currHashtable->getNeighbors(queryImage);

        //iterate through all images in same bucket (including self, so check)
        for(const auto& currCandidateImage: *currNeighbors) {
            //check if not self
            if(currCandidateImage->id != queryImage->id){
                //if no good neighbor found yet
                if(currBestImage == NULL){
                    //set this image as best
                    currBestImage = currCandidateImage;
                    //set the distance to it as best
                    currBestDistance = DIST_FUNC(queryImage, currBestImage);
                } else { // else check against best
                    // get distance to candidate
                    newDistance = DIST_FUNC(queryImage, currCandidateImage);

                    //if better, replace
                    if(newDistance < currBestDistance){
                        currBestImage = currCandidateImage;
                        currBestDistance = newDistance;
                    }
                }
            }
        }
    }

    return currBestImage;
}

vector<Image *> *LSH::kNearestNeighbors(Image *queryImage, unsigned int k) {
    vector<Image *> *currNeighbors = NULL;
    Hashtable *currHashtable = NULL;

    int currHashtableIndex;
    unsigned int counter = 0;

    //Initializing the final vector of query's nearest neighbors 
    vector<Image *> *nearestNeighbors = new vector<Image *>();

    //comparator lambda function
    auto comp = [](candidateImage a, candidateImage b) {
        //comparison logic
        return(a.score < b.score);
    };

    //Initializing priority queue
    priority_queue<candidateImage, vector<candidateImage>, decltype(comp)> pq(comp);

    //keep set of checked images for no duplicates
    set<unsigned int> *checkedImages = new set<unsigned int>();

    //for every hashtable
    for(currHashtableIndex = 0; currHashtableIndex < this->L; currHashtableIndex ++){
        currHashtable = this->hashtables->at(currHashtableIndex);
        //get neighbors of queryImage in such hashtable
        currNeighbors = currHashtable->getNeighbors(queryImage);

        //iterate through all images in same bucket (including self, so check)
        for(const auto& currCandidateImage: *currNeighbors) {

            //check if not self or already checked
            if(currCandidateImage->id != queryImage->id && checkedImages->find(currCandidateImage->id) == checkedImages->end()){
                
                //we checked now
                checkedImages->insert(currCandidateImage->id);

                //Initializing new struct
                candidateImage candidate;
                setCandidateImage(&candidate, currCandidateImage, queryImage);

                if(counter < k){
                    pq.push(candidate);
                    counter++;

                } else{
                    //checking the queue header
                    if(pq.top().score > candidate.score){
                        pq.pop();
                        pq.push(candidate);
                    }
                }         
            }
        }
    }


    delete checkedImages;

    while (!pq.empty()) {
        nearestNeighbors->push_back(pq.top().candidateImage);
        pq.pop();
    }

    reverse(nearestNeighbors->begin(), nearestNeighbors->end());

    return nearestNeighbors;
}


vector<Image *> *LSH::rangeSearch(Image *queryImage, double R, set<unsigned int> *ignoredPoints) {
    vector<Image *> *currNeighbors = NULL;
    Hashtable *currHashtable = NULL;

    unsigned int currDistance = 0;
    int currHashtableIndex;

    //Initializing the final vector of query's nearest neighbors 
    vector<Image *> *rangeNeighbors = new vector<Image *>();

    //keep set of checked images for no duplicates
    set<unsigned int> *checkedImages = new set<unsigned int>();

    //for every hashtable
    for(currHashtableIndex = 0; currHashtableIndex < this->L; currHashtableIndex ++){
        currHashtable = this->hashtables->at(currHashtableIndex);
        //get neighbors of queryImage in such hashtable
        currNeighbors = currHashtable->getNeighbors(queryImage);

        //iterate through all images in same bucket (including self, so check)
        for(const auto& currCandidateImage: *currNeighbors) {

            //if we have ignored points
            if(ignoredPoints != NULL){
                //Check if current candidate image has already been returned from older range search
                if(ignoredPoints->find(currCandidateImage->id) != ignoredPoints->end()){
                    continue;
                }
            }

            //check if not self or already checked
            if(currCandidateImage->id != queryImage->id && checkedImages->find(currCandidateImage->id) == checkedImages->end()){
                //we checked now
                checkedImages->insert(currCandidateImage->id);

                //get candidate's distance
                currDistance = DIST_FUNC(queryImage, currCandidateImage);
                //check if distance is in range
                if(currDistance < R){
                    rangeNeighbors->push_back(currCandidateImage);
                }
            }
        }
    }

    delete checkedImages;

    return rangeNeighbors;
}

//this takes some centroids and uses the range search to implement reverse assignment
vector<set <Image *> *> *LSH::clusterizeLSH(vector<Image *> *centroids, vector<Image *> *inputImages){
    //infer k from the centroids
    int k = centroids->size();

    //infer the number of images
    int numberOfImages = inputImages->size();
    int currentImageNumber;
    int currentCentroidNumber;

    unsigned int newDistance, centroidDistance;

    //initialize empty clusters
    vector<set <Image *> *> *clusters = new vector<set <Image *> *>(k, NULL);
    vector<set <Image *> *> *LloydClusters;
    vector<set <Image *> *> *finalClusters = new vector<set <Image *> *>(k, NULL);
    int clusterNumber;

    //make a new empty set for every cluster
    for(clusterNumber = 0; clusterNumber < k; clusterNumber++){
        clusters->at(clusterNumber) = new set<Image *>();
        finalClusters->at(clusterNumber) = new set<Image *>();
    }    

    vector<Image *> *fetchedImages;
    vector<Image *> *LloydImages = new vector<Image *>();


    double currentRadius;
    double radiusThreshold;

    //now we will find min and max distance between centroids in order to define starting radius and radius threshold
    int minDist = -1, maxDist = -1;
    int i, j;

    for(i = 0; i < k; i++){
        for(j = i + 1; j < k; j++){
            //calculate the distance between two centroids
            centroidDistance = DIST_FUNC(centroids->at(i), centroids->at(j));

            //update minimum distance
            if(minDist == -1){
                minDist = centroidDistance;
            } else if (centroidDistance < minDist){
                minDist = centroidDistance;
            }

            //update maximum distance
            if(maxDist == -1){
                maxDist = centroidDistance;
            } else if (centroidDistance > maxDist){
                maxDist = centroidDistance;
            }
        }
    }


    //set first radius and radius threshold accordingly
    currentRadius = ((double) minDist) / 2.0;
    radiusThreshold = ((double) maxDist) / 2.0;




    //we must initialize the structs for every image, that will hold data useful to the clustering process
    imageClusterData *imageClusteringData = (imageClusterData *) malloc(numberOfImages * sizeof(imageClusterData));

    set<unsigned int> *imagesToSetAfterConsidering;
    set<unsigned int> *imagesToIgnoreInRangeSearch = new set<unsigned int>();
    //this is used for stopping the loop
    int newImagesClustered;

    //initialize values in the array
    for(currentImageNumber = 0; currentImageNumber < numberOfImages; currentImageNumber ++){
        //not set to cluster yet
        imageClusteringData[currentImageNumber].assignedCluster = -1;
        //not considering yet
        imageClusteringData[currentImageNumber].considering = false;
        //also not set yet
        imageClusteringData[currentImageNumber].set = false;
        //also no min distance yet
        imageClusteringData[currentImageNumber].runningMinDist = -1;
    }


    //starting looping on radius
    do{
        currentCentroidNumber = 0;
        //initialize the batch of considering images to set
        imagesToSetAfterConsidering = new set<unsigned int>();
        //now, we must loop through all centroids
        for(Image *currentCentroid : *centroids){
            //we do a range search using the LSH
            fetchedImages = this->rangeSearch(currentCentroid, currentRadius, imagesToIgnoreInRangeSearch);

            //we must loop through all fetched images, and possibly assign some to current centroid's cluster
            for(Image *currentFetchedImage : *fetchedImages){
                //we must check if this image is not set yet
                if(!imageClusteringData[currentFetchedImage->id].set){
                    
                    //if not considering yet, virtually add to cluster and start considering
                    if(!imageClusteringData[currentFetchedImage->id].considering){
                        //start considering this image for assignment to cluster
                        imageClusteringData[currentFetchedImage->id].considering = true;
                        //set current (may change later) cluster of image
                        imageClusteringData[currentFetchedImage->id].assignedCluster = currentCentroidNumber;
                        //set minimum distance, in case we need to change cluster so we must compare
                        imageClusteringData[currentFetchedImage->id].runningMinDist = DIST_FUNC(currentFetchedImage, currentCentroid);
                        
                        //since we started considering, add to helper set in order to "set" the image after the current radius
                        imagesToSetAfterConsidering->insert(currentFetchedImage->id);

                    } else { //if not, see if current centroid is better than last
                        //check if current distance is better than last
                        newDistance = DIST_FUNC(currentFetchedImage, currentCentroid);
                        if(newDistance < imageClusteringData[currentFetchedImage->id].runningMinDist){
                            //update clustering data for this image

                            //set current (may change later) cluster of image
                            imageClusteringData[currentFetchedImage->id].assignedCluster = currentCentroidNumber;
                            //set minimum distance, in case we need to change cluster so we must compare
                            imageClusteringData[currentFetchedImage->id].runningMinDist = newDistance;
                        }
                        //else we do nothing
                    }

                }
            }

            delete fetchedImages;

            //advance centroid number
            currentCentroidNumber ++;
        }

        //count the new images
        newImagesClustered = imagesToSetAfterConsidering->size();

        //now, since we are done with current radius, we must turn all "considering" to "set"
        //seeing as how they will always be retrieved with an even bigger radius from now on
        for(unsigned int imageToSet : *imagesToSetAfterConsidering){
            imageClusteringData[imageToSet].considering = false;
            imageClusteringData[imageToSet].set = true;

            //add to ignore set for faster range search
            imagesToIgnoreInRangeSearch->insert(imageToSet);
        }

        delete imagesToSetAfterConsidering;

        //now double radius
        currentRadius *= 2;
    }while(newImagesClustered > numberOfImages / 100 || currentRadius < radiusThreshold); //we stop when we acquire too few new images into clusters



    //we should actually insert the set images to the clusters
    for(currentImageNumber = 0; currentImageNumber < numberOfImages; currentImageNumber ++){
        //if image is set
        if(imageClusteringData[currentImageNumber].set){
            //actually add to cluster
            clusters->at(imageClusteringData[currentImageNumber].assignedCluster)->insert(inputImages->at(currentImageNumber));
        } else {
            //add to batch that will be executed via Lloyd's
            LloydImages->push_back(inputImages->at(currentImageNumber));
        }
    }

    //after the approximate part of the algorithm, we must clean up
    LloydClusters = LloydsAlgorithm(centroids, LloydImages);



    //union all corresponding clusters
    for(clusterNumber = 0; clusterNumber < k; clusterNumber++){
        set_union(clusters->at(clusterNumber)->begin(), clusters->at(clusterNumber)->end(),
                  LloydClusters->at(clusterNumber)->begin(), LloydClusters->at(clusterNumber)->end(),
                  inserter(*finalClusters->at(clusterNumber), finalClusters->at(clusterNumber)->begin()));
    }    

    for(set<Image *> *cluster : *clusters){
        delete cluster;
    }
    delete clusters;

    for(set<Image *> *cluster : *LloydClusters){
        delete cluster;
    }
    delete LloydClusters;


    delete LloydImages;
    delete imagesToIgnoreInRangeSearch;
    free(imageClusteringData);

    return finalClusters;
}