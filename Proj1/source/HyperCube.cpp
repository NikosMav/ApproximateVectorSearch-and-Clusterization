#include "../headers/HyperCube.h"
#include "../headers/Hashtable.h"
#include "../headers/Aux.h"
#include "../headers/distFunc.h"
#include <cmath>
#include <random>
#include <algorithm>
#include <queue>
#include <set>

using namespace std;

HyperCube::HyperCube(unsigned int k, unsigned int hModulo, unsigned int d) {
    int numberOfBuckets;
    vector<unsigned char> *currentFunction;

    random_device rd;
    mt19937 gen(rd());

    this->k = k;
    this->d = d;
    this->hModulo = hModulo;

    //copy w from hashtable static
    this->sCeil = Hashtable::sCeil;

    //we must initialize the buckets needed, which will be 2^d in count
    numberOfBuckets = pow(2, k);
    this->bucketArray = new vector<Bucket *>(numberOfBuckets, NULL);

    //actually initialize the buckets
    int bucketNum;
    for(bucketNum = 0; bucketNum < numberOfBuckets; bucketNum++){
        this->bucketArray->at(bucketNum) = new Bucket();
    }

    //now we need to create the projection functions, which are just 
    //arrays with values 0 or 1
    this->projectionFunctions = new vector<vector <unsigned char>*>(k, NULL);

    //for every "function", we must create a vector of randomly distributed 0's and 1's
    int functionNumber;
    for(functionNumber = 0; functionNumber < k; functionNumber ++){
        //allocate the vector
        currentFunction = new vector<unsigned char>(hModulo);
        this->projectionFunctions->at(functionNumber) = currentFunction;

        //now fill half the vector with 1's
        fill_n(currentFunction->begin(), currentFunction->size() / 2, 1);

        //now shuffle the 1's
        shuffle(currentFunction->begin(), currentFunction->end(), gen);
    }





    //now we define the h function via the disturbance Vectors, as in hashtable
    //make random generator 

    uniform_real_distribution<> dist(0, this->sCeil);

    vector<float> *currentDisturbanceVector = NULL;

    //initialize the vector that will keep the disturbance vectors
    this->disturbanceVectors = new vector<vector<float> *>(this->k, NULL);

    //now we can call dist for s_i
    int disturbanceVectorCount;
    int currentDimension;
    //for every vector we need to make
    for(disturbanceVectorCount = 0; disturbanceVectorCount < this->k; disturbanceVectorCount++){
        //initialize the vector
        currentDisturbanceVector = new vector<float>(this->d, 0.0);

        //fill the values of the current disturbance vector with uniform distribution
        for(currentDimension = 0; currentDimension < this->d; currentDimension++){
            currentDisturbanceVector->at(currentDimension) = dist(gen);
        }

        //assign the created vector to the parent vector
        this->disturbanceVectors->at(disturbanceVectorCount) = currentDisturbanceVector;

    }


}

HyperCube::~HyperCube()
{
    for(Bucket *currBucket : *this->bucketArray){
        delete currBucket;
    }
    delete this->bucketArray;

    for(vector<unsigned char> *currFunction : *this->projectionFunctions){
        delete currFunction;
    }
    delete this->projectionFunctions;

    for(vector<float> *disturbanceVector : *this->disturbanceVectors){
        delete disturbanceVector;
    }
    delete this->disturbanceVectors;

}



void HyperCube::insert(Image *imageToInsert){
    unsigned int hashResult;
    //get the vertex where we will insert
    hashResult = this->hashFunction(imageToInsert);
    //insert to that bucket
    this->bucketArray->at(hashResult)->bucketImages->push_back(imageToInsert);
}


unsigned int HyperCube::hashFunction(Image *image){
    

    unsigned int currH;
    int functionNumber;
    vector<float> *disturbanceVector;
    unsigned char currBit;
    unsigned int concatenationOfBits = 0;


    //for every disturbance vector
    for(functionNumber = 0; functionNumber < this->k; functionNumber ++) {
        disturbanceVector = this->disturbanceVectors->at(functionNumber);
        currH = this->calculateH(image, disturbanceVector);
        //now we must pass the result through the corresponding projection function
        currBit = this->projectionFunctions->at(functionNumber)->at(currH);
        //now add this bit to the result
        concatenationOfBits |= currBit;
        concatenationOfBits <<= 1;
    }

    //last one was unnecessary
    concatenationOfBits >>= 1;

    return concatenationOfBits;

}



unsigned int HyperCube::calculateH(Image *image, vector<float> *const disturbanceVector){
    
    vector<unsigned int> *currentNormalizedVector = NULL;
    int currentPixel, currentDimension;
    unsigned int hashSum;
    int dimension;
    unsigned int d;
    unsigned int *mModulos = this->mPowerModulos;
    unsigned int currAmodM;
    unsigned int M;
    unsigned int moduloMask;
    unsigned int currH;
    float currImageValue;
    float currDisturbanceValue;
    float subtraction;
    vector<unsigned char> * imagePixels = image->imagePixels;

    //initialize vector that will store a_i
    currentNormalizedVector = new vector<unsigned int>(this->d, 0);

    //for every dimension
    for(currentPixel = 0; currentPixel < this->d; currentPixel++){
        currImageValue = (float) imagePixels->at(currentPixel);
        currDisturbanceValue = disturbanceVector->at(currentPixel);

        if(currImageValue >= currDisturbanceValue){
            subtraction = currImageValue - currDisturbanceValue;
        } else {
            subtraction = currDisturbanceValue - currImageValue;
        }

        float division = float(subtraction / this->sCeil);
        currentNormalizedVector->at(currentPixel) = floor(division);
    }

    //given the a vector (currentNormalizedVector), calculate the hash function with math tricks
    hashSum = 0;

    d = this->d;

    M = this->hModulo;
    moduloMask = M - 1;

    //for every dimension
    for(dimension = 0; dimension < d; dimension++){
        // get a_i mod M
        currAmodM = currentNormalizedVector->at(d - dimension - 1) & moduloMask;
        // add to hashSum keeping mod rule true
        hashSum = (hashSum + (currAmodM * mModulos[dimension]) & moduloMask) & moduloMask;
    }

    delete currentNormalizedVector;

    return hashSum;
}


Image *HyperCube::NearestNeighbor(Image *queryImage, unsigned int numberOfPoints, unsigned int maxVertProbed){
    unsigned int currentNode;
    vector<Image *> *currentNeighbors;
    vector<unsigned int> *newNodes;
    Image *currentNeighbor;
    Image *bestNeighbor = NULL;
    unsigned int bestDistance;
    unsigned int currDistance;
    unsigned int verticesCount = pow(2, this->k);
    //we will keep a queue for the BFS traversal of the nodes of the hypercube
    queue<unsigned int> *nodesToVisit = new queue<unsigned int>();
    //also keep a visited set, so BFS does not loop itself
    bool *visited = (bool *) malloc(verticesCount * sizeof(bool));
    unsigned int i;
    for(i = 0; i < verticesCount; i++){
        visited[i] = false;
    }

    //now we must start from the node that the queryImage hashes to
    //and BFS our way out, so long as the stopping criteria don't hold true
    
    unsigned int queryImageHash = this->hashFunction(queryImage);
    //start on hashed node
    visited[queryImageHash] = true;
    nodesToVisit->push(queryImageHash);

    //for stoppping criteria
    unsigned int verticesProbed = 0;
    unsigned int pointsChecked = 0;

    //while we should keep searching
    while(!nodesToVisit->empty()){
        //take the next node to visit
        currentNode = nodesToVisit->front();
        //remove from queue
        nodesToVisit->pop();

        //now look through the neighbors of current node
        currentNeighbors = this->bucketArray->at(currentNode)->bucketImages;

        for(const auto& currentNeighbor: *currentNeighbors) {
            //if image we check is not query image
            if(currentNeighbor != queryImage){
                //if first neighbor
                if(bestNeighbor == NULL){
                    bestNeighbor = currentNeighbor;
                    bestDistance = DIST_FUNC(queryImage, bestNeighbor);

                    pointsChecked++;
                } else {
                    //calculate new distance
                    currDistance = DIST_FUNC(queryImage, currentNeighbor);

                    //if better, update best to current
                    if(currDistance < bestDistance){
                        bestDistance = currDistance;
                        bestNeighbor = currentNeighbor;
                    }

                    //stopping criteria
                    pointsChecked++;
                    if(pointsChecked >= numberOfPoints){
                        return bestNeighbor;
                    }
                }

            }
        }

        //stopping criteria
        verticesProbed++;
        if(verticesProbed >= maxVertProbed) {
            return bestNeighbor;
        }

        //now we must add the next nodes to visit from here to the queue
        //only if they have not been visited yet
        newNodes = this->hammingNeighbors(currentNode);

        //add not visited nodes to queue and mark as visited
        for(auto const& newNode : *newNodes){
            if(!visited[newNode]){
                visited[newNode] = true;
                nodesToVisit->push(newNode);
            }
        }
        delete newNodes;
    }

    delete nodesToVisit;
    free(visited);


    return bestNeighbor;
}

vector<Image *> *HyperCube::kNearestNeighbors(Image *queryImage, unsigned int k, unsigned int numberOfPoints, unsigned int maxVertProbed){
    unsigned int currentNode;
    vector<Image *> *currentNeighbors;
    vector<unsigned int> *newNodes;
    Image *currentNeighbor;
    unsigned int currDistance;
    unsigned int verticesCount = pow(2, this->k);
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

    //we will keep a queue for the BFS traversal of the nodes of the hypercube
    queue<unsigned int> *nodesToVisit = new queue<unsigned int>();

    //also keep a visited set, so BFS does not loop itself
    bool *visited = (bool *) malloc(verticesCount * sizeof(bool));
    bool toBreak = false;
    unsigned int i;
    for(i = 0; i < verticesCount; i++){
        visited[i] = false;
    }

    //now we must start from the node that the queryImage hashes to
    //and BFS our way out, so long as the stopping criteria don't hold true
    
    unsigned int queryImageHash = this->hashFunction(queryImage);
    //start on hashed node
    visited[queryImageHash] = true;
    nodesToVisit->push(queryImageHash);

    //for stoppping criteria
    unsigned int verticesProbed = 0;
    unsigned int pointsChecked = 0;

    //while we should keep searching
    while(!nodesToVisit->empty()){
        if(toBreak){
            break;
        }

        //take the next node to visit
        currentNode = nodesToVisit->front();
        //remove from queue
        nodesToVisit->pop();

        //now look through the neighbors of current node
        currentNeighbors = this->bucketArray->at(currentNode)->bucketImages;

        for(const auto& currentNeighbor: *currentNeighbors) {
            //if image we check is not query image
            if(currentNeighbor != queryImage){
                
                //Initializing new struct
                candidateImage candidate;
                setCandidateImage(&candidate, currentNeighbor, queryImage);

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

                //stopping criteria
                pointsChecked++;
                if(pointsChecked >= numberOfPoints){
                    toBreak = true;
                    break;
                }
            }
        }

        //stopping criteria
        verticesProbed++;
        if(verticesProbed >= maxVertProbed) {
            break;
        }

        //now we must add the next nodes to visit from here to the queue
        //only if they have not been visited yet
        newNodes = this->hammingNeighbors(currentNode);

        //add not visited nodes to queue and mark as visited
        for(auto const& newNode : *newNodes){
            if(!visited[newNode]){
                visited[newNode] = true;
                nodesToVisit->push(newNode);
            }
        }
        delete newNodes;
    }

    delete nodesToVisit;
    free(visited);

    while (!pq.empty()) {
        nearestNeighbors->push_back(pq.top().candidateImage);
        pq.pop();
    }
    reverse(nearestNeighbors->begin(), nearestNeighbors->end());

    return nearestNeighbors;
}


vector<Image *> *HyperCube::rangeSearch(Image *queryImage, double R, unsigned int numberOfPoints, unsigned int maxVertProbed, set<unsigned int> *ignoredPoints) {
    unsigned int currentNode;
    vector<Image *> *currentNeighbors;
    vector<unsigned int> *newNodes;
    Image *currentNeighbor;
    unsigned int currDistance;
    unsigned int verticesCount = pow(2, this->k);

    //Initializing the final vector of query's nearest neighbors 
    vector<Image *> *rangeNeighbors = new vector<Image *>();

    //we will keep a queue for the BFS traversal of the nodes of the hypercube
    queue<unsigned int> *nodesToVisit = new queue<unsigned int>();

    //also keep a visited set, so BFS does not loop itself
    bool *visited = (bool *) malloc(verticesCount * sizeof(bool));
    bool toBreak = false;
    unsigned int i;
    for(i = 0; i < verticesCount; i++){
        visited[i] = false;
    }

    //now we must start from the node that the queryImage hashes to
    //and BFS our way out, so long as the stopping criteria don't hold true
    
    unsigned int queryImageHash = this->hashFunction(queryImage);
    //start on hashed node
    visited[queryImageHash] = true;
    nodesToVisit->push(queryImageHash);

    //for stoppping criteria
    unsigned int verticesProbed = 0;
    unsigned int pointsChecked = 0;

    //while we should keep searching
    while(!nodesToVisit->empty()){
        if(toBreak){
            break;
        }

        //take the next node to visit
        currentNode = nodesToVisit->front();
        //remove from queue
        nodesToVisit->pop();

        //now look through the neighbors of current node
        currentNeighbors = this->bucketArray->at(currentNode)->bucketImages;

        for(const auto& currentNeighbor: *currentNeighbors) {

            //if we have ignored points
            if(ignoredPoints != NULL){
                //Check if current candidate image has already been returned from older range search
                if(ignoredPoints->find(currentNeighbor->id) != ignoredPoints->end()){
                    continue;
                }
            }
            //if image we check is not query image
            if(currentNeighbor != queryImage){
                
                //get candidate's distance
                currDistance = DIST_FUNC(queryImage, currentNeighbor);
                //check if distance is in range
                if(currDistance < R){
                    rangeNeighbors->push_back(currentNeighbor);
                }

                //stopping criteria
                pointsChecked++;
                if(pointsChecked >= numberOfPoints){
                    toBreak = true;
                    break;
                }
            }
        }

        //stopping criteria
        verticesProbed++;
        if(verticesProbed >= maxVertProbed) {
            break;
        }

        //now we must add the next nodes to visit from here to the queue
        //only if they have not been visited yet
        newNodes = this->hammingNeighbors(currentNode);

        //add not visited nodes to queue and mark as visited
        for(auto const& newNode : *newNodes){
            if(!visited[newNode]){
                visited[newNode] = true;
                nodesToVisit->push(newNode);
            }
        }
        delete newNodes;
    }

    delete nodesToVisit;
    free(visited);


    return rangeNeighbors;
}

//produces the k numbers which have hamming distance 1 from the given
vector<unsigned int> *HyperCube::hammingNeighbors(unsigned int node){
    vector<unsigned int> *results = new vector<unsigned int>(this->k, 0);
    unsigned int i;

    //for every bit to flip
    for(i = 0; i < this->k; i++){
        //add to results the node but flip bit i
        results->push_back(node ^ (1 << i));
    }

    return results;

}


//this takes some centroids and uses the range search to implement reverse assignment
vector<set <Image *> *> *HyperCube::clusterizeHyperCube(vector<Image *> *centroids, vector<Image *> *inputImages, unsigned int numberOfPoints, unsigned int maxVertProbed){
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
            fetchedImages = this->rangeSearch(currentCentroid, currentRadius, numberOfPoints, maxVertProbed, imagesToIgnoreInRangeSearch);

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



unsigned int *HyperCube::mPowerModulos;