#include "../headers/Aux.h"
#include "../headers/Hashtable.h"
#include "../headers/Image.h"
#include "../headers/LSH.h"
#include "../headers/HyperCube.h"
#include "../headers/distFunc.h"
#include <bits/stdc++.h>
#include <cmath>
#include <arpa/inet.h>
#include <fstream>
#include <arpa/inet.h>
#include <string>
#include <vector>
#include <iterator>
#include <random>
#include <set>
#include <chrono>

//gives the argument after the one given
//e.g. getOptionValue("-p", ...) would give the file name
//that the user provided in the command line
string *getOptionValue(string option, int argc, char **argv){
    string *value = NULL;
    int i;

    for(i = 1; i < argc - 1; i++){
        if(option == argv[i]){
            value = new string(); 
            *value = argv[i + 1];
            break;
        }
    }
    return value;
}

bool getOption(string option, int argc, char **argv){
    int i;

    for(i = 1; i < argc; i++){
        if(option == argv[i]){
            return true;
        }
    }
    return false;
}


unsigned int parseInteger(fstream *binaryFile){
    unsigned int output;

    binaryFile->read((char *) &output, 4);

    output = htonl(output);

    return output;

}


Image *parseImage(fstream *binaryFile, unsigned int pixels){
    static int idOffset = 0;

    unsigned char imageArray[pixels];

    binaryFile->read((char *) imageArray, pixels);
    
    vector<unsigned char> *imageData = new vector<unsigned char>(imageArray, imageArray + pixels);
    
    Image *returnValue = new Image(idOffset, imageData);

    idOffset++;

    return returnValue;
}


vector<Image *> *parseImages(fstream *binaryFile, unsigned int numberOfImages, unsigned int pixels){
    //initalize vector which will keep pointers to images with NULL and enough slots
    vector<Image *> *images = new vector<Image *>(numberOfImages, NULL);
    
    int imageCount;

    for(imageCount = 0; imageCount < numberOfImages; imageCount++){
        images->at(imageCount) = parseImage(binaryFile, pixels);
    }

    return images;
}

//given initialized m, d and M, calculates and stores m^(i) mod M
void initializeModulosForHashTable(){
    unsigned int d = Hashtable::d;
    unsigned int m = Hashtable::m;
    unsigned int M = Hashtable::hModulo;
    unsigned int moduloMask = M - 1;
    unsigned int multFactor;
    unsigned int *modArray;
    int power;

    Hashtable::mPowerModulos = (unsigned int *) malloc(d * sizeof(unsigned int));

    modArray = Hashtable::mPowerModulos;

    //modArray[i] = m^i mod M
    modArray[0] = 1;
    multFactor = m & moduloMask;

    //for the rest slots in the array
    for(power = 1; power < d; power++){
        //we know m^i mod M = (m^(i - 1) * m) % M = ((m^(i - 1) % M) * (m % M)) % M 
        modArray[power] = (modArray[power - 1] * multFactor) & moduloMask;
    }
    //now the array is initialized

    HyperCube::mPowerModulos = Hashtable::mPowerModulos;
    
}


void setCandidateImage(candidateImage* candidate, Image* currCandidateImage, Image* queryImage) {
    
    candidate->candidateImage = currCandidateImage;
    //Calculating current distance and setting it to struct
    unsigned int currDistance = DIST_FUNC(queryImage, currCandidateImage);
    candidate->score = currDistance;
}


vector<Image *> *initialCentroidImages(vector<Image *> *images, unsigned int k){
    //make the random generators
    random_device rd;
    mt19937 e2(rd());

    //this is the return value
    vector<Image *> *chosenImages = new vector<Image *>(k, NULL);


    unsigned int imagesCount = images->size();

    //also remember the chosen ones
    bool *isChosen = (bool *) malloc(imagesCount * sizeof(bool));
    unsigned int chosenId;
    for(chosenId = 0; chosenId < imagesCount; chosenId ++){
        isChosen[chosenId] = false;
    }


    //initialize to 0, which means not set yet
    vector<unsigned int> *minDistances = new vector<unsigned int>(imagesCount, 0);
    //this is P(i)
    vector<double> *partialSums = new vector<double>(imagesCount + 1, 0.0);

    //choose a random first centroid
    srand(time(NULL));
    unsigned int firstCentroidId = rand() % imagesCount;

    //add to chosen images 
    chosenImages->at(0) = images->at(firstCentroidId);
    isChosen[firstCentroidId] = true;

    unsigned int currCentroidNumber;
    
    unsigned int maxDistance = 0;

    //loop for rest of centroids to pick
    for(currCentroidNumber = 1; currCentroidNumber < k; currCentroidNumber ++){
        
        //consider the changes to D(i) because of the last centroid chosen
        updateMinDistancesToCentroids(images, isChosen, chosenImages, currCentroidNumber, minDistances, &maxDistance);
        
        //calculate P(i) based on D(i)
        calculatePartialSums(partialSums, imagesCount, isChosen, minDistances, maxDistance);

        //pick the next centroid
        chooseNextCentroid(e2, partialSums, images, chosenImages, currCentroidNumber, isChosen);
    }

    free(isChosen);
    delete minDistances;
    delete partialSums;

    return chosenImages;

}


//this keeps the D(i) table up to date
void updateMinDistancesToCentroids(vector<Image *> *images, bool *isChosen, vector<Image *> *chosenImages, unsigned int currCentroidNumber, vector<unsigned int> *minDistances, unsigned int *maxDistance){
    unsigned int currImageId;
    unsigned int newDistance;

    unsigned int imagesCount = images->size();

    //possibly update min distances according to new centroid
    for(currImageId = 0; currImageId < imagesCount; currImageId++){
        //if not already a centroid
        if(!isChosen[currImageId]){

            //get distance to previously chosen centroid
            newDistance = DIST_FUNC(images->at(currImageId), chosenImages->at(currCentroidNumber - 1));

            //check if this distance is smaller than previously min, or if not set at all
            if(newDistance < minDistances->at(currImageId) || minDistances->at(currImageId) == 0){
                //replace distance
                minDistances->at(currImageId) = newDistance;

                //we should also consider updating maxDistance
                if(newDistance > *maxDistance){
                    *maxDistance = newDistance;
                }
            }
        }
    }
}


void calculatePartialSums(vector<double> *partialSums, unsigned int imagesCount, bool *isChosen, vector<unsigned int> *minDistances, unsigned int maxDistance){
    int partialSumNumber;
    double quotient;

    //now that we have updated the distances, we must calculate the partial sums
    for(partialSumNumber = 1; partialSumNumber <= imagesCount; partialSumNumber++){ //starting from 1 since P(0) = 0
        //P(i) = P(i - 1) + (D(i - 1) / maxD(j)) ^ 2
        if(isChosen[partialSumNumber - 1]){
            quotient = 0.0;
        } else{
            quotient = (((double) minDistances->at(partialSumNumber - 1)) / ((double) maxDistance));
        }
        partialSums->at(partialSumNumber) = partialSums->at(partialSumNumber - 1) + quotient * quotient;
    }
}



void chooseNextCentroid(mt19937 &e2, vector<double> *partialSums, vector<Image *> *images, vector<Image *> *chosenImages, unsigned int currCentroidNumber, bool *isChosen){
    unsigned int imagesCount = images->size();

    //now we must roll the dice and select the next centroid
    //define the range for the random number
    uniform_real_distribution<> dist(0, partialSums->at(imagesCount));

    //pick the random number in the range
    double diceRoll = dist(e2);        


    int left, right, binarySearchIndex;

    double previousValue, nextValue;

    //now we must see where this number falls in the partial sums so we pick the next centroid
    //start in the middle
    left = 0;
    right = imagesCount;

    do{
        binarySearchIndex = left + (right - left) / 2; 

        previousValue = partialSums->at(binarySearchIndex);
        nextValue = partialSums->at(binarySearchIndex + 1);
        //if we find where the dice falls
        if(previousValue < diceRoll && diceRoll <= nextValue){
            //get out of the loop, we found the next centroid
            chosenImages->at(currCentroidNumber) = images->at(binarySearchIndex);
            isChosen[binarySearchIndex] = true;
            break;
        }

        //else go on with the search

        //if we must go more left
        if(previousValue >= diceRoll){
            right = binarySearchIndex - 1;
        } else { //else more right
            left = binarySearchIndex + 1;
        }

    }while(true);
}



vector<set<Image *> *> * LloydsAlgorithm(vector<Image *> *centroids, vector<Image *> *images){
    Image *currentImage;
    Image *currentBestCentroid;
    
    unsigned int currentDistance;
    unsigned int bestDistance;
    unsigned int index, bestIndex;
    int numberOfCentroids = centroids->size();

    //initialize cluster vector to be returned
    vector<set<Image *> *> *clusters = new vector<set<Image *> *>(numberOfCentroids, NULL);
    for(index = 0; index < numberOfCentroids; index++){
        //initialize each cluster's set
        set<Image *> *clustersSet = new set<Image *>();
        //setting each cluster's set
        clusters->at(index) = clustersSet;
    }

    //iterating through images Array
    for(Image* currentImage: *images) {
        currentBestCentroid = NULL;
        //iterating through centroids Array
        for(index = 0; index < numberOfCentroids; index++) {
            //if no best found yet
            if(currentBestCentroid == NULL){
                //set first centroid as best
                bestDistance = DIST_FUNC(centroids->at(index), currentImage);
                currentBestCentroid = centroids->at(index);
                bestIndex = index;
            } else {
                //get candidte's distance
                currentDistance = DIST_FUNC(centroids->at(index), currentImage);

                //if better, replace
                if(currentDistance < bestDistance) {
                    bestDistance = currentDistance;
                    currentBestCentroid = centroids->at(index);
                    bestIndex = index;
                }
            }
        }
        //Now that we found the image's best centroid, we can update the centroid's cluster
        clusters->at(bestIndex)->insert(currentImage);
    }

    return clusters;
}

vector<Image *> *centroidize(vector<set<Image *> *> *clusters) {
    unsigned int imageSize, index, clusterSize, imageIndex, numberOfClusters, clusterIndex;
    //get images' size
    imageSize = (*clusters->at(0)->begin())->imagePixels->size();
    //Initialize transpose vector
    vector<vector<unsigned char> *> *transposeMatrix; 

    //Array of newCentroids has to have same size as the number of all clusters 
    numberOfClusters = clusters->size();
    vector<vector<unsigned char> *> *newCentroids = new vector<vector<unsigned char> *>(numberOfClusters, NULL);

    //Initializing final product with same size as newCentroids array
    vector<Image *> *centroidImages = new vector<Image *>(numberOfClusters, NULL);

    clusterIndex = 0;
    //For each cluster
    for(const auto& cluster: *clusters){
        //initialize 1st dimension 
        transposeMatrix = new vector<vector<unsigned char> *>(imageSize, NULL);
        //initialize 2nd dimension
        clusterSize = cluster->size();
        for(index = 0; index < imageSize; index++){
            transposeMatrix->at(index) = new vector<unsigned char>(clusterSize, 0);
        }

        imageIndex = 0;
        //iterate through it's set
        for(Image *setImage: *cluster){
            //for each image's dimension

            for(index = 0; index < imageSize; index++){
                transposeMatrix->at(index)->at(imageIndex) = setImage->imagePixels->at(index);

            }
            imageIndex++;
        }

        //sort every row of new transpose Matrix
        for(vector<unsigned char> *row: *transposeMatrix){
            sort(row->begin(), row->end());
        }

        newCentroids->at(clusterIndex) = new vector<unsigned char>(imageSize, 0);

        //get values for new centroid
        for(index = 0; index < imageSize; index++){
            if(clusterSize == 0){
                newCentroids->at(clusterIndex)->at(index) = 0;
            } else {
                //ith value is the median of the ith dimension of all images in cluster
                newCentroids->at(clusterIndex)->at(index) = transposeMatrix->at(index)->at(clusterSize / 2);
            }
        }

        clusterIndex++;

        for(vector<unsigned char> *transVector : *transposeMatrix){
            delete transVector;
        }
        delete transposeMatrix;
    }

    //creating final vector of newly made centroids
    for(index = 0; index < numberOfClusters; index++){
        centroidImages->at(index) = new Image(-1, newCentroids->at(index));
    }

    delete newCentroids;
    
    return centroidImages;
}

vector<Image *> *exhaustiveNearestNeighbors(Image *queryImage, vector<Image *> *images, unsigned int k){
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

    //iterate through images Array
    for(const auto& currCandidateImage: *images) {

        if(currCandidateImage->id != queryImage->id) {
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

    while (!pq.empty()) {
        nearestNeighbors->push_back(pq.top().candidateImage);
        pq.pop();
    }

    reverse(nearestNeighbors->begin(), nearestNeighbors->end());

    return nearestNeighbors;
}

void getCommonArguments(int argc, char **argv, string **inputFilePath, string **queryFilePath, string **outputFilePath, string *NStr, string *RStr, unsigned int &N, double &R){

    *inputFilePath = getOptionValue("-d", argc, argv);
    if(*inputFilePath == NULL){
        cout << "Please enter the input file path:" << endl;

        *inputFilePath = new string();
        getline(cin, **inputFilePath);
    }

    *queryFilePath = getOptionValue("-q", argc, argv);
    if(*queryFilePath == NULL){
        cout << "Please enter the query file path:" << endl;

        *queryFilePath = new string();
        getline(cin, **queryFilePath);
    }


    *outputFilePath = getOptionValue("-o", argc, argv);
    if(*outputFilePath == NULL){
        cout << "Please enter the output file path:" << endl;

        *outputFilePath = new string();
        getline(cin, **outputFilePath);
    }

    NStr = getOptionValue("-N", argc, argv);
    if(NStr == NULL){
        N = 1;
    }else{
        N = stoi(*NStr);
        //we don't need NStr anymore
        delete NStr;
    }

    RStr = getOptionValue("-R", argc, argv);
    if(RStr == NULL){
        R = 10000.0;
    }else{
        R = stof(*RStr);
        //we don't need Rstr anymore
        delete RStr;
    }
}


vector<Image *> *openTrainFile(string *inputFilePath, unsigned int &pixels, unsigned int &numberOfImages, fstream **trainFile){
    unsigned int rows, columns;
    vector<Image *> *images;
    
    //open train file in read mode
    *trainFile = new fstream(inputFilePath->c_str(), ios::in | ios::binary);

    //get rid of metadata
    parseInteger(*trainFile);
    //except the number of images, rows, columns
    numberOfImages = parseInteger(*trainFile);
    rows = parseInteger(*trainFile);
    columns = parseInteger(*trainFile);
    pixels = rows * columns;

    //load the train images into memory
    images = parseImages(*trainFile, numberOfImages, pixels);

    return images;
}

vector<Image *> *openQueryFile(string *queryFilePath, unsigned int &numberOfQueries, unsigned int pixels){
    vector<Image *> *queryImages;

    //open query file in read mode
    fstream *queryFile = new fstream(queryFilePath->c_str(), ios::in | ios::binary);

    //get rid of metadata
    parseInteger(queryFile);
    //except the number of queries
    numberOfQueries = parseInteger(queryFile);
    parseInteger(queryFile);
    parseInteger(queryFile);

    //load the query images into memory
    queryImages = parseImages(queryFile, numberOfQueries, pixels);

    delete queryFile;

    return queryImages;
}

LSH* initializeLSH(unsigned int k, unsigned int pixels, unsigned int L, vector<Image *> *images){
    //set static members of Hashtable
    Hashtable::sCeil = W;
    Hashtable::k = k;
    Hashtable::hModulo = pow(2, 32 / k);

    // less than M / 2 but not a power of 2, so that the mPowerModulos won't be 0
    Hashtable::m = Hashtable::hModulo / 2 - 1;
    Hashtable::d = pixels;

    //now initialize m modulos 
    initializeModulosForHashTable();

    //initialize an LSH class
    LSH *myLSH = new LSH(L, images->size() / BUCKET_FACTOR + 1);

    //insert all images
    for(const auto& image: *images) {
        cout << image->id << endl;
        myLSH->insert(image);
    }

    return myLSH;
}

HyperCube* initializeHyperCube(unsigned int k, unsigned int pixels, vector<Image *> *images){
    //set static members of Hashtable
    Hashtable::sCeil = W;
    Hashtable::k = k;
    Hashtable::hModulo = pow(2, 32 / k);

    // less than M / 2 but not a power of 2, so that the mPowerModulos won't be 0
    Hashtable::m = Hashtable::hModulo / 2 - 1;
    Hashtable::d = pixels;

    //now initialize m modulos 
    initializeModulosForHashTable();

    //initialize an HyperCube class
    HyperCube *myHyperCube = new HyperCube(k, Hashtable::hModulo, pixels);


    //insert all images
    for(const auto& image: *images) {
        cout << image->id << endl;
        myHyperCube->insert(image);
    }

    return myHyperCube;
}

//Needed to measure time execution
using namespace std::chrono;

void calculateOutput(string* outputFilePath, unsigned int pixels, unsigned int numberOfImages, vector<Image *> *images, vector<Image *> *queryImages, void* myHasher, bool hasherBoolean, unsigned int k, unsigned int N, double R, unsigned int numberOfPoints, unsigned int maxVertProbed){
    unsigned int numberOfQueries, index;
    int offset;
    vector<Image *> *kNearestNeighbors;
    vector<Image *> *kActualNearestNeighbors;
    vector<Image *> *rangeNeighbors;

    //open output file in write mode
    ofstream *outputFile = new ofstream(outputFilePath->c_str(), ios::out | ios::app);

    //calculating offset for query images' id
    offset = images->size();

    //For every query Image
    for(const auto& queryImage: *queryImages){
        *outputFile << "Query: " << queryImage->id - offset + 1 << endl;

        //hasherBoolean = true, means LSH

        // Get starting timepoint 
        auto start1 = high_resolution_clock::now();
        
        //get the results from kNearestNeighbors
        if(hasherBoolean){
            kNearestNeighbors = ((LSH *) myHasher)->kNearestNeighbors(queryImage, N);
        }else{ //else, HyperCube
            kNearestNeighbors = ((HyperCube *) myHasher)->kNearestNeighbors(queryImage, N, numberOfPoints, maxVertProbed);
        }

        // Get ending timepoint 
        auto stop1 = high_resolution_clock::now(); 


        int neighborCount = kNearestNeighbors->size();

        if(neighborCount == 0){
            delete kNearestNeighbors;
            continue;
        }


        // Get starting timepoint 
        auto start2 = high_resolution_clock::now();

        
        //get the results from exhaustiveNearestNeighbors only for as many as we found
        kActualNearestNeighbors = exhaustiveNearestNeighbors(queryImage, images, neighborCount);

        // Get ending timepoint 
        auto stop2 = high_resolution_clock::now();
    
        //Calculate durations
        auto durationLSH = duration_cast<microseconds>(stop1 - start1);
        auto durationTrue = duration_cast<microseconds>(stop2 - start2);

        //for each neighbor
        for(index = 0; index < neighborCount; index++){
            *outputFile << "Nearest neighbor-" << index + 1 << ": " << kNearestNeighbors->at(index)->id << endl;

            if(hasherBoolean){
                *outputFile << "distanceLSH: " << DIST_FUNC(queryImage, kNearestNeighbors->at(index)) << endl;
            }else{ //else, HyperCube
                *outputFile << "distanceHyperCube: " << DIST_FUNC(queryImage, kNearestNeighbors->at(index)) << endl;
            }
            
            *outputFile << "distanceTrue: " << DIST_FUNC(queryImage, kActualNearestNeighbors->at(index)) << endl;            
        }

        //we don't need the results anymore
        delete kNearestNeighbors;
        delete kActualNearestNeighbors;

        if(hasherBoolean){
            *outputFile << "tLSH: " << durationLSH.count() * 1.0 * pow(10, -6) << endl;
        }else{ //else, HyperCube
            *outputFile << "tHyperCube: " << durationLSH.count() * 1.0 * pow(10, -6) << endl;
        }

        *outputFile << "tTrue: " << durationTrue.count() * 1.0 * pow(10, -6) << endl; 
        
        //get the results from rangeSearch
        if(hasherBoolean){
            rangeNeighbors = ((LSH *) myHasher)->rangeSearch(queryImage, R);
        }else{ //else, HyperCube
            rangeNeighbors = ((HyperCube *) myHasher)->rangeSearch(queryImage, R, numberOfPoints, maxVertProbed);
        }

        *outputFile << "R-near neighbors:" << endl;
        //for each neighbor
        for(const auto& rangeNeighbor: *rangeNeighbors){
            *outputFile << rangeNeighbor->id << endl;
        }


        delete rangeNeighbors;
    }


    outputFile->close();
    delete outputFile;
}


//this returns [s1, s2, s3, ... ,sk, sTotal]
vector<double> *calculateSilhouettes(vector<Image *> *centroids, vector<set<Image *> *> *clusters){
    //count clusters
    unsigned int clusterCount = centroids->size();
    //we keep a vector where we will store the return values
    vector<double> *silhouettes = new vector<double>(clusterCount + 1, 0.0);

    unsigned int clusterIndex;

    //we have a memoization distance matrix to avoid recalculation of distances
    //we need to have this array be 0 + 1 + 2 + 3 + ... maxClusterSize - 1 =  in order to recalculate it
    unsigned int maxClusterSize = 0;
    for(set<Image *> *cluster : *clusters){
        if(cluster->size() > maxClusterSize){
            maxClusterSize = cluster->size();
        }
    }

    // make the array contiguous so as to be able to use one memset for instant clearing
    int *distanceMatrix = (int*) malloc(maxClusterSize * maxClusterSize * sizeof(int));
    
    double currClusterSilhouette;
    double sTotal = 0.0;

    clusterIndex = 0;
    //we have to calculate stuff for every cluster
    for(set<Image *> *cluster : *clusters){

        currClusterSilhouette = calculateClusterSilhouette(cluster, clusterIndex, clusters, centroids, distanceMatrix, maxClusterSize);

        //log to vector
        silhouettes->at(clusterIndex) = currClusterSilhouette;

        //add to sTotal
        sTotal += currClusterSilhouette / ( (double) clusterCount);

        clusterIndex ++;
    }

    //finally add sTotal to result
    silhouettes->at(clusterIndex) = sTotal;

    free(distanceMatrix);

    return silhouettes;

}

//this takes an image and a cluster and calculates the average distance of the image to 
//the cluster's images
//uses a memoization matrix for better performance, only if it is given
double calculateMeanDistanceToCluster(Image *imageToEvaluate, int imageToEvaluateIndex, set<Image *> *cluster, unsigned int maxClusterSize, int *distanceMatrix){
    double distanceMean = 0.0;
    unsigned int currentDistance;
    double denominator;
    unsigned int otherImageIndex;

    //the case where we check against a different cluster, so don't memoize
    if(distanceMatrix == NULL){
        //the number of "neighbors" is just clusterSize
        denominator = (double) cluster->size();

        //just calculate and sum the distances to other images
        for(Image *otherImage : *cluster){
            currentDistance = DIST_FUNC(imageToEvaluate, otherImage);
            distanceMean += ((double) currentDistance) / denominator;
        }
    } else { //in this case we will utilize the distance matrix

        //neighbors are cluster size - 1
        denominator = (double) (cluster->size() - 1);

        if(denominator == 0.0){
            denominator ++;
        }

        //we need this index for access to the distance matrix cells
        otherImageIndex = 0;

        for(Image *otherImage : *cluster){
            //check if we have the distance already
            if(distanceMatrix[imageToEvaluateIndex * maxClusterSize + otherImageIndex] != - 1){
                //get memoized distance
                currentDistance = distanceMatrix[imageToEvaluateIndex * maxClusterSize + otherImageIndex];
            } else {
                //calculate and save distance
                currentDistance = DIST_FUNC(imageToEvaluate, otherImage);
                distanceMatrix[imageToEvaluateIndex * maxClusterSize + otherImageIndex] = currentDistance;
            }

            distanceMean += ( (double) currentDistance) / denominator;
     
        }

    }

    //now we have the distanceMean, just return it
    return distanceMean;
}


int getNextBestCluster(Image *imageToEvaluate, vector<Image *> *centroids, int currentClusterIndex){
    unsigned int bestClusterDistance, candidateClusterDistance;
    int currOtherClusterIndex, clusterCount, nextBestClusterIndex;

    clusterCount = centroids->size();

    //we must find next best cluster based on distance to centroid
    bestClusterDistance = 0;
    for(currOtherClusterIndex = 0; currOtherClusterIndex < clusterCount; currOtherClusterIndex ++){
        //ignore current cluster obviously
        if(currOtherClusterIndex == currentClusterIndex){
            continue;
        }

        //minimize distance
        if(bestClusterDistance == 0){
            bestClusterDistance = DIST_FUNC(imageToEvaluate, centroids->at(currOtherClusterIndex));
            nextBestClusterIndex = currOtherClusterIndex;
        } else {
            //if distance is lower, update
            candidateClusterDistance = DIST_FUNC(imageToEvaluate, centroids->at(currOtherClusterIndex));
            if(candidateClusterDistance < bestClusterDistance){
                bestClusterDistance = candidateClusterDistance;
                nextBestClusterIndex = currOtherClusterIndex;
            }
        }
    }

    return nextBestClusterIndex;
}


double calculateClusterSilhouette(set<Image *> *cluster, unsigned int clusterIndex, vector<set<Image *> *> *clusters, vector<Image *> *centroids, int *distanceMatrix, unsigned int maxClusterSize){
    double currClusterSilhouette;
    double currA, currB;
    int nextBestClusterIndex;
    int imageIndex;


    if(cluster->size() == 0){
        return 0.0;
    }

    currClusterSilhouette = 0.0;

    //we must initialise the distance array
    memset(distanceMatrix, -1, sizeof(distanceMatrix[0]) * maxClusterSize * maxClusterSize);
    
    imageIndex = 0;
    //we have to calculate the average distance for every image in the cluster
    for(Image *imageToEvaluate : *cluster){
        
        //get ai 
        currA = calculateMeanDistanceToCluster(imageToEvaluate, imageIndex, cluster, maxClusterSize, distanceMatrix);

        //we must find next best cluster
        nextBestClusterIndex = getNextBestCluster(imageToEvaluate, centroids, clusterIndex);

        //get bi
        currB = calculateMeanDistanceToCluster(imageToEvaluate, imageIndex, clusters->at(nextBestClusterIndex), maxClusterSize, NULL);
        

        //we can now use the mathematical formula to add to si
        if(currA < currB){
            currClusterSilhouette += (1.0 - currA / currB) / ( (double) cluster->size());
        } else {
            currClusterSilhouette += (currB / currA - 1.0) / ( (double) cluster->size());
        }

        imageIndex ++;
    }

    return currClusterSilhouette;
}

bool isPrefix(string *s1, string *s2) {
    const char*p = s1->c_str();
    const char*q = s2->c_str();
    while (*p&&*q)
        if (*p++!=*q++)
            return false;
    return true;
}

void readConfigFile(string *configFilePath, int *numberOfClusters, int *L, int *kLSH, int *MHyperCube, int *kHyperCube, int *probesHyperCube){ 
    //open the config file
    ifstream *configFile = new ifstream(configFilePath->c_str(), ios::in);

    string *currLine = new string();

    string *dummy = new string();

    //insert default values beforehand
    *L = 3;
    *kLSH = 4;
    *MHyperCube = 10;
    *kHyperCube = 3;
    *probesHyperCube = 2;
    
    //number of clusters MUST be initialised, insert 0 in order to check
    *numberOfClusters = 0;

    size_t last_index;

    //specified option prefixes
    string *options[6] = {
                          new string("number_of_clusters"),
                          new string("number_of_vector_hash_tables"),
                          new string("number_of_vector_hash_functions"),
                          new string("max_number_M_hypercube"),
                          new string("number_of_hypercube_dimensions"),
                          new string("number_of_probes")
                          };

    int *optionValues[6] = { 
                            numberOfClusters,
                            L,
                            kLSH,
                            MHyperCube,
                            kHyperCube,
                            probesHyperCube
                            };

    int optionIndex;

    //loop through all lines in file
    while(getline(*configFile, *currLine)){
        //check all possible options
        for(optionIndex = 0; optionIndex < 6; optionIndex ++){
            if(isPrefix(options[optionIndex], currLine)){
                //we must update the value from the file
                stringstream os(*currLine);

                //find where the number starts
                last_index = currLine->find_last_not_of("0123456789");

                //drop all text before number and parse to variable
                *(optionValues[optionIndex]) = stoi(currLine->substr(last_index + 1));
            }
        }
    }

    if(*numberOfClusters == 0){
        //if we didn't get the number of clusters, ask user
        cout << "No number of clusters specified in config file. Please enter one:" << endl;
        cin >> *numberOfClusters;
    }


    delete configFile;
    delete currLine;
    delete dummy;

    for(optionIndex = 0; optionIndex < 6; optionIndex ++){
        delete options[optionIndex];
    }

    return;

}





//this will run the specified clustering algorithm and produce results on the output file
void runClusteringAlgorithm(vector<Image *> *images, unsigned int numberOfClusters, string *method, ofstream *outputFile, LSH *myLSH, HyperCube *myHyperCube, int MHyperCube, int probesHyperCube, bool complete){
    int clusteringIterations = 0;
    int clusterIndex;
    const auto* sep = "";
    
    bool classic = false, LSH = false, HyperCube = false;

    classic = *method == "Classic";
    LSH = *method == "LSH";
    HyperCube = *method == "Hypercube";


    //now the clustering begins
    auto clusteringStart = high_resolution_clock::now();

    //get the first centroids
    vector<Image *> *centroids = initialCentroidImages(images, numberOfClusters);
    vector<set<Image *> *> *clusters = NULL;

    do{

        if(clusters != NULL){
            for(set<Image *> *cluster : *clusters){
                delete cluster;
            }
            delete clusters;
        }

        if(classic){
            clusters = LloydsAlgorithm(centroids, images);
        } else if(LSH){
            clusters = myLSH->clusterizeLSH(centroids, images);
        } else if(HyperCube){
            clusters = myHyperCube->clusterizeHyperCube(centroids, images, MHyperCube, probesHyperCube);
        } else {
            //no method given, abort
            cout << "No clustering method given, aborting!" << endl;
            return;
        }

        //we did one iteration
        clusteringIterations ++;

        //delete old centroids only if they are not the first ones
        if(clusteringIterations > 1){
            for(Image *centroid : *centroids){
                delete centroid;
            }
        }
        delete centroids;

        //now get the new centroids 
        centroids = centroidize(clusters);
        
        //if new centroids are pretty close to old ones, we should stop
        //TODO: implement above check

        //if we reached max iterations, stop
        if(clusteringIterations == MAX_CLUSTERING_ITERATIONS){
            break;
        }





    }while(true);


    //now clustering is over, calculate and output the time
    auto clusteringStop = high_resolution_clock::now();

    auto durationClustering = duration_cast<microseconds>(clusteringStop - clusteringStart);


    //after clustering, output to file 
    //for every cluster
    clusterIndex = 0;
    for(set<Image *> *cluster : *clusters){

        *outputFile << "CLUSTER-" << clusterIndex + 1 << " {size: " << cluster->size() << ", centroid: [";

        sep = "";
        for(unsigned char pixel : *centroids->at(clusterIndex)->imagePixels) {
            *outputFile << sep << (int) pixel;
            sep = ", ";
        }


        *outputFile << "]}" << endl;

        clusterIndex ++;
    }

    *outputFile << "clustering_time: " << durationClustering.count() * 1.0 * pow(10, -6) << endl; 

    vector<double> *silhouettes = calculateSilhouettes(centroids, clusters);

    *outputFile << "Silhouette: [";

    sep = "";
    for(double currSil : *silhouettes) {
        *outputFile << sep << currSil;
        sep = ", ";
    }

    *outputFile << "]" << endl;

    delete silhouettes;


    //if complete, print cluster info
    if(complete){

        //add padding
        *outputFile << endl << endl;

        clusterIndex = 0;
        for(set<Image *> *cluster : *clusters){

            *outputFile << "CLUSTER-" << clusterIndex + 1 << " {[";

            //output centroid
            sep = "";
            for(unsigned char pixel : *centroids->at(clusterIndex)->imagePixels) {
                *outputFile << sep << (int) pixel;
                sep = ", ";
            }

            *outputFile << "], ";


            sep = "";
            for(Image *currClusterImage : *cluster) {
                *outputFile << sep << currClusterImage->id;
                sep = ", ";
            }

            *outputFile << "}" << endl;

            clusterIndex ++;
        }
    }


    //add some padding
    *outputFile << endl << endl;


    for(Image *centroid : *centroids){
        delete centroid;
    }
    delete centroids;

    for(set<Image *> *cluster : *clusters){
        delete cluster;
    }
    delete clusters;

    return;
}