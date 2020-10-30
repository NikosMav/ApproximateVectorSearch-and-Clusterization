#ifndef AUX
#define AUX

#define MAX_CLUSTERING_ITERATIONS 10

using namespace std;

#include <string>
#include <vector>
#include <set>
#include <random>
#include "../headers/Image.h"
#include "../headers/LSH.h"

class LSH;
class HyperCube;

string *getOptionValue(string option, int argc, char **argv);
bool getOption(string option, int argc, char **argv);
unsigned int parseInteger(fstream *binaryFile);
Image *parseImage(fstream *binaryFile, unsigned int pixels);
vector<Image *> *parseImages(fstream *binaryFile, unsigned int numberOfImages, unsigned int pixels);
void initializeModulosForHashTable();
int manhattanDistance(Image *image1, Image *image2);

struct candidateImage{
    Image *candidateImage;
    unsigned int score;
};
typedef struct candidateImage candidateImage;

struct imageClusterData{
    int assignedCluster;
    bool considering;
    bool set;
    int runningMinDist;
};
typedef struct imageClusterData imageClusterData;

void setCandidateImage(candidateImage* candidate, Image* currCandidateImage, Image* queryImage);
vector<Image *> *initialCentroidImages(vector<Image *> *images, unsigned int k);
vector<set<Image *> *> * LloydsAlgorithm(vector<Image *> *centroids, vector<Image *> *images);
vector<Image *> *centroidize(vector<set<Image *> *> *clusters);
vector<Image *> *exhaustiveNearestNeighbors(Image *queryImage, vector<Image *> *images, unsigned int k);

void getCommonArguments(int argc, char **argv, string **inputFilePath, string **queryFilePath, string **outputFilePath, string *NStr, string *RStr, unsigned int &N, double &R);
vector<Image *> *openTrainFile(string *inputFilePath, unsigned int &pixels, unsigned int &numberOfImages, fstream **trainFile);
vector<Image *> *openQueryFile(string *queryFilePath, unsigned int &numberOfQueries, unsigned int pixels);
void calculateOutput(string *outputFilePath, unsigned int pixels, unsigned int numberOfImages, vector<Image *> *images, vector<Image *> *queryImages, void* myHasher, bool hasherBoolean, unsigned int k, unsigned int N, double R, unsigned int numberOfPoints = 0, unsigned int maxVertProbed = 0);
LSH* initializeLSH(unsigned int k, unsigned int pixels, unsigned int L, vector<Image *> *images);
HyperCube* initializeHyperCube(unsigned int k, unsigned int pixels, vector<Image *> *images);

double calculateMeanDistanceToCluster(Image *imageToEvaluate, int imageToEvaluateIndex, set<Image *> *cluster, unsigned int maxClusterSize, int *distanceMatrix);
vector<double> *calculateSilhouettes(vector<Image *> *centroids, vector<set<Image *> *> *clusters);
int getNextBestCluster(Image *imageToEvaluate, vector<Image *> *centroids, int currentClusterIndex);
double calculateClusterSilhouette(set<Image *> *cluster, unsigned int clusterIndex, vector<set<Image *> *> *clusters, vector<Image *> *centroids, int *distanceMatrix, unsigned int maxClusterSize);

void updateMinDistancesToCentroids(vector<Image *> *images, bool *isChosen, vector<Image *> *chosenImages, unsigned int currCentroidNumber, vector<unsigned int> *minDistances, unsigned int *maxDistance);

void calculatePartialSums(vector<double> *partialSums, unsigned int imagesCount, bool *isChosen, vector<unsigned int> *minDistances, unsigned int maxDistance);

void chooseNextCentroid(mt19937 &e2, vector<double> *partialSums, vector<Image *> *images, vector<Image *> *chosenImages, unsigned int currCentroidNumber, bool *isChosen);

bool isPrefix(string *s1, string *s2);

void readConfigFile(string *configFilePath, int *numberOfClusters, int *L, int *kLSH, int *MHyperCube, int *kHyperCube, int *probesHyperCube);

void runClusteringAlgorithm(vector<Image *> *images, unsigned int numberOfClusters, string *method, ofstream *outputFile, LSH *myLSH, HyperCube *myHyperCube, int MHyperCube, int probesHyperCube, bool complete);
#endif