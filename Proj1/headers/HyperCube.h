#ifndef HYPERCUBE
#define HYPERCUBE

#include <vector>
#include <set>
#include "../headers/Bucket.h"
#include "../headers/Image.h"
#include "../headers/Aux.h"


using namespace std;

class HyperCube
{
private:
    unsigned int k;
    unsigned int d;
    unsigned int sCeil;
    unsigned int hModulo;
    vector<Bucket *> *bucketArray;
    vector<vector <unsigned char> *> *projectionFunctions;
    vector<vector<float> *> *disturbanceVectors; 
public:
    HyperCube(unsigned int k, unsigned int hModulo, unsigned int d);
    ~HyperCube();

    void insert(Image *);
    unsigned int hashFunction(Image *image);
    unsigned int calculateH(Image *image, vector<float> *const disturbanceVector);
    static unsigned int *mPowerModulos; // this keeps the values m^(i) mod M, which we reuse again and again
    Image *NearestNeighbor(Image *queryImage, unsigned int numberOfPoints, unsigned int maxVertProbed);
    vector<Image *> *kNearestNeighbors(Image *queryImage, unsigned int k, unsigned int numberOfPoints, unsigned int maxVertProbed);
    vector<Image *> *rangeSearch(Image *queryImage, double R, unsigned int numberOfPoints, unsigned int maxVertProbed,  set<unsigned int> *ignoredPoints = NULL);
    vector<unsigned int> *hammingNeighbors(unsigned int node);
    vector<set <Image *> *> *clusterizeHyperCube(vector<Image *> *centroids, vector<Image *> *inputImages, unsigned int numberOfPoints, unsigned int maxVertProbed);
};
#endif