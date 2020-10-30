#ifndef HASHTABLE
#define HASHTABLE
#define W 40
#define BUCKET_FACTOR 128

#include <vector>
#include "../headers/Image.h"
#include "../headers/Bucket.h"

using namespace std;

class Hashtable
{
private:
    vector<vector<float> *> *disturbanceVectors; 

public:
    int numberOfBuckets;
    vector<Bucket *> *bucketArray;

    static float sCeil; //specifies the upper limit of s_i for all disturbance vectors s
    static unsigned int k; //specifies the number of disturbance vectors s, also hi and M (it is k)
    static unsigned int hModulo; //this is M, the number we modulo after calculating the polynomial of H
    static unsigned int d; //this is d, the dimension of the images
    static unsigned int m;
    static unsigned int *mPowerModulos; // this keeps the values m^(i) mod M, which we reuse again and again

    Hashtable(int numberOfBuckets);
    ~Hashtable();

    unsigned int hashFunction(Image *image);
    unsigned int calculateH(Image* image, vector<float> *const disturbanceVector);
    void insert(Image* image);
    vector<Image *> *getNeighbors(Image *image);
};

#endif