#ifndef LSH_H
#define LSH_H

#include <vector>
#include "../headers/Hashtable.h"
#include "../headers/Aux.h"
#include <queue>
#include <set>

class LSH
{
public:
    unsigned int L; // number of different hashtables that LSH uses
    vector<Hashtable *> *hashtables; // vector of different hashtables used for LSH technique


    LSH(unsigned int L, int numberOfBuckets);
    ~LSH();
    void insert(Image *image);
    Image *NearestNeighbor(Image *queryImage);
    vector<Image *> *kNearestNeighbors(Image *queryImage, unsigned int k);
    vector<Image *> *rangeSearch(Image *queryImage, double R, set<unsigned int> *ignoredPoints = NULL);
    vector<set <Image *> *> *clusterizeLSH(vector<Image *> *centroids, vector<Image *> *inputImages);
};

#endif