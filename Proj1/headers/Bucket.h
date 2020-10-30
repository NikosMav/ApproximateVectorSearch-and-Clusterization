#ifndef BUCKET
#define BUCKET

#include "../headers/Image.h"

class Bucket {
public:
    vector<Image *> *bucketImages;

    Bucket();
    ~Bucket();
};

#endif