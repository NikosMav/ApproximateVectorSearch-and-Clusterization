#include <vector>
#include "../headers/Bucket.h"

Bucket::Bucket(){
    this->bucketImages = new vector<Image *>();
}

Bucket::~Bucket(){
    delete this->bucketImages;
}