#include <vector>
#include <random>

#include "../headers/Hashtable.h"
#include "../headers/Image.h"
#include "../headers/Bucket.h"


using namespace std;

Hashtable::Hashtable(int numberOfBuckets){
    //make random generator 
    random_device rd;
    mt19937 e2(rd());
    uniform_real_distribution<> dist(0, this->sCeil);

    vector<float> *currentDisturbanceVector = NULL;

    this->numberOfBuckets = numberOfBuckets;

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
            currentDisturbanceVector->at(currentDimension) = dist(e2);
        }

        //assign the created vector to the parent vector
        this->disturbanceVectors->at(disturbanceVectorCount) = currentDisturbanceVector;

    }

    this->bucketArray = new vector<Bucket *>(numberOfBuckets, NULL);

    int bucketNum;
    for(bucketNum = 0; bucketNum < numberOfBuckets; bucketNum++){
        this->bucketArray->at(bucketNum) = new Bucket();
    }


}

Hashtable::~Hashtable(){
    for(vector<float> *disturbanceVector : *this->disturbanceVectors){
        delete disturbanceVector;
    }
    delete this->disturbanceVectors;

    for(Bucket *currBucket : *this->bucketArray){
        delete currBucket;
    }

    delete this->bucketArray;

}


unsigned int Hashtable::hashFunction(Image *image){
    

    unsigned int currH;
    unsigned int hConcatenation = 0;
    unsigned int shiftVal = 32 / this->k;

    //for every disturbance vector
    for(const auto& disturbanceVector: *(this->disturbanceVectors)) {
        currH = this->calculateH(image, disturbanceVector);
        //shift concatenation left by 32/k bits
        hConcatenation <<= shiftVal;
        //add the bits of current h function to end of concatenation
        hConcatenation |= currH;
    }


    return hConcatenation % this->numberOfBuckets;
}

unsigned int Hashtable::calculateH(Image *image, vector<float> *const disturbanceVector){
    
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

    //make the mask instead of %
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

void Hashtable::insert(Image* image){
    unsigned int goingToBucket;
    goingToBucket = this->hashFunction(image);

    this->bucketArray->at(goingToBucket)->bucketImages->push_back(image);
}


//returns the bunch of buckets that the image hashes to
vector<Image *> *Hashtable::getNeighbors(Image *image){
    unsigned int goingToBucket;
    goingToBucket = this->hashFunction(image);

    //return the raw vector which has the images
    return this->bucketArray->at(goingToBucket)->bucketImages;
}








float Hashtable::sCeil;
unsigned int Hashtable::k;
unsigned int Hashtable::hModulo;
unsigned int Hashtable::d;
unsigned int Hashtable::m;
unsigned int *Hashtable::mPowerModulos;