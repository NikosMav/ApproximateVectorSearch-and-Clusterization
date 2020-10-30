#include <vector>
#include "../headers/Image.h"

Image::Image(int id, vector<unsigned char> * imagePixels){
    this->id = id;
    this->imagePixels = imagePixels;
}

Image::~Image(){
    delete this->imagePixels;
}