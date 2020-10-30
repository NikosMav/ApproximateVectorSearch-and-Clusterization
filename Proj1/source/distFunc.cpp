#include "../headers/Image.h"

int manhattanDistance(Image *image1, Image *image2) {
    int dimension, i;
    dimension = image1->imagePixels->size();

    int sum = 0;

    unsigned char image1Value, image2Value;
    for(i = 0; i < dimension; i++){
        image1Value = image1->imagePixels->at(i);
        image2Value = image2->imagePixels->at(i);

        if(image1Value > image2Value) {
            sum += image1Value - image2Value;
        }else{
            sum += image2Value - image1Value;
        }
    }
    
    return sum;
}