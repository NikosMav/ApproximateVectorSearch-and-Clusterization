#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <arpa/inet.h>
#include "../headers/Aux.h"
#include "../headers/Hashtable.h"
#include "../headers/Image.h"
#include "../headers/LSH.h"
#include "../headers/HyperCube.h"
#include <cmath>

using namespace std;

int main(int argc, char **argv){
    string *inputFilePath, *kStr, *NStr, *RStr, *queryFilePath, *outputFilePath, *MStr, *probesStr;
    unsigned int numberOfImages, pixels, k, N, numberOfQueries, numberOfPoints, maxVertProbed;
    double R;
    vector<Image *> *images, *queryImages;

    getCommonArguments(argc, argv, &inputFilePath, &queryFilePath, &outputFilePath, NStr, RStr, N, R);

    MStr = getOptionValue("-M", argc, argv);
    if(MStr == NULL){
        numberOfPoints = 10;
    } else{
        numberOfPoints = stoi(*MStr);
        delete MStr;
    }

    probesStr = getOptionValue("-probes", argc, argv);
    if(probesStr == NULL){
        maxVertProbed = 2;
    } else{
        maxVertProbed = stoi(*probesStr);
        delete probesStr;
    }

    kStr = getOptionValue("-k", argc, argv);
    if(kStr == NULL){
        k = 14;
    } else{
        k = stoi(*kStr);
        delete kStr;
    }

    fstream *trainFile;
    images = openTrainFile(inputFilePath, pixels, numberOfImages, &trainFile);
    delete inputFilePath;

    HyperCube* myHyperCube = initializeHyperCube(k, pixels, images);
    bool hasherBoolean;
    string answer;

    do{
        queryImages = openQueryFile(queryFilePath, numberOfQueries, pixels);
        delete queryFilePath;

        hasherBoolean = false;
        //Output function
        calculateOutput(outputFilePath, pixels, numberOfImages, images, queryImages, myHyperCube, hasherBoolean, k, N, R, numberOfPoints, maxVertProbed);

        //clean query images    
        for(Image *image : *queryImages){
            delete image;
        }
        delete queryImages;

        cout << "Repeat process with different query file? (Y/N)" << endl;
        getline(cin, answer);

        if(answer == "Y"){
            cout << "Please enter new query file name: " << endl;
            queryFilePath = new string();
            getline(cin, *queryFilePath);

        }else{
            break;
        }

    }while(true);

    delete outputFilePath;

    delete myHyperCube;

    //clean actual images    
    for(Image *image : *images){
        delete image;
    }
    delete images;

    trainFile->close();

    delete trainFile;

    //clear mPowerModulos
    free(HyperCube::mPowerModulos);
}