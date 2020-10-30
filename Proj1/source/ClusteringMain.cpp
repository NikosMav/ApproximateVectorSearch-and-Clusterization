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
    string *inputFilePath, *configFilePath, *outputFilePath, *methodStr;
    bool complete;
    unsigned int magic_number, numberOfImages, pixels, k, N, numberOfQueries;
    int numberOfClusters, kLSH, MHyperCube, kHyperCube, probesHyperCube, L;
    double R;
    vector<Image *> *images, *queryImages;

    LSH *myLSH = NULL;
    HyperCube *myHyperCube = NULL;

    inputFilePath = getOptionValue("-i", argc, argv);
    if(inputFilePath == NULL){
        cout << "Please enter the input file path:" << endl;

        inputFilePath = new string();
        getline(cin, *inputFilePath);
    }

    configFilePath = getOptionValue("-c", argc, argv);
    if(configFilePath == NULL){
        cout << "Please enter the config file path:" << endl;

        configFilePath = new string();
        getline(cin, *configFilePath);
    }

    outputFilePath = getOptionValue("-o", argc, argv);
    if(outputFilePath == NULL){
        cout << "Please enter the output file path:" << endl;

        outputFilePath = new string();
        getline(cin, *outputFilePath);
    }

    complete = getOption("-complete", argc, argv);

    //get method, if it is NULL, do all
    methodStr = getOptionValue("-m", argc, argv);

    //get info from config file
    readConfigFile(configFilePath, &numberOfClusters, &L, &kLSH, &MHyperCube, &kHyperCube, &probesHyperCube);

    delete configFilePath;

    //read the input file
    fstream *trainFile;
    images = openTrainFile(inputFilePath, pixels, numberOfImages, &trainFile);

    delete inputFilePath;

    trainFile->close();

    delete trainFile;

    //now open the output file in order to write to it
    ofstream *outputFile = new ofstream(outputFilePath->c_str(), ios::out);
    delete outputFilePath;

    //see what methdos we need to run
    if(methodStr == NULL){ // all algorithms

        //actually assign memory for the method
        methodStr = new string();

        //set the method
        methodStr->assign("Classic");

        //first output for classic
        *outputFile << "Algorithm: Lloyds" << endl;

        runClusteringAlgorithm(images, numberOfClusters, methodStr, outputFile, myLSH, myHyperCube, MHyperCube, probesHyperCube, complete);

        //build LSH database
        myLSH = initializeLSH(kLSH, pixels, L, images);


        //set the method
        methodStr->assign("LSH");

        //output
        *outputFile << "Algorithm: Range Search LSH" << endl;

        runClusteringAlgorithm(images, numberOfClusters, methodStr, outputFile, myLSH, myHyperCube, MHyperCube, probesHyperCube, complete);

        delete myLSH;

        free(Hashtable::mPowerModulos);

        //build HyperCube database
        myHyperCube = initializeHyperCube(kHyperCube, pixels, images);


        //set the method
        methodStr->assign("Hypercube");

        //output
        *outputFile << "Algorithm: Range Search Hypercube" << endl;

        runClusteringAlgorithm(images, numberOfClusters, methodStr, outputFile, myLSH, myHyperCube, MHyperCube, probesHyperCube, complete);

        delete myHyperCube;

    } else if(*methodStr == "LSH"){ //lsh only
        //build LSH database
        myLSH = initializeLSH(kLSH, pixels, L, images);

        //output
        *outputFile << "Algorithm: Range Search LSH" << endl;

        runClusteringAlgorithm(images, numberOfClusters, methodStr, outputFile, myLSH, myHyperCube, MHyperCube, probesHyperCube, complete);

        delete myLSH;
    } else if(*methodStr == "Hypercube"){ //hypercube only
        //build HyperCube database
        myHyperCube = initializeHyperCube(kHyperCube, pixels, images);

        //output
        *outputFile << "Algorithm: Range Search Hypercube" << endl;

        runClusteringAlgorithm(images, numberOfClusters, methodStr, outputFile, myLSH, myHyperCube, MHyperCube, probesHyperCube, complete);

        delete myHyperCube;

    } else { //classic only
        *outputFile << "Algorithm: Lloyds" << endl;

        runClusteringAlgorithm(images, numberOfClusters, methodStr, outputFile, myLSH, myHyperCube, MHyperCube, probesHyperCube, complete);
    }


    if(methodStr != NULL){
        delete methodStr;
    }

    outputFile->close();
    delete outputFile;


    for(Image *image : *images){
        delete image;
    }
    delete images;

    free(Hashtable::mPowerModulos);

    return 0;

}