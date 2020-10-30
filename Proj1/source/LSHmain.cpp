#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <arpa/inet.h>
#include "../headers/Aux.h"
#include "../headers/Hashtable.h"
#include "../headers/Image.h"
#include "../headers/LSH.h"
#include <cmath>

using namespace std;

int main(int argc, char **argv){
    string *inputFilePath, *kStr, *LStr, *NStr, *RStr, *queryFilePath, *outputFilePath;
    unsigned int numberOfImages, pixels, k, L, N, numberOfQueries;
    double R;
    vector<Image *> *images, *queryImages;

    getCommonArguments(argc, argv, &inputFilePath, &queryFilePath, &outputFilePath, NStr, RStr, N, R);

    LStr = getOptionValue("-L", argc, argv);
    if(LStr == NULL){
        L = 5;
    } else{
        L = stoi(*LStr);
        delete LStr;
    }

    kStr = getOptionValue("-k", argc, argv);
    if(kStr == NULL){
        k = 4;
    } else{
        k = stoi(*kStr);
        delete kStr;
    }

    fstream *trainFile;
    images = openTrainFile(inputFilePath, pixels, numberOfImages, &trainFile);
    delete inputFilePath;

    LSH* myLSH = initializeLSH(k, pixels, L, images);
     bool hasherBoolean;
    string answer;

    do{
        queryImages = openQueryFile(queryFilePath, numberOfQueries, pixels);
        delete queryFilePath;

        hasherBoolean = true;
        //Output function
        calculateOutput(outputFilePath, pixels, numberOfImages, images, queryImages, myLSH, hasherBoolean, k, N, R);

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

    delete myLSH;

    //clean actual images    
    for(Image *image : *images){
        delete image;
    }
    delete images;

    trainFile->close();

    delete trainFile;

    //clear mPowerModulos
    free(Hashtable::mPowerModulos);
}