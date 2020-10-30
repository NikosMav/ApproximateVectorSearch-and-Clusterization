#ifndef IMAGE
#define IMAGE

#include <vector>

using namespace std;

class Image
{
public:
    int id;
    vector<unsigned char> * imagePixels;

    Image(int id, vector<unsigned char> * imagePixels);
    ~Image();
};


#endif