
#ifndef PDFANALYSIS_H
#define PDFANALYSIS_H
#define MAXWEIGHT 110

#include "core/Enums.h"
#include "../../../../Smurf/Core/SmurfTree.h"
#include "core/SmurfSample.h"

// C++ includes
#include <iostream>
#include <vector>

static const unsigned int set_ = 2;
static const unsigned int genset_ = 1;

class MyScanChain {
    public:
        MyScanChain(float analysis, Option option, std::string pdfName, unsigned int pdfSubset);
        ~MyScanChain() {};
        int ScanChain(SmurfSample *sample, std::string pdfName);

    private:
        bool cuts(SmurfTree *tree, DataType dataType, const unsigned int jetbin, const float& dymva);
        float analysis_;
        Option option_; 
};

#endif

