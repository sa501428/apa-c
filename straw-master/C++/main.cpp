/*
  The MIT License (MIT)

  Copyright (c) 2011-2016 Broad Institute, Aiden Lab

 Permission is hereby granted, free of charge, to any person obtaining a copy
 of this software and associated documentation files (the "Software"), to deal
 in the Software without restriction, including without limitation the rights
 to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 copies of the Software, and to permit persons to whom the Software is
 furnished to do so, subject to the following conditions:

 The above copyright notice and this permission notice shall be included in
 all copies or substantial portions of the Software.

 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 THE SOFTWARE.
*/
#include <iostream>
#include <string>
#include <algorithm>  // for transform
#include <cctype>    // for tolower
#include "straw.h"
#include "hic_slice.h"
using namespace std;

// Helper function to check argument strings
bool containsIgnoreCase(const string& str, const string& substr) {
    string strLower = str;
    string substrLower = substr;
    transform(strLower.begin(), strLower.end(), strLower.begin(), ::tolower);
    transform(substrLower.begin(), substrLower.end(), substrLower.begin(), ::tolower);
    return strLower.find(substrLower) != string::npos;
}

int main(int argc, char *argv[])
{
    // Check if this is a dump command
    if (argc > 1 && string(argv[1]) == "dump") {
        if (argc < 9 || argc > 10) {
            cerr << "Incorrect arguments for dump command" << endl;
            cerr << "Usage: straw dump <observed/oe/expected> <NONE/VC/VC_SQRT/KR> <hicFile> <BP/FRAG> <binsize> <outputFile> <compressed> [-intra-short|-intra-long|-inter|-intra]" << endl;
            exit(1);
        }
        string matrixType = argv[2];
        string norm = argv[3];
        string fname = argv[4];
        string unit = argv[5];
        int32_t binsize = stoi(argv[6]);
        string outputPath = argv[7];
        bool compressed = (string(argv[8]) == "1" || string(argv[8]) == "true" || string(argv[8]) == "compressed");
        ContactFilter filter = ContactFilter::ALL;
        
        // Parse optional filter argument more flexibly
        if (argc == 10) {
            string arg = argv[9];
            if (containsIgnoreCase(arg, "inter")) {
                filter = ContactFilter::INTER;
            } else if (containsIgnoreCase(arg, "intra")) {
                if (containsIgnoreCase(arg, "short")) {
                    filter = ContactFilter::INTRA_SHORT;
                } else if (containsIgnoreCase(arg, "long")) {
                    filter = ContactFilter::INTRA_LONG;
                } else {
                    filter = ContactFilter::INTRA;
                }
            }
        }

        dumpGenomeWideDataAtResolution(matrixType, norm, fname, unit, binsize, outputPath, compressed, filter);
        return 0;
    }

    // Original functionality
    if (argc != 7 && argc != 8) {
        cerr << "Incorrect arguments" << endl;
        cerr << "Usage: straw [observed/oe/expected] <NONE/VC/VC_SQRT/KR> <hicFile(s)> <chr1>[:x1:x2] <chr2>[:y1:y2] <BP/FRAG/MATRIX> <binsize>" << endl;
        cerr << "   or: straw dump <observed/oe/expected> <NONE/VC/VC_SQRT/KR> <hicFile> <BP/FRAG> <binsize> <outputFile>" << endl;
        exit(1);
    }
    int offset = 0;
    string matrixType = "observed";
    if(argc == 8){
        offset = 1;
        matrixType = argv[1];
    }
    string norm = argv[1 + offset];
    string fname = argv[2 + offset];
    string chr1loc = argv[3 + offset];
    string chr2loc = argv[4 + offset];
    string unit = argv[5 + offset];
    string size = argv[6 + offset];
    int32_t binsize = stoi(size);

    if(unit == "MATRIX"){
        vector<vector<float>> matrix = strawAsMatrix(matrixType, norm, fname, chr1loc, chr2loc, "BP", binsize);
        for(int i = 0; i < matrix.size(); i++){
            for(int j = 0; j < matrix[i].size(); j++){
                cout << matrix[i][j] << "\t";
            }
            cout << endl;
        }
    } else {
        vector<contactRecord> records;
        records = straw(matrixType, norm, fname, chr1loc, chr2loc, unit, binsize);
        for (int i = 0; i < records.size(); i++) {
            printf("%d\t%d\t%.14g\n", records[i].binX, records[i].binY, records[i].counts);
        }
    }
}
