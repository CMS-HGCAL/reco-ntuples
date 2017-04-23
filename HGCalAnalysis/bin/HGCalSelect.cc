// c++ headers
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cstdlib>

// ROOT headers
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TClonesArray.h>
#include <TChain.h>
#include <TFileCollection.h>
#include <TSelector.h>
#include <TStopwatch.h>

// HGCAL headers
#include "RecoNtuples/HGCalAnalysis/interface/HGCalSelector.h"

using namespace std;

int main(int argc, char **argv) 
{
    int showTiming = 1;
    TStopwatch timer;

    gROOT->SetBatch();
    //gSystem->Load("libFWCoreFWLite.so");
    //AutoLibraryLoader::enable();

    // Get arguments
    if (argc != 3) {
        cout << "HGCalSelect <number of events> <input file(s)>" << endl;
        throw runtime_error("bad input");
    } 
    string sMaxEvents  = argv[1];
    string inputFiles  = argv[2];
    string option      = "";
    
    // Handle input files
    TChain* myChain = new TChain("Events");
    string inputFileExt = "root"; //get_file_extension(inputFiles);
    if (inputFileExt == "root") {
        if (myChain->AddFile(inputFiles.c_str())) {
            cout << "Trying to open " << inputFiles << "." << endl;
        } else {
            cout << "Failed to open " << inputFiles << endl;
            throw runtime_error("bad input");
        }
    } else if (inputFileExt == "txt") {
        TFileCollection fc;
        if (fc.AddFromFile(inputFiles.c_str()) && myChain->AddFileInfoList((TCollection*) fc.GetList()) ) {
            cout << "Trying to open " << inputFiles << "." << endl;
        } else {
            cout << "Failed to open " << inputFiles << endl;
            throw runtime_error("bad input");
        }
    } else {
        cout << "Unrecognized file extension: " << inputFileExt << endl;
        throw runtime_error("bad input");
    }

    cout << "Successfully openned input file; loading tree..." << endl;
    if (myChain->LoadTree(0) < 0) {  // needed because TChain =/= TTree
        cout << "Failed to load the entries in the file." << endl;
        throw runtime_error("bad input");
    }
    cout << "Successfully loaded tree..." << endl;

    if (showTiming) {
        timer.Start();
    }

    // Run
    long int maxEvents = stol(sMaxEvents);
    if (maxEvents < 0) 
        maxEvents = myChain->GetEntriesFast();

    HGCalSelector *selector = new HGCalSelector();
    myChain->Process(selector, option.c_str(), maxEvents, 0);

    // Done
    if (showTiming) {
        timer.Stop();
        cout << "CPU  time: " << timer.CpuTime() << endl;
        cout << "Real time: " << timer.RealTime() << endl;
    }
    return 0;
}
