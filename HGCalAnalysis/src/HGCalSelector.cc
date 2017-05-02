#define HGCalSelector_cxx

#include "RecoNtuples/HGCalAnalysis/interface/HGCalSelector.h"
#include <TH2.h>
#include <TStyle.h>

using namespace std;

void HGCalSelector::Begin(TTree * /*tree*/)
{
   // The Begin() function is called at the start of the query.
   // When running with PROOF Begin() is only called on the client.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();
}

void HGCalSelector::SlaveBegin(TTree * /*tree*/)
{
   // The SlaveBegin() function is called after the Begin() function.
   // When running with PROOF SlaveBegin() is called on each slave server.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();
}

Bool_t HGCalSelector::Process(Long64_t entry)
{
   //fReader.SetEntry(entry); // TTreeReader seems to be broken
   GetEntry(entry);

   cout << event->run << ", " << event->evn << endl;

   std::ostringstream oss;
   oss << "events/event_" << event->evn << ".csv";
   testFile.open(oss.str());

   // rechits
   cout << rechits->size() << endl;
   for (unsigned i = 0; i < rechits->size(); ++i) {
       testFile << rechits->at(i).x << "," 
                << rechits->at(i).y << "," 
                << rechits->at(i).z  << "," 
                << rechits->at(i).energy 
                << "\n"; 
   }

   testFile.close();

   return kTRUE;
}

void HGCalSelector::SlaveTerminate()
{
   // The SlaveTerminate() function is called after all entries or objects
   // have been processed. When running with PROOF SlaveTerminate() is called
   // on each slave server.

}

void HGCalSelector::Terminate()
{
   // The Terminate() function is the last function to be called during
   // a query. It always runs on the client, it can be used to present
   // the results graphically or save the results to file.

}

int HGCalSelector::main() {
    return 1;
}
