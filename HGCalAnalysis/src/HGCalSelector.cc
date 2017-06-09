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
   testFile.open("rechits.csv");
   testFile2.open("genparticles.csv");
   testFile3.open("multicluster.csv");

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
   cout << "\n ------------------------- \n";
   cout << "This event has run# and event# "<< event->run << ", " << event->evn << endl;
   //cout << event->run << ", " << event->evn << endl;
   //std::ostringstream oss;
   //oss << "events/event_" << event->evn << ".csv";
   //testFile.open(oss.str());
   //testFile.open("test.csv");
   

   // rechits
   //cout << "...Filling the event rechits, count of rechits is  " << rechits->size() << endl;
   for (unsigned i = 0; i < rechits->size(); ++i) {
       testFile << event->evn << ","
		        << rechits->at(i).x << "," 
                << rechits->at(i).y << "," 
                << rechits->at(i).z  << "," 
                << rechits->at(i).time  << "," 
                << rechits->at(i).energy 
                << "\n"; 
                
       cout << rechits->at(i).time << endl;

   }
   
   //cout << "...Filling the event genparticles, count of genparticles is  " << genparticles->size() << endl;
   for (unsigned i = 0; i < genparticles->size(); ++i) {
       testFile2 << event->evn << ","
		         << genparticles->at(i).pid << ","
		         << genparticles->at(i).eta << ","
                 << genparticles->at(i).phi << ","
                 << genparticles->at(i).pt  << ","
                 << genparticles->at(i).energy
                 << "\n";
   }

   //cout << "...Filling the event multicluster, count of multicluster is  " << multicluster->size() << endl;
   for (unsigned i = 0; i < multicluster->size(); ++i) {
       testFile3 << event->evn << "," 
		         << multicluster->at(i).eta << ","
                 << multicluster->at(i).phi << ","
                 << multicluster->at(i).z   << ","
		         << multicluster->at(i).pt  << ","
                 << multicluster->at(i).energy
                 << "\n";
   }
  // testFile.close();

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
   testFile.close();
   testFile2.close();
   testFile3.close();
}

int HGCalSelector::main() {
    return 1;
}
