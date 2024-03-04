#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <algorithm>

#include <TH1.h>
#include <TH2.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TF1.h>
#include <TF2.h>
#include <TTree.h>
#include <TTreeCache.h>
#include <TChain.h>
#include <TLorentzVector.h>
#include <TFile.h>
#include <TDirectory.h>
#include <TString.h>
#include <TMath.h>
#include <TStopwatch.h>
#include <TTimeStamp.h>
#include <TSystem.h>

void DEBUG(int i, TString str = "")
{
    cout << "DEBUG: " << i << "\t" << str << endl;
}

void printRunTime(TStopwatch timer_)
{
    Double_t cpuTime = timer_.CpuTime();
    Double_t realTime = timer_.RealTime();

    cout << endl;
    cout << "************************************************" << endl;
    cout << "Total real time: " << realTime << " (seconds)" << endl;
    cout << "Total CPU time:  " << cpuTime << " (seconds)" << endl;
    cout << "  CPU time / real time = " << cpuTime / realTime << endl;
    cout << "************************************************" << endl;
}

void printRunTimeShort(TStopwatch timer_)
{
    // Double_t cpuTime = timer_.CpuTime();
    Double_t realTime = timer_.RealTime();

    cout << "Total real time: " << realTime << " (seconds)" << endl;
    // cout << "Total CPU time:  " << cpuTime << " (seconds)" << endl;
}

static inline void printMemory( TString tab = "" )
{
    ifstream proc_status("/proc/self/status");
    string buffer;
    while (proc_status.peek() != EOF) {
        getline(proc_status, buffer);
        TString str = buffer;
        if(str.Contains("RSS")) {
            cout << tab << str << endl;
            break;
        }
    }
}

static inline void loadBar(int x, int n, int r, int w)
{
    // Only update r times.
    if( x == n )
        cout << endl;

    if ( x % (n/r +1) != 0 ) return;

    // Calculuate the ratio of complete-to-incomplete.
    float ratio = x/(float)n;
    int   c     = ratio * w;

    // Show the percentage complete.
    printf("%3d%% [", (int)(ratio*100) );

    // Show the load bar.
    for (int x=0; x<c; x++) cout << "=";

    for (int x=c; x<w; x++) cout << " ";

    // ANSI Control codes to go back to the
    // previous line and clear it.
    cout << "]\r" << flush;
}

// -- bin width of the histograms for res (line 764)
const static int n_eta_bins = 15;
double eta_bins[n_eta_bins] = {
  -2.4, -2.1, -1.6, -1.2, -0.9,
  -0.3, -0.2,  0.0,  0.2,  0.3,
   0.9,  1.2,  1.6,  2.1,  2.4
};

const static int n_pt_bins = 9;
double pt_bins[n_pt_bins] = {
     0, 2, 4, 6, 8, 
     12, 16, 20, 30
};

struct sort_by_pt
{
    inline bool operator() (const Object& a, const Object& b)
    {
        return (a.pt > b.pt);
    }
};

const double MU_MASS = 0.1056583745;

bool acceptance(Object obj)
{
    return ( fabs(obj.eta) < 2.4 );
}

double invMass(Object obj1, Object obj2)
{
    TLorentzVector mu0, mu1;
    mu0.SetPtEtaPhiM(obj1.pt,
                     obj1.eta,
                     obj1.phi,
                     MU_MASS);
    mu1.SetPtEtaPhiM(obj2.pt,
                     obj2.eta,
                     obj2.phi,
                     MU_MASS);
    return (mu0+mu1).M();
}