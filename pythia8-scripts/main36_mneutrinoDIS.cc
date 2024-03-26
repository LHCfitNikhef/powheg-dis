// main36.cc is a part of the PYTHIA event generator.
// Copyright (C) 2023 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Keywords: basic usage; DIS

// Basic setup for Deeply Inelastic Scattering at HERA.

#include "Pythia8/Pythia.h"
#include "Pythia8/PartonDistributions.h"
using namespace Pythia8;

int main() {

  // Beam energies, minimal Q2, number of events to generate.

  // HERA proton
  double eProton   = 0; // GeV
  double eElectron = 1e9 ; // GeV
  double Q2min     = 1.; // Only perturbative events
  int    nEvent    = 1e4;

  bool NC = false;

  // Generator. Shorthand for event.
  Pythia pythia;
  Event& event = pythia.event;

  // Set up incoming beams, for frame with unequal beam energies.
  pythia.readString("Beams:frameType = 2");
  // BeamA = proton.
  pythia.readString("Beams:idA = 2212");
  pythia.settings.parm("Beams:eA", eProton);
  // BeamB = electron.
  //pythia.readString("Beams:idB = 11");
  pythia.readString("Beams:idB = 12");
  pythia.settings.parm("Beams:eB", eElectron);

  // Set up DIS process within some phase space.
  // Neutral current (with gamma/Z interference).
  if(NC){
    pythia.readString("WeakBosonExchange:ff2ff(t:gmZ) = on");
  }
  // Charged current (W exchange)
  if(!NC){
    pythia.readString("WeakBosonExchange:ff2ff(t:W) = on");
  }
  // Phase-space cut: minimal Q2 of process.
  pythia.settings.parm("PhaseSpace:Q2Min", Q2min);

  // Set dipole recoil on. Necessary for DIS + shower.
  pythia.readString("SpaceShower:dipoleRecoil = on");

  // Allow emissions up to the kinematical limit,
  // since rate known to match well to matrix elements everywhere.
  pythia.readString("SpaceShower:pTmaxMatch = 2");

  // QED radiation off lepton not handled yet by the new procedure.
  pythia.readString("PDF:lepton = off");
  pythia.readString("TimeShower:QEDshowerByL = off");

  // Initialize.
  pythia.init();

  // Define bins
  // double xmin = 0.0002;
  // double xmax = 0.00032;
  // double Q2min_bin = 18;
  // double Q2max_bin = 20;

  // double xmin = 0.1;
  // double xmax = 0.11;
  // double Q2min_bin = 190;
  // double Q2max_bin = 230;

  //double xmin = 0.3;
  //double xmax = 0.32;
  //double Q2min_bin = 1e4;
  // double Q2max_bin = 1.4e4;

  double xmin = 1e-9;
  double xmax = 1;
  double Q2min_bin = 1;
  double Q2max_bin = 1e8;

  int event_count = 0;
  int event_bin_count = 0;

  // Begin event loop.
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {
    if (!pythia.next()) continue;

    // Four-momenta of proton, incoming electron, outgoing electron, virtual photon/Z^0/W^+-.
    Vec4 pProton = event[1].p();
    Vec4 peIn    = event[4].p();
    Vec4 peOut   = event[6].p();
    Vec4 pPhoton = peIn - peOut;

    // Q2, W2, Bjorken x, y.
    double Q2    = - pPhoton.m2Calc();
    double W2    = (pProton + pPhoton).m2Calc();
    double x     = Q2 / (2. * pProton * pPhoton);
    double y     = (pProton * pPhoton) / (pProton * peIn);
    
    // Count events in bin
    event_count ++;
    if( xmin < x && x < xmax){
      if( Q2min_bin < Q2 && Q2 < Q2max_bin){
	//	cout << "x = " << x << endl;
	// cout << "Q2 = " << Q2 << endl;
	event_bin_count++;
      }
    }


  // End of event loop. Statistics and histograms.
  }
  pythia.stat();
  
  // print
  cout<<endl;
  cout <<"Event fraction in bin = "<< double(event_bin_count)/event_count <<" ,Total number of bins = "<<event_count<<endl;
  cout <<"Events in bin = "<< event_bin_count <<endl;

  // xsec_incl
  // Would be good to get automatically from Pythia8
  double xsec_tot =  9.820e-06; // mb // Now it is hardwired
  double xsec_bin = xsec_tot * double(event_bin_count)/event_count; // mb
  // conversion to GeV
  // 1 GeV-2 = 0.389 mb
  // 1 mb = 2.5707 GeV^-2
  xsec_bin = xsec_bin *2.5707; // GeV^-2
  cout <<"\n xsec_acceptance [GeV^-2] = "<<xsec_bin <<endl;
  // 0.389379 pb =	1×10−9 GeV−2
  // 1 GeV^-2 = 0.389379e9 pb
  xsec_bin *=  0.389379e9;
  cout <<"\n xsec_acceptance [pb] = "<<xsec_bin <<endl;

  // Bin width
  xsec_bin = xsec_bin/(xmax-xmin)/(Q2max_bin-Q2min_bin); // GeV^-4

  cout <<"\n xsec_bin [GeV^-4] = "<<xsec_bin <<endl;

  // Neutral Current scattering
  if(NC){
    cout<<"\n NC scattering with electrons \n"<<endl;
    // Reduced xsec => NC
    double x = (xmax+xmin)/2;
    double Q2 = (Q2max_bin+Q2min_bin)/2; 
    double alpha_qed = 1/137.035999084;
    double prefactor = ( x*Q2*Q2 ) / ( 2*3.14159*pow(alpha_qed,2) );
    double s = 4 * eProton * eElectron;
    double y = Q2 / x / s;
    prefactor = prefactor / (1+pow(1-y,2) );
    cout<<"prefactor = "<<prefactor<<endl;
    
    double xsec_bin_red = xsec_bin * prefactor;  // Dimensionless

    cout<<"x = "<<x<<endl;
    cout<<"Q2 (GeV) = "<<Q2<<endl;
    cout<<"\n [Q2min, Q2max ] "<<Q2min_bin<<" , "<<Q2max_bin<<" GeV "<<endl;
    cout<<"\n [xmin, xmax ] "<<xmin<<" , "<<xmax<<endl;
    cout <<"  xsec_bin_red [Pythia8] = "<<xsec_bin_red <<endl;
    
  }

  // Charged Current scattering
  if(!NC){
    //cout<<"\n CC scattering with electrons \n"<<endl;
    cout<<"\n CC scattering with nuetrinos \n"<<endl;
    // Reduced xsec => CC
    double x = (xmax+xmin)/2;
    double Q2 = (Q2max_bin+Q2min_bin)/2;
    double s = 4 * eProton * eElectron;
    double y = Q2 / x / s;
    double const GF = 1.166e-5; // GeV^-2
    double const mw = 80.3692; // GeV
    double prefactor = ( 2 * 3.14159 * x )/pow(GF,2.0);
    prefactor =  prefactor * pow( Q2 + pow(mw,2.0) , 2.0)/ pow(mw, 4.0);
    cout<<"prefactor = "<<prefactor<<endl; // GeV2
    
    double xsec_bin_red = xsec_bin * prefactor;  // Dimensionless

    cout<<"x = "<<x<<endl;
    cout<<"y = "<<y<<endl;
    cout<<"Q2 (GeV) = "<<Q2<<endl;
    cout<<"\n [Q2min, Q2max ] "<<Q2min_bin<<" , "<<Q2max_bin<<" GeV "<<endl;
    cout<<"\n [xmin, xmax ] "<<xmin<<" , "<<xmax<<endl;
    cout <<"  xsec_bin_red [Pythia8] = "<<xsec_bin_red <<endl;
    
  }
  
  
  // Done.
  return 0;
}
