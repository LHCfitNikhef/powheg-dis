
// Basic setup for Deeply Inelastic Scattering at HERA.

#include "Pythia8/Pythia.h"
#include <iostream>
#include <vector>
#include <fstream>
//#include "TH1F.h"
//#include "TFile.h"

using namespace Pythia8;

int main() {

  // Beam energies, minimal Q2, number of events to generate.
  double eProton   = 0.938;
  double eElectron = 1000;
  double Q2min     = 4;
  string pdfSet = "LHAPDF6:NNPDF40_nnlo_as_01180";
  string pdfPath;
  int    nEvent    = 1e6;

  // Generator. Shorthand for event.
  Pythia pythia;
  Event& event = pythia.event;
  pdfPath = pythia.settings.word("xmlPath") + "../NNPDF40_nnlo_as_01180";

  // Set up PDF set
  pythia.readString("PDF:pSet = " + pdfSet);

  // Set up incoming beams, for frame with unequal beam energies.
  pythia.readString("Beams:frameType = 2");
  // BeamA = proton.
  pythia.readString("Beams:idA = 2212");
  pythia.settings.parm("Beams:eA", eProton);
  // BeamB = electron.
  pythia.readString("Beams:idB = 12");
  pythia.settings.parm("Beams:eB", eElectron);

  // Set up DIS process within some phase space.
  // Neutral current (with gamma/Z interference).
  //pythia.readString("WeakBosonExchange:ff2ff(t:gmZ) = on");
  // Charged current.
  pythia.readString("WeakBosonExchange:ff2ff(t:W) = on");
  // Phase-space cut: minimal Q2 of process.
  pythia.settings.parm("PhaseSpace:Q2Min", Q2min);
  //pythia.settings.parm("PhaseSpace:xmin", 0.0002);
  //pythia.settings.parm("PhaseSpace:xmax", 0.00032);
  pythia.settings.parm("PhaseSpace:ymin", 0.04);
  pythia.settings.parm("PhaseSpace:ymax", 0.96);  


  // Set dipole recoil on. Necessary for DIS + shower.
  pythia.readString("SpaceShower:dipoleRecoil = on");

  // Allow emissions up to the kinematical limit,
  // since rate known to match well to matrix elements everywhere.
  pythia.readString("SpaceShower:pTmaxMatch = 2");

  // QED radiation off lepton not handled yet by the new procedure.
  pythia.readString("PDF:lepton = off");
  pythia.readString("TimeShower:QEDshowerByL = off");

  pythia.readString("PartonLevel:FSR = off"); 
  pythia.readString("TimeShower:QCDshower = off");
  pythia.readString("TimeShower:QEDshowerByQ = off");
  pythia.readString("TimeShower:QEDshowerByL = off"); 
  pythia.readString("TimeShower:QEDshowerByOther = off");
  pythia.readString("TimeShower:QEDshowerByGamma = off"); 
  pythia.readString("PartonLevel:ISR = off");
  pythia.readString("SpaceShower:QCDshower = off"); 
  pythia.readString("SpaceShower:QEDshowerByQ = off"); 
  pythia.readString("SpaceShower:QEDshowerByL = off"); 
  pythia.readString("PartonLevel:MPI = off");
  pythia.readString("HadronLevel:Hadronize = off");

  // Initialize.
  pythia.init();

  // Histograms.
  double Wmax = sqrt(2.* eProton * eElectron);
  Hist Q2hist("Q [GeV]", 25, 4, 1000);
  Hist Whist("W [GeV]", 100, 0., Wmax);
  Hist xhist("x", 1, 0.0002, 0.00032);
  Hist yhist("y", 100, 0., 1.);
  Hist pTehist("pT of scattered electron [GeV]", 100, 0., 1000.);
  Hist pTrhist("pT of radiated parton [GeV]", 100, 0., 50.);
  Hist pTdhist("ratio pT_parton/pT_electron", 100, 0., 5.);
  Hist Elhist("E_l [GeV]", 25, 25., 1000.);
  Hist thetalhist("theta_l [rad]", 25, 0., 0.025);


  // Begin event loop.
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {
    if (!pythia.next()) continue;

    // Four-momenta of proton, electron, virtual photon/Z^0/W^+- (charged current case W boson, but still pythia still calls it pPhoton)
    Vec4 pProton = event[1].p();
    Vec4 peIn    = event[4].p();
    Vec4 peOut   = event[6].p();
    Vec4 pPhoton = peIn - peOut;

    // Q2, W2, Bjorken x, y.
    double Q2    = - pPhoton.m2Calc();
    double W2    = (pProton + pPhoton).m2Calc();
    double x     = Q2 / (2. * pProton * pPhoton);
    double y     = (pProton * pPhoton) / (pProton * peIn);
    double E_l   = eElectron - y*eElectron;
    double theta = 2*asin(sqrt(Q2 / 4 / eElectron / E_l));


    // Fill kinematics histograms.
    Q2hist.fill( Q2 );
    Whist.fill( sqrt(W2) );
    xhist.fill( x );
    yhist.fill( y );
    pTehist.fill( event[6].pT() );
    Elhist.fill( E_l );
    thetalhist.fill( theta );
    

    // Normalize the histogram to be the cross sections in pb
    //double scaleFactor = 1.0 / nEvent;
    //Q2hist /= Q2hist / scaleFactor * 4.223 / Q2bin;

    // pT spectrum of partons being radiated in shower.
    for (int i = 0; i < event.size(); ++i) if (event[i].statusAbs() == 43) {
      pTrhist.fill( event[i].pT() );
      pTdhist.fill( event[i].pT() / event[6].pT() );
    }

  // End of event loop. Statistics and histograms.
  }
  pythia.stat();
  
  std::ofstream outputFile("output-CC-nu-1tev.txt"); // create a new output file or overwrite an existing one

  if (outputFile.is_open()) { // check if the file was opened successfully
    outputFile << Elhist << thetalhist << Q2hist << Whist << xhist << yhist << pTehist << pTrhist << pTdhist; // write data to the file
    outputFile.close(); // close the file when done
    std::cout << "Data was written to output.txt\n";
  }
  else {
    std::cerr << "Error opening file\n";
  }

  //cout Elhist << thetalhist << Q2hist << Whist << xhist << yhist << pTehist << pTrhist << pTdhist;
  std::cout << Q2hist;  
  
  HistPlot hpl("plotE_l");
  hpl.frame( "El-CC-nu-1tev-ycuts", "Number of events for final state lepton energy E_l", "$E_l$ (GeV)", "$N_{\\mathrm{event}}$");
  hpl.add( Elhist, ".", "$E_l$");
  // Use logarithmic y scale.
  hpl.plot( true);

  HistPlot hpl2("plotQ");
  hpl2.frame( "Q2-CC-nu-1tev-ycuts", "Number of events for Q2", "$Q^2$ (GeV^2)", "$N_{\\mathrm{event}}$");
  hpl2.add( Q2hist, ".", "$Q^2$");
  // Use logarithmic y scale.
  hpl2.plot( true);

  HistPlot hpl3("plottheta");
  hpl3.frame( "theta-CC-nu-1tev-ycuts", "Number of events for scattering angle", "$\theta_l$ (pb)", "$N_{\\mathrm{event}}$");
  hpl3.add( thetalhist, ".", "$\theta_l$");
  // Use logarithmic y scale.
  hpl3.plot( true);

  // Done.
  return 0;
}
