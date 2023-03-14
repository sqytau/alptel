//
/// \brief Implementation of the HistoManager class
//

#include "HistoManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "HistoMessenger.hh"
#include "LXSetUp.hh"


HistoManager::HistoManager()
  : fFileName("lxphoton_out_vg1"), fHMessanger(0),
    fTreeCutX(1.0*mm), fTreeCutY(1.0*mm),
    fTreeParticle("gamma"), fvHitTrackList(0), fvTracks(0)
{
  Book();
  fHMessanger = new HistoMessenger(this);
  AllocateVolumeTracksMap();
}



HistoManager::~HistoManager()
{
  delete G4AnalysisManager::Instance();
  delete fHMessanger;
}


void  HistoManager::SetTracksVtx(const std::vector<G4ThreeVector> &vht) 
{
  fVtxx.clear();
  fVtxy.clear();
  fVtxz.clear();
  std::for_each(vht.begin(), vht.end(), [this](const G4ThreeVector v){ fVtxx.push_back(v.x()); 
                                                                   fVtxy.push_back(v.y()); 
                                                                   fVtxz.push_back(v.z());}); 
}


void  HistoManager::SetHitTrackPosition(const std::vector<G4ThreeVector> &vht) 
{
  ftrackx.clear();
  ftracky.clear();
  ftrackz.clear();
  std::for_each(vht.begin(), vht.end(), [this](const G4ThreeVector v){ ftrackx.push_back(v.x()); 
                                                                   ftracky.push_back(v.y()); 
                                                                   ftrackz.push_back(v.z());}); 
}

                  
void  HistoManager::SetTracksMomentum(const std::vector<G4ThreeVector> &vht) 
{
  fPx.clear();
  fPy.clear();
  fPz.clear();
  std::for_each(vht.begin(), vht.end(), [this](const G4ThreeVector v){ fPx.push_back(v.x()); 
                                                                   fPy.push_back(v.y()); 
                                                                   fPz.push_back(v.z());});
}



void HistoManager::AllocateVolumeTracksMap()
{
  const G4int nd = 16;
  const G4int id = 6;
  G4int did[nd] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};
  G4int iid[id] = {16, 17, 18, 19, 20, 21};
  for (G4int ii = 0; ii < nd; ++ii) fvolTrackDMap[did[ii]] = std::vector<G4double>();
  for (G4int ii = 0; ii < id; ++ii) fvolTrackIMap[iid[ii]] = std::vector<G4int>();
}



void HistoManager::SetVolumeTracks(const std::vector<VolumeTrackTuple> &tvtracks)
{
  for (auto itr = fvolTrackDMap.begin(); itr != fvolTrackDMap.end(); ++itr) {
    std::vector<G4double> &vv  =  itr->second;
    vv.clear();
  }
  for (auto itr = fvolTrackIMap.begin(); itr != fvolTrackIMap.end(); ++itr) {
    std::vector<G4int> &vv  =  itr->second;
    vv.clear();
  }

  for (const auto &ttrck : tvtracks) {
    fvolTrackDMap[0].push_back(std::get<0>(ttrck));
    fvolTrackDMap[1].push_back(std::get<1>(ttrck));
    fvolTrackDMap[2].push_back(std::get<2>(ttrck));
    fvolTrackDMap[3].push_back(std::get<3>(ttrck));
    fvolTrackDMap[4].push_back(std::get<4>(ttrck));
    fvolTrackDMap[5].push_back(std::get<5>(ttrck));
    fvolTrackDMap[6].push_back(std::get<6>(ttrck));
    fvolTrackDMap[7].push_back(std::get<7>(ttrck));
    fvolTrackDMap[8].push_back(std::get<8>(ttrck));
    fvolTrackDMap[9].push_back(std::get<9>(ttrck));
    fvolTrackDMap[10].push_back(std::get<10>(ttrck));
    fvolTrackDMap[11].push_back(std::get<11>(ttrck));
    fvolTrackDMap[12].push_back(std::get<12>(ttrck));
    fvolTrackDMap[13].push_back(std::get<13>(ttrck));
    fvolTrackDMap[14].push_back(std::get<14>(ttrck));
    fvolTrackDMap[15].push_back(std::get<15>(ttrck));
    fvolTrackIMap[16].push_back(std::get<16>(ttrck));
    fvolTrackIMap[17].push_back(std::get<17>(ttrck));
    fvolTrackIMap[18].push_back(std::get<18>(ttrck));
    fvolTrackIMap[19].push_back(std::get<19>(ttrck));
    fvolTrackIMap[20].push_back(std::get<20>(ttrck));
    fvolTrackIMap[21].push_back(std::get<21>(ttrck));
    fvolTrackDMap[22].push_back(std::get<22>(ttrck));
  }
}



void HistoManager::SetTrajectory(const std::vector<TrajectoryTuple> &trjectory)
{
  for (auto itr = fTrajectoryD.begin(); itr != fTrajectoryD.end(); ++itr) {
    std::vector<G4double> &vv  =  itr->second;
    vv.clear();
  }
  
  for (auto itr = fTrajectoryI.begin(); itr != fTrajectoryI.end(); ++itr) {
    std::vector<G4int> &vv  =  itr->second;
    vv.clear();
  }
  
  fTrajectoryS.clear();
  
  for (const auto &trjp : trjectory) {
    fTrajectoryD[0].push_back(std::get<0>(trjp));
    fTrajectoryD[1].push_back(std::get<1>(trjp));
    fTrajectoryD[2].push_back(std::get<2>(trjp));
    fTrajectoryD[3].push_back(std::get<3>(trjp));
    fTrajectoryD[4].push_back(std::get<4>(trjp));
    fTrajectoryD[5].push_back(std::get<5>(trjp));
    fTrajectoryD[6].push_back(std::get<6>(trjp));
    fTrajectoryD[7].push_back(std::get<7>(trjp));
    fTrajectoryD[8].push_back(std::get<8>(trjp));
    fTrajectoryI[0].push_back(std::get<9>(trjp));
    fTrajectoryI[1].push_back(std::get<10>(trjp));
  }
}



void HistoManager::Book()
{
  // Create or get analysis manager
  // The choice of analysis technology is done via selection of a namespace
  // in HistoManager.hh
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  analysisManager->SetFileName(fFileName);
  analysisManager->SetVerboseLevel(1);
  analysisManager->SetActivation(true);    // enable inactivation of histograms
  analysisManager->SetNtupleMerging(true);
  analysisManager->SetHistoDirectoryName("hist");

  // Define histograms start values
  const G4int kMaxHisto = 71;
  const G4String id[] = { "0", "1", "2", "3", "4", "5", "6", "7", "8", "9",
                         "10","11","12","13","14","15","16","17","18","19",
                         "20","21","22","23","24","25","26","27","28","29",
                         "30","31","32","33","34","35","36","37","38","39",
                         "40","41","42","43","44","45","46","47","48","49",
                         "50","51","52", "53","54","55","56","57","58","59",
                         "60","61","62", "63","64","65","66","67","68","69",
                         "70"
                        };
                        
  const G4String title[] =
                { "number of events processed",                           //0
                  "energy deposit in absorber",                           //1
                  "energy of charged secondaries at creation",            //2
                  "energy of neutral secondaries at creation",            //3
                  "energy of charged at creation (log scale)",            //4
                  "energy of neutral at creation (log scale)",            //5
                  "x_vertex of charged secondaries (all)",                //6
                  "x_vertex of charged secondaries (not absorbed)",       //7
                  "dummy","dummy",                                        //8-9
                  "(transmit, charged) : kinetic energy at exit",         //10
                  "(transmit, charged) : ener fluence: dE(MeV)/dOmega",   //11
                  "(transmit, charged) : space angle: dN/dOmega",         //12
                  "(transmit, charged) : projected angle at exit",        //13
                  "(transmit, charged) : projected position at exit",     //14
                  "(transmit, charged) : radius at exit",                 //15
                  "energy of Auger e- at creation",                       //16
                  "energy of fluorescence gamma at creation",             //17
                  "energy of Auger e- at creation (log scale)",           //18
                  "energy of fluorescence gamma at creation (log scale)", //19
                  "(transmit, neutral) : kinetic energy at exit",         //20
                  "(transmit, neutral) : ener fluence: dE(MeV)/dOmega",   //21
                  "(transmit, neutral) : space angle: dN/dOmega",         //22
                  "(transmit, neutral) : projected angle at exit",        //23
                  "gamma: x-position at exit",                            //24
                  "gamma: radius at exit",                                //25
                  "dummy","dummy","dummy","dummy",       //26-29
                  "(reflect , charged) : kinetic energy at exit",         //30
                  "(reflect , charged) : ener fluence: dE(MeV)/dOmega",   //31
                  "(reflect , charged) : space angle: dN/dOmega",         //32
                  "(reflect , charged) : projected angle at exit",        //33
                  "dummy","dummy","dummy","dummy","dummy","dummy",       //34-39
                  "(reflect , neutral) : kinetic energy at exit",         //40
                  "(reflect , neutral) : ener fluence: dE(MeV)/dOmega",   //41
                  "(reflect , neutral) : space angle: dN/dOmega",         //42
                  "(reflect , neutral) : projected angle at exit",        //43
                  "energy of PIXE Auger e- at creation",                  //44
                  "energy of PIXE gamma at creation",                     //45
                  "energy of PIXE Auger e- at creation (log scale)",      //46
                  "energy of PIXE gamma at creation (log scale)",         //47
                  "dummy","dummy",                                        //48-49
                  "dummy", 
                  "Z vertex position of charged produced after the magnet", //51
                  "Z vertex position of neutral produced after the magnet", //52

                  "X vertex position of primary particles", //53
                  "Y vertex position of primary particles", //54
                  "Z vertex position of primary particles", //55
                  
                  "PX of primary particles", //56
                  "PY of primary particles", //57
                  "PZ of primary particles", //58
                  
                  "X at focus of primary particles", //59
                  "Y at focus of primary particles", //60

                  "electron time at exit", //61
                  "gamma time at exit", //62

                  "dummy","dummy",
                  "electron spectrum",  //65
                  "photon spectrum",    //66
                  "positron spectrum",  //67
                  
                  "electron x",  //68
                  "photon x",    //69
                  "positron x"  //70
                    
                };

  // Default values (to be reset via /analysis/h1/set command)
  G4int nbins = 100;
  G4double vmin = 0.;
  G4double vmax = 100.;
  LXSetUp *lxs = LXSetUp::Instance();

  // Create all histograms as inactivated 
  // as we have not yet set nbins, vmin, vmax
  for (G4int k=0; k<kMaxHisto; k++) {
    G4int ih = analysisManager->CreateH1("h"+id[k], title[k], nbins,vmin,vmax);
    analysisManager->SetH1Activation(ih, false);
  }
  
// ob  
  const G4int kMaxHisto2D = 13;
  const G4String id2D[] = { "0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12"
                        };
                        
  const G4String title2D[] =
                { "(transmit, charged) : projected position at exit vs energy",  //0
                  "(transmit, charged vs neutral): kinetic energy e-gamma at exit", //1
                  "(transmit, charged): x, y position at exit", //2
                  "(transmit, neutral): x, y position at exit", //3
                  "Z vertex position vs E of charged produced after the magnet", //4
                  "Z vertex position vs E of neutral produced after the magnet", //5
                  "electron polar anvgle vs E at exit",     // 6
                  "photons polar anvgle vs E at exit",      // 7
                  "photons polar anvgle vs E at exit, urad scale",      // 8
//                   "photons radial position vs momentum theta at exit"   //9
                  "photons position x, y at exit",   //9
                  "primary electron x, vs dx/dz",   //10
//                   "primary electron y, vs dy/dz"   //11
                  "electron x, vs dx/dz at IP",   //11
                  "positron polar anvgle vs E at exit"   //12
                };
  
  for (G4int k=0; k<kMaxHisto2D; k++) {
    G4int ih = analysisManager->CreateH2("h2d"+id2D[k], title2D[k], nbins,vmin,vmax, nbins,vmin,vmax);
    analysisManager->SetH2Activation(ih, false);
    G4cout << "2D histogram id: " << ih << G4endl; }

    
      if(lxs->ScintCerenkovPhysics){
    analysisManager->CreateH1("h71", "Optical photons produced in Brem. Cerenkov Channels (Channel ID)", lxs->CerenkovChannels, 0, lxs->CerenkovChannels);
    analysisManager->CreateH1("h72", "Optical photons produced in IP Cerenkov Channels (Channel ID)", lxs->CerenkovChannels, 0, lxs->CerenkovChannels);
    analysisManager->CreateH1("h73", "Optical photons produced in gamma spectrometer e+ side Cerenkov Channels (Channel ID)", lxs->CerenkovChannels, 0, lxs->CerenkovChannels);
    analysisManager->CreateH1("h74", "Optical photons produced in gamma spectrometer e- side Cerenkov Channels (Channel ID)", lxs->CerenkovChannels, 0, lxs->CerenkovChannels);
  
    analysisManager->CreateH1("h75", "Optical photons produced in Brem. Scintillator (mm)", 8000, -lxs->ScintX/2.*mm, lxs->ScintX/2.*mm);
    analysisManager->CreateH1("h76", "Optical photons produced in IP Scintillator (mm)", 8000, -lxs->ScintX/2.*mm, lxs->ScintX/2.*mm);
    analysisManager->CreateH1("h77", "Optical photons produced in gamma spectrometer e+ side Scintillator (mm)", 8000, -lxs->ComptonLysoX/2.*mm, lxs->ComptonLysoX/2.*mm);
    analysisManager->CreateH1("h78", "Optical photons produced in gamma spectrometer e- side Scintillator (mm)", 8000, -lxs->ComptonLysoX/2.*mm, lxs->ComptonLysoX/2.*mm);
    
    analysisManager->CreateH1("h79", "Optical photons registered in Brem. Scintillator Cameras (Camera ID)", 2, 0, 2);
    analysisManager->CreateH1("h80", "Optical photons registered in IP Scintillator Cameras (Camera ID)", 2, 0, 2);

    // gamma spectrometer cameras do not exist yet

    //    analysisManager->CreateH1("h81", "Optical photons registered in gamma spectrometer e+ side Scintillator Cameras", 2, 0, 2);
    // analysisManager->CreateH1("h82", "Optical photons registered in gamma spectrometer e- side Scintillator Cameras", 2, 0, 2);

    // for (int i= 71; i<83; i++){
    //   analysisManager->SetH1Activation(ih, true);}

    // H2D with mm^2 precision, 1-D histos have nominal 125 um pixel resolution. ROOT / G4 crashes (for me) for 8000 * 800 H2D histos. - John
    
    analysisManager->CreateH2("h2d13", "Optical photons produced in Brem. Scintillator (mm / mm)", 1000, -lxs->ScintX/2.*mm, lxs->ScintX/2.*mm, 100, -lxs->ScintY/2.*mm, lxs->ScintY/2.*mm);
    analysisManager->CreateH2("h2d14", "Optical photons produced in IP Scintillator (mm / mm)", 1000, -lxs->ScintX/2.*mm, lxs->ScintX/2.*mm, 100, -lxs->ScintY/2.*mm, lxs->ScintY/2.*mm);
    analysisManager->CreateH2("h2d15", "Optical photons produced in gamma spectrometer e+ side Scintillator (mm / mm)", 1000, -lxs->ComptonLysoX/2.*mm, lxs->ComptonLysoX/2.*mm, 50, -lxs->ComptonLysoY/2.*mm, lxs->ComptonLysoY/2.*mm);
    analysisManager->CreateH2("h2d16", "Optical photons produced in gamma spectrometer e- side Scintillator (mm / mm)", 1000, -lxs->ComptonLysoX/2.*mm, lxs->ComptonLysoX/2.*mm, 100, -lxs->ComptonLysoY/2.*mm, lxs->ComptonLysoY/2.*mm);

    
    analysisManager->CreateH2("h2d17", "Optical photons produced in Brem. Cherenkov (mm / mm)", 1000, -lxs->CerenkovTotalBoxLength/2., lxs->CerenkovTotalBoxLength/2., 100, -lxs->CerenkovTotalBoxHeight/2., lxs->CerenkovTotalBoxHeight/2.);
    analysisManager->CreateH2("h2d18", "Optical photons produced in IP Cherenkov (mm / mm)", 1000, -lxs->CerenkovTotalBoxLength/2., lxs->CerenkovTotalBoxLength/2., 100, -lxs->CerenkovTotalBoxHeight/2., lxs->CerenkovTotalBoxHeight/2.);
    analysisManager->CreateH2("h2d19", "Optical photons produced in gamma spectrometer e+ side Cherenkov (mm / mm)", 1000, -lxs->CerenkovTotalBoxLength/2., lxs->CerenkovTotalBoxLength/2., 100, -lxs->CerenkovTotalBoxHeight/2., lxs->CerenkovTotalBoxHeight/2.);
    analysisManager->CreateH2("h2d20", "Optical photons produced in gamma spectrometer e- side Cherenkov (mm / mm)", 1000, -lxs->CerenkovTotalBoxLength/2., lxs->CerenkovTotalBoxLength/2., 100, -lxs->CerenkovTotalBoxHeight/2., lxs->CerenkovTotalBoxHeight/2.);


    analysisManager->CreateH2("h2d21", "Energy all non-scintillation particles in contact w// Brem. Scintillator (GeV / mm)", 1000, -lxs->ScintX/2.*mm, lxs->ScintX/2.*mm, 180, 0., 18.*GeV);
    analysisManager->CreateH2("h2d22", "Energy all non-scintillation particles in contact w// IP Scintillator (GeV / mm)", 1000, -lxs->ScintX/2.*mm, lxs->ScintX/2.*mm, 180, 0., 18.*GeV);
    analysisManager->CreateH2("h2d23", "Energy all non-scintillation particles in contact w// gamma spectrometer e+ side Scintillator (GeV / mm)", 1000, -lxs->ComptonLysoX/2.*mm, lxs->ComptonLysoX/2.*mm, 180, 0., 18.*GeV);
    analysisManager->CreateH2("h2d24", "Energy all non-scintillation particles in contact w// gamma spectrometer e- side Scintillator (GeV / mm)", 1000, -lxs->ComptonLysoX/2.*mm, lxs->ComptonLysoX/2.*mm, 180, 0., 18.*GeV);
    
    analysisManager->CreateH2("h2d25", "Energy all non-Cherenkov particles in contact w// Brem. Cherenkov Detector (GeV / mm)", 1000, -lxs->ScintX/2.*mm, lxs->ScintX/2.*mm, 180, 0., 18.*GeV);
    analysisManager->CreateH2("h2d26", "Energy all non-Cherenkov particles in contact w// IP Cherenkov Detector (GeV / mm)", 1000, -lxs->ScintX/2.*mm, lxs->ScintX/2.*mm, 180, 0., 18.*GeV);
    analysisManager->CreateH2("h2d27", "Energy all non-Cherenkov particles in contact w// gamma spectrometer e+ side  Cherenkov Detector (GeV / mm)", 1000, -lxs->CerenkovTotalBoxLength/2.*mm, lxs->CerenkovTotalBoxLength/2.*mm,  180, 0., 18.*GeV);
    analysisManager->CreateH2("h2d28", "Energy all non-Cherenkov particles in contact w// gamma spectrometer e- side Cherenkov Detector (GeV / mm)", 1000, -lxs->CerenkovTotalBoxLength/2.*mm, lxs->CerenkovTotalBoxLength/2.*mm,  180, 0., 18.*GeV);
    
    // for (int i= 13; i<17; i++){
    //   analysisManager->SetH2Activation(ih, true);}

    
      }

    

  // Creating ntuple
  //
  analysisManager->CreateNtuple("lxtsim", "Bremsstrahlung photons at IP");
  analysisManager->CreateNtupleDColumn(0, "E");
  analysisManager->CreateNtupleDColumn(0, "x");
  analysisManager->CreateNtupleDColumn(0, "y");
  analysisManager->CreateNtupleDColumn(0, "z");
  analysisManager->CreateNtupleDColumn(0, "px");
  analysisManager->CreateNtupleDColumn(0, "py");
  analysisManager->CreateNtupleDColumn(0, "pz");
  analysisManager->CreateNtupleDColumn(0, "theta");
  analysisManager->CreateNtupleDColumn(0, "phi");
  analysisManager->CreateNtupleIColumn(0, "pdg");
  analysisManager->CreateNtupleIColumn(0, "physproc");
  analysisManager->CreateNtupleDColumn(0, "weight");
  analysisManager->FinishNtuple(0);

  analysisManager->CreateNtuple("Tracks", "Tracks hitting volumes marked for track interception");
  analysisManager->CreateNtupleIColumn(1, "eventid");                               //0
  analysisManager->CreateNtupleIColumn(1, "trackid",  fvolTrackIMap[16]);           //1
  analysisManager->CreateNtupleIColumn(1, "detid",    fvolTrackIMap[17]);           //2
  analysisManager->CreateNtupleIColumn(1, "pdg",      fvolTrackIMap[18]);           //3
  analysisManager->CreateNtupleIColumn(1, "physproc", fvolTrackIMap[19]);           //4
  analysisManager->CreateNtupleDColumn(1, "E", fvolTrackDMap[0]);                   //5
  analysisManager->CreateNtupleDColumn(1, "x", fvolTrackDMap[1]);                   //6
  analysisManager->CreateNtupleDColumn(1, "y", fvolTrackDMap[2]);                   //7
  analysisManager->CreateNtupleDColumn(1, "z", fvolTrackDMap[3]);                   //8
  analysisManager->CreateNtupleDColumn(1, "t", fvolTrackDMap[4]);                   //9
  analysisManager->CreateNtupleDColumn(1, "vtxx", fvolTrackDMap[5]);                //10
  analysisManager->CreateNtupleDColumn(1, "vtxy", fvolTrackDMap[6]);                //11
  analysisManager->CreateNtupleDColumn(1, "vtxz", fvolTrackDMap[7]);                //12
  analysisManager->CreateNtupleDColumn(1, "px",   fvolTrackDMap[8]);                //13
  analysisManager->CreateNtupleDColumn(1, "py",   fvolTrackDMap[9]);                //14
  analysisManager->CreateNtupleDColumn(1, "pz",   fvolTrackDMap[10]);               //15
  analysisManager->CreateNtupleDColumn(1, "theta",  fvolTrackDMap[11]);             //16
  analysisManager->CreateNtupleDColumn(1, "phi",    fvolTrackDMap[12]);             //17
  analysisManager->CreateNtupleDColumn(1, "xlocal", fvolTrackDMap[13]);             //18
  analysisManager->CreateNtupleDColumn(1, "ylocal", fvolTrackDMap[14]);             //19
  analysisManager->CreateNtupleDColumn(1, "zlocal", fvolTrackDMap[15]);             //20
  analysisManager->CreateNtupleDColumn(1, "weight");                                //21
  analysisManager->CreateNtupleIColumn(1, "ptrackid", fvolTrackIMap[20]);           //22
  analysisManager->CreateNtupleIColumn(1, "nsecondary", fvolTrackIMap[21]);         //23
  analysisManager->CreateNtupleDColumn(1, "esecondary", fvolTrackDMap[22]);         //24
  analysisManager->CreateNtupleIColumn(1, "runid");                                 //25
  analysisManager->FinishNtuple(1);

  analysisManager->CreateNtuple("Hits", "Hits in sensitive detectors");
  analysisManager->CreateNtupleIColumn(2, "eventid");
  analysisManager->CreateNtupleIColumn(2, "detid");
  analysisManager->CreateNtupleIColumn(2, "layerid");
  analysisManager->CreateNtupleIColumn(2, "cellx");
  analysisManager->CreateNtupleIColumn(2, "celly");
  analysisManager->CreateNtupleDColumn(2, "edep");
  analysisManager->CreateNtupleIColumn(2, "hitid");
  analysisManager->CreateNtupleIColumn(2, "track_list", fvHitTrackList);
  analysisManager->CreateNtupleDColumn(2, "trackx", ftrackx);
  analysisManager->CreateNtupleDColumn(2, "tracky", ftracky);
  analysisManager->CreateNtupleDColumn(2, "trackz", ftrackz);
  analysisManager->CreateNtupleDColumn(2, "trackt", ftrackt);
  analysisManager->CreateNtupleDColumn(2, "trackedep", ftrackedep);
  analysisManager->CreateNtupleDColumn(2, "weight");
  analysisManager->CreateNtupleIColumn(2, "runid");
  analysisManager->FinishNtuple(2);

  analysisManager->CreateNtuple("HitTracks", "Tracks which produced hits in sensitive detectors");
  analysisManager->CreateNtupleIColumn(3, "eventid");
  analysisManager->CreateNtupleIColumn(3, "trackid", fvTracks);
  analysisManager->CreateNtupleDColumn(3, "vtxx", fVtxx);
  analysisManager->CreateNtupleDColumn(3, "vtxy", fVtxy);
  analysisManager->CreateNtupleDColumn(3, "vtxz", fVtxz);
  analysisManager->CreateNtupleDColumn(3, "px", fPx);
  analysisManager->CreateNtupleDColumn(3, "py", fPy);
  analysisManager->CreateNtupleDColumn(3, "pz", fPz);
  analysisManager->CreateNtupleDColumn(3, "E", fE);
  analysisManager->CreateNtupleIColumn(3, "pdg", fPDG);
  analysisManager->CreateNtupleIColumn(3, "pproc", fPhysProc);
  analysisManager->CreateNtupleIColumn(3, "ptid", fPTId);
  analysisManager->CreateNtupleDColumn(3, "weight");
  analysisManager->CreateNtupleIColumn(3, "runid");
  analysisManager->FinishNtuple(3);

  analysisManager->CreateNtuple("DetSettings", "Sensitive detector settings");
  analysisManager->CreateNtupleIColumn(4, "detid");
  analysisManager->CreateNtupleSColumn(4, "det_name");
  analysisManager->CreateNtupleIColumn(4, "layerid");
  analysisManager->CreateNtupleIColumn(4, "n_cell_x");
  analysisManager->CreateNtupleIColumn(4, "n_cell_y");
  analysisManager->CreateNtupleDColumn(4, "size_x");
  analysisManager->CreateNtupleDColumn(4, "size_y");
  analysisManager->CreateNtupleDColumn(4, "size_z");
  analysisManager->CreateNtupleDColumn(4, "translation_x");
  analysisManager->CreateNtupleDColumn(4, "translation_y");
  analysisManager->CreateNtupleDColumn(4, "translation_z");
  analysisManager->CreateNtupleDColumn(4, "e_phi");
  analysisManager->CreateNtupleDColumn(4, "e_theta");
  analysisManager->CreateNtupleDColumn(4, "e_psi");
  analysisManager->CreateNtupleDColumn(4, "mass");
  analysisManager->CreateNtupleSColumn(4, "material");
  analysisManager->CreateNtupleDColumn(4, "density");
  analysisManager->FinishNtuple(4);
  
  analysisManager->CreateNtuple("Trajectory", "Track trajectories");
  analysisManager->CreateNtupleIColumn(5, "eventid");
  analysisManager->CreateNtupleIColumn(5, "trackid");
  analysisManager->CreateNtupleIColumn(5, "type");
  analysisManager->CreateNtupleIColumn(5, "pdg");
  analysisManager->CreateNtupleIColumn(5, "ptid");
  analysisManager->CreateNtupleDColumn(5, "weight");
  analysisManager->CreateNtupleDColumn(5, "x", fTrajectoryD[0]);
  analysisManager->CreateNtupleDColumn(5, "y", fTrajectoryD[1]);
  analysisManager->CreateNtupleDColumn(5, "z", fTrajectoryD[2]);
  analysisManager->CreateNtupleDColumn(5, "t", fTrajectoryD[3]);
  analysisManager->CreateNtupleDColumn(5, "E", fTrajectoryD[4]);
  analysisManager->CreateNtupleDColumn(5, "px", fTrajectoryD[5]);
  analysisManager->CreateNtupleDColumn(5, "py", fTrajectoryD[6]);
  analysisManager->CreateNtupleDColumn(5, "pz", fTrajectoryD[7]);
  analysisManager->CreateNtupleDColumn(5, "edep", fTrajectoryD[8]);
  analysisManager->CreateNtupleIColumn(5, "pproc",   fTrajectoryI[0]);
  analysisManager->CreateNtupleIColumn(5, "nsecond", fTrajectoryI[1]);
  analysisManager->CreateNtupleIColumn(5, "runid");
  analysisManager->FinishNtuple(5);

  analysisManager->CreateNtuple("RunInfo", "Run Information");
  analysisManager->CreateNtupleSColumn(6, "parameter");
  analysisManager->CreateNtupleSColumn(6, "value");
  analysisManager->FinishNtuple(6);

}
