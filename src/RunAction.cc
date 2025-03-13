//
/// \brief Implementation of the RunAction class
//

#include <stdexcept>
#include <iomanip>
#include <vector>
#include <string>
#include <sys/stat.h>

#include "RunAction.hh"
#include "DetectorConstruction.hh"
#include "PrimaryGeneratorAction.hh"
#include "HistoManager.hh"
#include "Run.hh"

#include "G4AnalysisManager.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4EmCalculator.hh"

#include "Randomize.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcommand.hh"
#include "G4GDMLParser.hh"


#include "RunActionMessenger.hh"
#include "PrimarySpectra.hh"


RunAction::RunAction(DetectorConstruction* det, PrimaryGeneratorAction* kin)
:G4UserRunAction(),fDetector(det), fPrimary(kin), fRun(0), fHistoManager(0), fDumpGeometry(false),
 fSaveTrajectory(false), fExternalRunID(-1)
{ 
  // Book predefined histograms
  fHistoManager = new HistoManager();
  fRunActionMessenger = new RunActionMessenger(this);

}



RunAction::~RunAction()
{ 
  delete fHistoManager;
  delete fRunActionMessenger;
}



G4Run* RunAction::GenerateRun()
{ 
  fRun = new Run(fDetector);
  fRun->SetHistoManager(fHistoManager);
  fRun->SetInterceptVolumes(fInterceptVol);
  fRun->SetInterceptVolECut(fInterceptVolECut);
  fRun->SetSensitiveVolumes(fSensitiveVol);
  fRun->SetSaveTrajectory(fSaveTrajectory);
  return fRun;
}



void RunAction::BeginOfRunAction(const G4Run*)
{
  // save Rndm status
  ////  G4RunManager::GetRunManager()->SetRandomNumberStore(true);
  if (isMaster) G4Random::showEngineStatus();
     
  // keep run condition
  if ( fPrimary ) { 
    G4ParticleDefinition* particle = fPrimary->GetParticleGun()->GetParticleDefinition();
    G4double energy = fPrimary->GetParticleGun()->GetParticleEnergy();
    fRun->SetPrimary(particle, energy);
    fRun->SetSkipEvents(fPrimary->GetSkipEvents());

    G4int btype = fPrimary->GetBeamType();
      if (btype == PrimaryGeneratorAction::beamMC || btype == PrimaryGeneratorAction::beamMCh5) {
      std::string fmcname = fPrimary->GetMCfile();
      std::hash<std::string> mcfhash;
      fRun->SetRunID(mcfhash(fmcname));
    } else if (btype == PrimaryGeneratorAction::beamMono) {
        if (fExternalRunID > 0) fRun->SetRunID(fExternalRunID);
        else fRun->SetRunID(1);
    } else if (btype == PrimaryGeneratorAction::beamGauss) {
        if (fExternalRunID > 0) fRun->SetRunID(fExternalRunID);
        else fRun->SetRunID(2);
    }
  }
  
  //histograms
  //        
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  if ( analysisManager->IsActive() ) {
    PrimaryGeneratorAction *gen = fPrimary;
    if ( gen ) { 
//      G4double zmax = fabs(0.5 * fDetector->GetWorldSizeZ() - gen->GetZ0()); 
      G4double zmax = fabs(gen->GetBeamFocus() - gen->GetZ0()); 
   
      G4double beta = gen->GetBetaX();
      G4double range = 5.0*gen->GetSigmaX() * sqrt(1.0 + pow(zmax/beta, 2.0) ); 
      G4double drange = range*zmax/(zmax*zmax + beta*beta) + 5.0*sqrt((gen->GetEmittanceX()*beta)/(zmax*zmax + beta*beta));
      analysisManager->SetH2Activation(10, true);
      analysisManager->SetH2(10, 500, -range, range, 500, -drange, drange);  
      analysisManager->SetH1Activation(53, true);
      analysisManager->SetH1(53, 1000, -range, range);  
      G4cout << "x range : " << range << "   " << drange << G4endl;

      beta = gen->GetBetaY();
      range = 5.0*gen->GetSigmaY() * sqrt(1.0 + pow(zmax/beta, 2.0) ); 
      drange = range*zmax/(zmax*zmax + beta*beta) + 5.0*sqrt((gen->GetEmittanceY()*beta)/(zmax*zmax + beta*beta));
      analysisManager->SetH2Activation(11, true);
      analysisManager->SetH2(11, 500, -range, range, 500, -drange, drange);  
      analysisManager->SetH1Activation(54, true);
      analysisManager->SetH1(54, 1000, -range, range);  
      G4cout << "y range : " << range << "   " << drange << G4endl;
      
      range = zmax/c_light;
      drange = 10.0*gen->GetSigmaZ()/c_light;
      analysisManager->SetH1Activation(61, true);
      analysisManager->SetH1(61, 1000, range - drange, range + drange);  
      analysisManager->SetH1Activation(62, true);
      analysisManager->SetH1(62, 1000, range - drange, range+drange);  
      G4cout << "t range : " << range - drange << "   " << range + drange << G4endl;
    }
    analysisManager->OpenFile();
  }

  if (!isMaster) {
    fRun->ScanGeometry(fDetector->GetphysiWorld(), 1);
  }

  if (isMaster) {
    if (fDumpGeometry) {
      G4GDMLParser parser;
      G4String fname("lxgeomdump.gdml");
      struct stat finfo;
      int i = 0;
      while (stat(fname.c_str(), &finfo) == 0) {
        fname = G4String("lxgeomdump_") + G4UIcommand::ConvertToString(i++) + G4String(".gdml");
      }
      parser.Write(fname, fDetector->GetphysiWorld());
    }
  }

  if (!isMaster) {
    SaveRunInfo();
    PrintRunInfo();
  }
}



void RunAction::EndOfRunAction(const G4Run*)
{  
  // print Run summary
  //
  if (isMaster) fRun->EndOfRun();    
      
  // save histograms
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();  
  if ( analysisManager->IsActive() ) {    
    analysisManager->Write();
    analysisManager->CloseFile();
  }  

  // show Rndm status
  if (isMaster) G4Random::showEngineStatus();
}



void RunAction::AddInterceptVolume (const G4String volumename)
{
  G4String vname;
  G4int    vid = -1;
  std::stringstream vinfstrm(volumename.data());

  std::vector<G4int>  voldepthvec;
  vinfstrm >> vname >> vid;
  voldepthvec.push_back(vid);
  G4int  voldepth;
  while (vinfstrm >> voldepth) voldepthvec.push_back(voldepth);

  G4int vdblid = -1;
  for (const auto &kv : fInterceptVol) { if (kv.second.at(0) == vid) vdblid = kv.second.at(0); }
  if (vdblid > 0) {
    G4String msgstr("Warning: Volume index ");
    msgstr += G4UIcommand::ConvertToString(vid);
    msgstr += G4String(" is already assigned to another volume! Automatic index will be used.");
    G4Exception("Run::", "AddInterceptVolume(volname)", JustWarning, msgstr.c_str());
    vid = -1;
  }

  if (fInterceptVol.find(vname) != fInterceptVol.end()) {
    G4String msgstr("Warning adding volume to intercept tracks. Volume ");
    msgstr += vname;
    msgstr += G4String(" is already in the list");
    G4Exception("Run::", "AddInterceptVolume(volname)", JustWarning, msgstr.c_str());
  } else {
    if (vid < 0) {
      G4int vmaxid = 0;
      for (const auto &kv : fInterceptVol) { if (kv.second.at(0) > vmaxid) vmaxid = kv.second.at(0); }
      G4int vinc(1);
      while (vmaxid/vinc > 10) vinc *= 10;
      vid = vmaxid + vinc;
      voldepthvec[0] = vid;
    }
    G4cout << "Volume " << vname << " will be added for track interception with id = " << vid << G4endl;
    if (voldepthvec.size()>1) {
      G4cout << "Depth levels up for volume copy numbers with respect volume level: ";
      std::for_each(voldepthvec.begin()+1, voldepthvec.end(), [](const G4int dpt) {G4cout << "  " << dpt; });
      G4cout << G4endl;
    }
    fInterceptVol[vname] = voldepthvec;
  }
}



void RunAction::AddSensitiveVolume (const G4String volumename)
{
  G4String volname(volumename);
  G4String  delim(":");
  G4int vinc(1000);
  G4String  vname("");
  G4int     vid(-1), dsd(1), dl(0);

  size_t  dpos, dpos0(0);
  std::vector <size_t> vdelimpos;
  G4StrUtil::strip(volname);
  do {
    dpos = volname.find(delim, dpos0);
    if (dpos != std::string::npos) {
      vdelimpos.push_back(dpos);
      dpos0 = dpos + 1;
    }
  } while (dpos != std::string::npos);
  vdelimpos.push_back(volname.length());

  G4String tmp;
  switch (vdelimpos.size()) {
    case 4 :
      tmp = volname.substr(vdelimpos[2]+1, vdelimpos[3]-vdelimpos[2]-1);
      G4StrUtil::strip(tmp);
      dl = G4UIcmdWithAnInteger::GetNewIntValue(tmp);
      [[fallthrough]];
    case 3 :
      tmp = volname.substr(vdelimpos[1]+1, vdelimpos[2]-vdelimpos[1]-1);
      G4StrUtil::strip(tmp);
      dsd = G4UIcmdWithAnInteger::GetNewIntValue(tmp);
      [[fallthrough]];
    case 2 :
      tmp = volname.substr(vdelimpos[0]+1, vdelimpos[1]-vdelimpos[0]-1);
      G4StrUtil::strip(tmp);
      vid = G4UIcmdWithAnInteger::GetNewIntValue(tmp);
      [[fallthrough]];
    case 1:
      vname = volname.substr(0, vdelimpos[0]);
      G4StrUtil::strip(vname);
      break;
    default :
      G4String msgstr("Warning: AddSensitiveVolume: Line ");
      msgstr += volumename;
      msgstr += G4String(" can not be processed! Ignore!");
      G4Exception("Run::", "AddInterceptVolume(volname)", JustWarning, msgstr.c_str());
  }

  G4int vmax = 0;
  bool matchid(false);
  for (const auto &svt : fSensitiveVol) {
    if (vmax < std::get<0>(svt.second))  vmax = std::get<0>(svt.second);
    if (vid == std::get<0>(svt.second))  matchid = true;
  }
  if (matchid) {
    G4String msgstr("Warning: AddSensitiveVolume: Detector ID ");
    msgstr += G4UIcommand::ConvertToString(vid);
    msgstr += G4String(" is already used! Automatic number will be assigned.");
    G4Exception("Run::", "AddInterceptVolume(volname)", JustWarning, msgstr.c_str());
  }
  if (vid == 0) {
    G4String msgstr("Warning: AddSensitiveVolume: Detector ID ");
    msgstr += G4UIcommand::ConvertToString(vid);
    msgstr += G4String(" is reserved! Automatic number will be assigned.");
    G4Exception("Run::", "AddInterceptVolume(volname)", JustWarning, msgstr.c_str());
  }
  if (vid <= 0 || matchid) {
    vid = vmax + vinc;
  }

  if (vname.size()>0) {
    fSensitiveVol[vname] = std::make_tuple (vid, dsd, dl);
    G4cout << "Sensitive Volume: " << vname << " will be added with id = " << vid
           << " subdetector depth: " << dsd << " and layer depth: " << dl << G4endl;
  }

  try {
    const auto &sseg = fDetector->GetSensorSegmentation(vname);
  } catch (const std::out_of_range& oor)
  {
    G4cout << "Segmentation is not defined for sensitive detector " << vname << "! Use default." << G4endl;
    fDetector->AddSensorSegmentation(vname, fDetector->GetWorldSizeXY(), fDetector->GetWorldSizeXY(), 1, 1);
  }
}



void RunAction::AddInterceptVolumeECut (const G4String volname, const G4double ecutt)
{
  const auto itr = fInterceptVol.find(volname);
  if (itr != fInterceptVol.end()) {
    G4cout << "Setting energy cut for volume " << volname << "  with ID " << itr->second.at(0)
           << " to " << ecutt/GeV << " GeV" << G4endl;
    fInterceptVolECut[itr->second.at(0)] = ecutt;
  } else {
    G4String msgstr("Warning: AddInterceptVolumeECut: Volume name ");
    msgstr += volname;
    msgstr += G4String(" is not specified for track inteception. Ignore this setting!");
    G4Exception("Run::", "AddInterceptVolume(volname)", JustWarning, msgstr.c_str());
  }
}


void RunAction::SaveRunInfo()
{
  if (!fPrimary) return;
  const G4ParticleGun *pgun = fPrimary->GetParticleGun();

  G4int btype = fPrimary->GetBeamType();
  if (btype == PrimaryGeneratorAction::beamMC || btype == PrimaryGeneratorAction::beamMCh5
      || btype == PrimaryGeneratorAction::beamMCTupleG4) {
    SaveParameterInfo("BeamType", "MC");
    SaveParameterInfo("MCFile", fPrimary->GetMCfile());

  } else if (btype == PrimaryGeneratorAction::beamMono) {
    SaveParameterInfo("BeamType", "Mono");
    SaveParameterInfo("BeamPosition",
       G4UIcommand::ConvertToString(G4ThreeVector(fPrimary->GetX0(), fPrimary->GetY0(), fPrimary->GetZ0()), "mm"));
    SaveParameterInfo("Particle", pgun->GetParticleDefinition()->GetParticleName());
    SaveParameterInfo("Energy", G4UIcommand::ConvertToString(pgun->GetParticleEnergy(), "GeV"));

  } else if (btype == PrimaryGeneratorAction::beamGauss) {
    SaveParameterInfo("BeamType", "Gauss");
    SaveParameterInfo("BeamPosition",
        G4UIcommand::ConvertToString(G4ThreeVector(fPrimary->GetX0(), fPrimary->GetY0(), fPrimary->GetZ0()), "mm"));
    SaveParameterInfo("BeamSigmaX", G4UIcommand::ConvertToString(fPrimary->GetSigmaX(), "mm"));
    SaveParameterInfo("BeamSigmaY", G4UIcommand::ConvertToString(fPrimary->GetSigmaY(), "mm"));
    SaveParameterInfo("BeamBetaX", G4UIcommand::ConvertToString(fPrimary->GetBetaX()));
    SaveParameterInfo("BeamBetaY", G4UIcommand::ConvertToString(fPrimary->GetBetaY()));
    SaveParameterInfo("BeamEmittanceX", G4UIcommand::ConvertToString(fPrimary->GetEmittanceX()));
    SaveParameterInfo("BeamEmittanceY", G4UIcommand::ConvertToString(fPrimary->GetEmittanceY()));
    SaveParameterInfo("Particle", pgun->GetParticleDefinition()->GetParticleName());
    SaveParameterInfo("Energy", G4UIcommand::ConvertToString(pgun->GetParticleEnergy(), "GeV"));

  } else if (fPrimary->GetSpectraSettings()) {
    SaveParameterInfo("BeamType", "Spectra");
    SaveParameterInfo("BeamPosition",
        G4UIcommand::ConvertToString(G4ThreeVector(fPrimary->GetX0(), fPrimary->GetY0(), fPrimary->GetZ0()), "mm"));
    SaveParameterInfo("Particle", pgun->GetParticleDefinition()->GetParticleName());
  }

  SaveParameterInfo("RunId", G4UIcommand::ConvertToString(fRun->GetRunID()));
#ifdef LX_TRACK_TIMING_CUT
  SaveParameterInfo("TimingCut", G4UIcommand::ConvertToString(LX_TRACK_TIMING_CUT, "ns"));
#endif
#ifdef LX_FAST_SIM
  SaveParameterInfo("FastSim", "ON");
#endif

}



void RunAction::SaveParameterInfo(const G4String paramid, const G4String paramval)
{
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  analysisManager->FillNtupleSColumn(6, 0, paramid);
  analysisManager->FillNtupleSColumn(6, 1, paramval);
  analysisManager->AddNtupleRow(6);
}


void RunAction::PrintRunInfo()
{
#ifdef LX_TRACK_TIMING_CUT
  G4cout << "******************** LX_TRACK_TIMING_CUT ********************" << G4endl;
  G4cout << "This run uses timing cut for particle simulation: " << LX_TRACK_TIMING_CUT << G4endl;
  G4cout << "************************************************************" << G4endl;
#endif
#ifdef LX_FAST_SIM
  G4cout << "******************** LX_FAST_SIM ********************" << G4endl;
  G4cout << "This run uses fast simulation option!" << G4endl;
  G4cout << "Particles simulation wil stop when they hit the following volumes:" << G4endl;
  G4cout << "BeamDumpAssembly, Shielding, GammaBeamDumpAssembly, HICSDumpAssembly, HICSDump2T, " << G4endl;
  G4cout << "ComptShieldingConcrete, ComptShieldingPlate" << G4endl;
  G4cout << "************************************************************" << G4endl;
#endif

}
