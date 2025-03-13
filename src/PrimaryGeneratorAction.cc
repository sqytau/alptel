//
/// \brief Implementation of the PrimaryGeneratorAction class
//

#include <algorithm>
#include <tuple>
#include <string>

#include "PrimaryGeneratorAction.hh"

#include "DetectorConstruction.hh"
#include "PrimaryGeneratorMessenger.hh"

#include "G4Event.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "G4UnitsTable.hh"
#include "PrimarySpectra.hh"
#include "LuxeTestGenerator.hh"
#include "G4RunManager.hh"
#include "LXSetUp.hh"

#include "EventAction.hh"

#include "G4AutoLock.hh"

namespace { G4Mutex MCReadMutex = G4MUTEX_INITIALIZER; }

LuxeMCGenerator *PrimaryGeneratorAction::lxgen = 0;


PrimaryGeneratorAction::PrimaryGeneratorAction(DetectorConstruction* DC)
:G4VUserPrimaryGeneratorAction(),
 fParticleGun(0),fDetector(DC),fBeamType(beamGauss),
 fx0(0.0), fy0(0.0), fz0(0.0), fsigmax(5.0*um), fsigmay(5.0*um), fsigmaz(24.0*um),
 femittancex(0.0), femittancey(0.0), fbetax(0.0), fbetay(0.0), fbeamfocus(0.0), fSpectra(0),
 fMCfile(""), fMCList(false), fnfixparticles(0), fMCWeightScale(1), fMCEventInfo(0),
 fselectpdgvec(0), fselectpdg(false), fSkipEvents(0), fMCPrimeParentId(-1), fMCPrimeNDescendant(-1),
 fMCPrimeXi(-1.0), fMCPhysProc(0), fGunMessenger(0)
{
  G4int n_particle = 1;
  fParticleGun  = new G4ParticleGun(n_particle);
  SetDefaultKinematic();
  fMCEventInfo = new EventInfo(1.0, -1); // set event weight

  //create a messenger for this class
  fGunMessenger = new PrimaryGeneratorMessenger(this);
}



PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete fParticleGun;
  delete fGunMessenger;
  if (fSpectra) delete fSpectra;
  G4AutoLock mcfilemutex(&MCReadMutex);
  if (lxgen) {
    delete lxgen;
    lxgen = 0;
  }
  mcfilemutex.unlock();
  if (fMCEventInfo) delete fMCEventInfo;
}



void PrimaryGeneratorAction::SetDefaultKinematic()
{
  // default particle kinematic
  // defined by the beam parameters: emittancex, emittancey, fsigmax, fsigmay, fsigmaz and drift distance to IP fz0.
  //
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* particle = particleTable->FindParticle("e-");
  fParticleGun->SetParticleDefinition(particle);
  fParticleGun->SetParticleEnergy(16.5*GeV);

  fz0 = fDetector->GetzstartAbs() - 20.0*cm;

  fsigmax = 5.0*um;
  fsigmay = 5.0*um;
  fsigmaz = 24.0*um;

  G4double lf = fParticleGun->GetParticleEnergy() / particle->GetPDGMass();

  femittancex = 1.4e-3 * mm / lf;
  femittancey = 1.4e-3 * mm / lf;

  fbetax = fsigmax*fsigmax/femittancex;
  fbetay = fsigmay*fsigmay/femittancey;

  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));
  fParticleGun->SetParticlePosition(G4ThreeVector(fx0, fy0, fz0));
}



void PrimaryGeneratorAction::SetBeamType(G4String val)
{
  if (val == "gaussian") {
    fBeamType = beamGauss;
  } else if (val == "mono") {
    fBeamType = beamMono;
  } else if (val == "mc") {
    fBeamType = beamMC;
  } else if (val == "mchdf5") {
    fBeamType = beamMCh5;
  } else if (val == "mctupleg4") {
    fBeamType = beamMCTupleG4;
  } else if (val == "source") {
    fBeamType = beamSource;
  } else if (val == "cosmic") {
    fBeamType = beamCosmic;
  } else {
    G4cout << "PrimaryGeneratorAction::SetBeamType: <" << val << ">"
           << " is not defined. Using default: gaussian."
           << G4endl;
  }
}



void PrimaryGeneratorAction::SetSigmaX(const G4double sigma)
{
  fsigmax = sigma;
  fbetax = fsigmax*fsigmax/femittancex;
  std::cout << "PrimaryGeneratorAction: Set beam sigmaX at IP to : " << G4BestUnit(fsigmax, "Length") << std::endl;
//G4cout does not print anything...   G4cout << "Set beam sigmaX at IP to : " << G4BestUnit(fsigmax, "Length") << G4endl;
}



void PrimaryGeneratorAction::SetSigmaY(const G4double sigma)
{
  fsigmay = sigma;
  fbetay = fsigmay*fsigmay/femittancey;
  std::cout << "PrimaryGeneratorAction: Set beam sigmaY at IP to : " << G4BestUnit(fsigmay, "Length") << std::endl;
//   G4cout << "Set beam sigmaY at IP to : " << G4BestUnit(fsigmay, "Length") << G4endl;
}



void PrimaryGeneratorAction::SetBeamPosZ(const G4double posz)
{
  fz0 = posz;
  std::cout << "PrimaryGeneratorAction: Set beam Z position : " << G4BestUnit(fz0, "Length") << std::endl;
//   G4cout << "PrimaryGeneratorAction: Set beam Z position : " << G4BestUnit(fz0, "Length") << G4endl;
}



void PrimaryGeneratorAction::SetSpectraFile(G4String val)
{
  if (fSpectra) delete fSpectra;
  fSpectra = new PrimarySpectra();
//  fSpectra->SetVerbose(1);
  int ndata = fSpectra->LoadData(val);
  fSpectra->SetScale(GeV);
  G4cout << "PrimaryGeneratorAction: Primary spectra with " << ndata
            << " data points was loaded from the file " << val << G4endl;
}



void PrimaryGeneratorAction::SetMCParticleFile(G4String val, const G4bool list)
{
  G4AutoLock mcfilemutex(&MCReadMutex);
  if (lxgen) { delete lxgen;  lxgen = 0; }
  mcfilemutex.unlock();

  fMCfile = val;
  fMCList = list;
  if (fMCList) G4cout << "File with the list of files with primary MC particles " << fMCfile << G4endl;
  else         G4cout << "File with primary MC particles " << fMCfile << G4endl;
}



void PrimaryGeneratorAction::AddMCSelectParticle(const G4int pdgid)
{
  if (std::find(fselectpdgvec.begin(), fselectpdgvec.end(), pdgid) == fselectpdgvec.end() ) {
    G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
    G4ParticleDefinition *particle = particleTable->FindParticle(pdgid);
    if (particle) {
      fselectpdgvec.push_back(pdgid);
      fselectpdg = true;
      G4cout << "Particle with PDG_ID = " << pdgid << " was added to the list of selected MC patricles."  << G4endl;
    } else {
      G4cout << "Attempt to add particle with PDG_ID = " << pdgid
             << " to the list of selected MC patricles failed. PDG_ID is not valid! Ignore."  << G4endl;
    }
  }
}



void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  // this function is called at the begining of event
  if (fSpectra) fParticleGun->SetParticleEnergy(fSpectra->GetRandom());

  if (fBeamType == beamGauss) {
    GenerateGaussian(anEvent);
  } else if (fBeamType == beamMono) {
    GenerateMono(anEvent);
  } else if (fBeamType == beamSource) {
    GenerateSource(anEvent);
  } else if (fBeamType == beamCosmic) {
    GenerateCosmic(anEvent);
  }
    else if (fBeamType == beamMC || fBeamType == beamMCh5 || fBeamType == beamMCTupleG4) {
    GeneratefromMC(anEvent);
  }
}



void PrimaryGeneratorAction::GenerateGaussian(G4Event* anEvent)
{
  G4double z0 = G4RandGauss::shoot(fz0, fsigmaz);
  G4double zdrift = z0 - fbeamfocus;  // This is needed to have correct drift distance for x, y distribution.

  G4double sigmax = fsigmax * sqrt(1.0 + pow(zdrift/fbetax, 2.0));
  G4double x0 = G4RandGauss::shoot(fx0, sigmax);
  G4double meandx = x0*zdrift / (zdrift*zdrift + fbetax*fbetax);
  G4double sigmadx = sqrt( femittancex*fbetax / (zdrift*zdrift + fbetax*fbetax) );
  G4double dx0 = G4RandGauss::shoot(meandx, sigmadx);

  G4double sigmay = fsigmay * sqrt(1.0 + pow(zdrift/fbetay, 2.0));
  G4double y0 = G4RandGauss::shoot(fy0, sigmay);
  G4double meandy = y0*zdrift / (zdrift*zdrift + fbetay*fbetay);
  G4double sigmady = sqrt( femittancey*fbetay / (zdrift*zdrift + fbetay*fbetay) );
  G4double dy0 = G4RandGauss::shoot(meandy, sigmady);

  G4double mass = fParticleGun->GetParticleDefinition()->GetPDGMass();
  G4double E = fParticleGun->GetParticleEnergy();

  G4double pz = sqrt( (E*E - mass*mass)/ (dx0*dx0 + dy0*dy0 + 1.0) );

  fParticleGun->SetParticlePosition(G4ThreeVector(x0, y0, z0));
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(dx0*pz, dy0*pz, pz));
  fParticleGun->GeneratePrimaryVertex(anEvent);
}



void PrimaryGeneratorAction::GenerateSource(G4Event* anEvent)
{
  LXSetUp *lxs = LXSetUp::Instance();

  fz0 = -lxs->Sr90GapZ - lxs->BTargetZ;

  G4double rphi = 2.0*M_PI*G4UniformRand();
  G4double sreffectiver = lxs->Sr90SourceR;    // lxs->CollimatorLeadHoleR;
  G4double rr = sqrt(sreffectiver*sreffectiver*G4UniformRand());
  G4double x0 = rr*cos(rphi);
  G4double y0 = rr*sin(rphi);
  G4double z0 = fz0;   // fz0 + fz0+Sr90SourceZ*G4UniformRand();

  fParticleGun->SetParticlePosition(G4ThreeVector(x0, y0, z0));
  fParticleGun->SetParticleEnergy(0*eV);
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.0, 0.0, 1.0));

  G4int Z = 38, A = 90;
  G4double ionCharge   = 0.*eplus;
  G4double excitEnergy = 0.*keV;

  G4ParticleDefinition* ion = G4IonTable::GetIonTable()->GetIon(Z,A,excitEnergy);
  fParticleGun->SetParticleDefinition(ion);
  fParticleGun->SetParticleCharge(ionCharge);

  fParticleGun->GeneratePrimaryVertex(anEvent);
}



void PrimaryGeneratorAction::GenerateCosmic(G4Event* anEvent)
{
  LXSetUp *lxs = LXSetUp::Instance();

  G4int pid = 13;
  G4double z0 = lxs->TelescopeSensorZpos.at(0) - 1.0*mm;
  G4double tlx = lxs->OPPPSensorX;
  G4double tly = lxs->OPPPSensorY;
  G4double tlz = lxs->TelescopeSensorZpos.at(2) - z0;
  G4double vmax = 1 - pow(cos(atan2(sqrt(tlx*tlx+tly*tly),tlz)), 3.0);
  G4ThreeVector dir_p, posmu, rzl;
  do {
    posmu.set(tlx * (G4UniformRand() - 0.5), tly * (G4UniformRand() - 0.5), z0);
    dir_p.setRThetaPhi(1.0, acos(pow(1.0-vmax*G4UniformRand(), 1.0/3.0)), 2.0*M_PI*G4UniformRand());
    rzl = posmu + dir_p * tlz;
  } while ( (abs(rzl.x()) > tlx/2.0) || (abs(rzl.y()) > tly/2.0));

  fParticleGun->SetParticlePosition(posmu);
  fParticleGun->SetParticleMomentumDirection(dir_p);

  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* particle = particleTable->FindParticle(pid);
  fParticleGun->SetParticleDefinition(particle);

  fParticleGun->GeneratePrimaryVertex(anEvent);
}



void PrimaryGeneratorAction::GenerateMono(G4Event* anEvent)
{
  G4double mass = fParticleGun->GetParticleDefinition()->GetPDGMass();
  G4double E = fParticleGun->GetParticleEnergy();
//  G4double pz = sqrt(E*E - mass*mass);
  G4double pz = sqrt(E*E + 2.0*E*mass);  // E is kinetic energy

  fParticleGun->SetParticlePosition(G4ThreeVector(fx0, fy0, fz0));
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.0, 0.0, pz));
  fParticleGun->GeneratePrimaryVertex(anEvent);
}



void PrimaryGeneratorAction::InitReader()
{
  if (fBeamType == beamMC) {
    G4AutoLock mcfilemutex(&MCReadMutex);
    LuxeTestGenerator *gipr = new LuxeTestGenerator();
    if (fMCList) gipr->SetFileList(fMCfile);
    else gipr->AddEventFile(fMCfile);
    std::cout << "Processing file: " << fMCfile << std::endl;
//    gipr->SetFileType("out", 9, 9);
    gipr->SetFileType("out", 10, 10);
    lxgen = gipr;
    std::cout << "MC file type is set "  << std::endl;
    mcfilemutex.unlock();
    
  } else if (fBeamType == beamMCh5) {
    G4AutoLock mcfilemutex(&MCReadMutex);
    lxgen = new LxHDF5Reader(fMCfile);
    mcfilemutex.unlock();
    
  } else if (fBeamType == beamMCTupleG4) {
    G4AutoLock mcfilemutex(&MCReadMutex);
    lxgen = new LxNTupleReader(fMCfile);
    mcfilemutex.unlock();

  } else {
    G4String msgstr("Value of fBeamType: ");
    msgstr += std::to_string(fBeamType) + G4String(" is not supported\n");
    G4Exception("PrimaryGeneratorAction::", "InitReader()", FatalException, msgstr.c_str());
  }

  if (fSkipEvents > 0) {
    G4AutoLock mcfilemutex(&MCReadMutex);
    lxgen->SkipEvents(fSkipEvents);
    std::cout << "Skipping " << fSkipEvents << " events form the file " << fMCfile << std::endl;
    fSkipEvents = 0;   // This is to ensure that only one thread does it
    mcfilemutex.unlock();
  }
}



void PrimaryGeneratorAction::GeneratefromMC(G4Event* anEvent)
{
  if (!lxgen) {
    InitReader();
  }

  if (fnfixparticles > 0) {
    --fnfixparticles;
    fParticleGun->GeneratePrimaryVertex(anEvent);
    anEvent->SetUserInformation(new EventInfo(fMCEventInfo->GetWeight(), fMCEventInfo->GetMCTrackId()));
//     std::cout << "---------> New event with the same particle. Weight: " << fMCEventInfo->GetWeight() << std::endl;
    SavePrimaryTrack();
    return;
  }

  std::vector < std::vector <double> > ptcls;

  G4AutoLock mcfilemutex(&MCReadMutex);

  int nscat = lxgen->GetEventFromFile(ptcls);
  mcfilemutex.unlock();

//   std::cout << "GetEventFromFile "  << nscat << std::endl;
  if (nscat > 1) {
      G4String msgstr("Error reading particle from a file! More than one were read, it is not supported!\n");
      G4Exception("PrimaryGeneratorAction::", "GeneratefromMC(Event)", FatalException, msgstr.c_str());
  }
  if (nscat <= 0) {
//     mcfilemutex.lock();
//     delete lxgen;       // it will be destroied in destructor
//     lxgen = 0;          // and keeping it will prevent other threads from opening the same file again
//     mcfilemutex.unlock();
    if (nscat == -1) {
      G4cout << "All particle from the file " << fMCfile << " are processed." << G4endl;  
      fParticleGun->SetNumberOfParticles(0);
      ptcls.push_back(std::vector<double>(10, 90.0));
      G4RunManager::GetRunManager()->AbortRun(true);
//       G4RunManager::GetRunManager()->RunTermination();
    } else {
      G4String msgstr("Error reading MC particle from the file!\n");
      G4Exception("PrimaryGeneratorAction::", "GeneratefromMC(Event)", FatalException, msgstr.c_str());
    }
  }

  std::vector <double> pdata = ptcls[0];
  int pid = static_cast<int>(pdata[7]);

  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* particle = particleTable->FindParticle(pid);
  if (!particle) {
    if (pid == 90)  {
      particle = particleTable->FindParticle("geantino");
      G4String msgstr("Initial particle with pdg_id=");
      msgstr += std::to_string(pid) + G4String(" is in the file! Ignore it\n");
      G4Exception("PrimaryGeneratorAction::", "GeneratefromMC(Event)", JustWarning, msgstr.c_str());
    }
    else {
      G4String msgstr("Error setting initial particle from a file!\n");
      G4Exception("PrimaryGeneratorAction::", "GeneratefromMC(Event)", FatalException, msgstr.c_str());
    }
  }

  if (pdata.size()>12) {
    fMCPrimeParentId = pdata[10];
    fMCPrimeNDescendant = pdata[11];
    fMCPrimeXi = pdata[12];
  }
  if (pdata.size() > 14)  {
    fMCPhysProc = pdata[14];
  }

  std::vector <G4double> pp(4);
  G4double  vtx[3];
  G4double  wghtf, wght = 1.0;

  wghtf = pdata[8];
  fnfixparticles = 1;
  if (wghtf >= fMCWeightScale)   fnfixparticles = fMCWeightScale;
  else if (wghtf > 1.0)  fnfixparticles = static_cast<G4int>(floor(wghtf));
  wght = wghtf/static_cast<G4double>(fnfixparticles);
  fMCEventInfo->SetWeight(wght);
  fMCEventInfo->SetMCTrackId(pdata[9]);

  std::copy(pdata.begin()+4, pdata.begin()+7, pp.begin());
  std::copy(pdata.begin()+1, pdata.begin()+4, vtx);
  std::for_each(vtx, vtx+3, [=](G4double &x){x *= um;});
//  vtx[2] -= 7.0*m;//this is to get to IP
  pp[3] = pdata[0];

  fParticleGun->SetParticleDefinition(particle);
  fParticleGun->SetParticleEnergy(pp[3]*GeV - particle->GetPDGMass());
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(pp[0], pp[1], pp[2]));
  fParticleGun->SetParticlePosition(G4ThreeVector(vtx[0], vtx[1], vtx[2]));
  if (pdata.size() > 13)  fParticleGun->SetParticleTime(pdata[13]);

  if (    (!fselectpdg && (pid != 90))
       || (fselectpdg && std::find(fselectpdgvec.begin(), fselectpdgvec.end(), pid) != fselectpdgvec.end()) ) {
    if (true) {
      fParticleGun->SetNumberOfParticles(1);
      fnfixparticles -= 1;
//       std::cout << "Generating " << fnfixparticles+1 << " particles (weight " << wghtf
//                 << ", scale: " << wghtf/wght << ") with weight " << fMCEventInfo->GetWeight()
//                 << "  and following mC data:\n";
      // Save MC primary track
      SavePrimaryTrack();
    } else {
//       std::cout << "Particle does not hit the target, " << wghtf << " particle with following mC data:\n";
      fParticleGun->SetNumberOfParticles(0);
      fnfixparticles = 0;
    }
//     std::for_each(pdata.begin(), pdata.end(), [](const double x){std::cout << x << "  ";});
//     std::cout << std::endl;
  } else {
    fParticleGun->SetNumberOfParticles(0);
    fnfixparticles = 0;
  }

  fParticleGun->GeneratePrimaryVertex(anEvent);
  anEvent->SetUserInformation(new EventInfo(fMCEventInfo->GetWeight(), fMCEventInfo->GetMCTrackId()));
}



void PrimaryGeneratorAction::SavePrimaryTrack() const
{
  EventAction *evaction = dynamic_cast<EventAction*>(G4EventManager::GetEventManager()->GetUserEventAction());
  if (evaction) {
    G4int pid = fParticleGun->GetParticleDefinition()->GetPDGEncoding();
    G4double enegk = fParticleGun->GetParticleEnergy();
    G4ThreeVector momentum(fParticleGun->GetParticleMomentumDirection());
    G4double pmass = fParticleGun->GetParticleDefinition()->GetPDGMass();
    momentum *= std::sqrt(enegk*enegk + 2.0*enegk*pmass);
    G4ThreeVector vtx =	fParticleGun->GetParticlePosition();
    G4double t = fParticleGun->GetParticleTime();

    evaction->AddVolumeTrack(std::make_tuple(enegk/GeV, vtx.x(), vtx.y(), vtx.z(), t, vtx.x(), vtx.y(), vtx.z(),
              momentum.x()/GeV, momentum.y()/GeV, momentum.z()/GeV, momentum.getTheta(), momentum.getPhi(),
              0.0, 0.0, 0.0, fMCEventInfo->GetMCTrackId(), -1, pid, fMCPhysProc, fMCPrimeParentId, fMCPrimeNDescendant, fMCPrimeXi));

    //Copy primary particle to the HitTrack tree
    G4DynamicParticle *primeParticle = new G4DynamicParticle(fParticleGun->GetParticleDefinition(),
                                                             fParticleGun->GetParticleMomentumDirection(), fParticleGun->GetParticleEnergy());
    G4Track *trk = new G4Track(primeParticle, t, vtx);
    trk->SetTrackID(0);
    trk->SetParentID(fMCEventInfo->GetMCTrackId());
    trk->SetVertexPosition(vtx);

    G4StepPoint *prstp = new G4StepPoint();
    prstp->SetPosition(vtx);
    prstp->SetGlobalTime(t);
    prstp->SetMass(pmass);
    prstp->SetKineticEnergy(enegk);
    prstp->SetMomentumDirection(momentum.unit());
//     G4StepPoint *psstp = new G4StepPoint(*prstp);
    G4StepPoint *psstp = new G4StepPoint();
    psstp->SetPosition(vtx);
    psstp->SetGlobalTime(t);
    psstp->SetMass(pmass);
    psstp->SetKineticEnergy(enegk);
    psstp->SetMomentumDirection(momentum.unit());

    G4Step *trkstp = new G4Step();
    trkstp->SetTrack(trk);
    trkstp->SetTotalEnergyDeposit(0.0);
    trkstp->SetPreStepPoint(prstp);
    trkstp->SetPostStepPoint(psstp);

    trk->SetStep(trkstp);

    evaction->AddOpppTrackerCellDepposit(0, 0, 0, 0, fMCPrimeXi, trk);

    delete trk;     // deletes also primeParticle
    delete trkstp;  // deletes also PreStepPoint and PostStepPoint
  }
}
