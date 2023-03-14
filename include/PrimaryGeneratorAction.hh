//
/// \brief Definition of the PrimaryGeneratorAction class
//


#ifndef PrimaryGeneratorAction_h
#define PrimaryGeneratorAction_h 1

#include <vector>

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"
#include "globals.hh"

class G4Event;
class DetectorConstruction;
class PrimaryGeneratorMessenger;
class PrimarySpectra;
class LuxeMCGenerator;
class EventInfo;

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
    PrimaryGeneratorAction(DetectorConstruction*);
   ~PrimaryGeneratorAction();

  public:
    void SetDefaultKinematic();
    void SetBeamType(G4String val);
    void SetSpectraFile(G4String val);
    void SetMCParticleFile(G4String val, const G4bool list = false);
    virtual void GeneratePrimaries(G4Event*);
    const G4ParticleGun* GetParticleGun() const {return fParticleGun;}

    void GenerateGaussian(G4Event* anEvent);
    void GenerateMono(G4Event* anEvent);
    void GeneratefromMC(G4Event* anEvent);
    void GenerateSource(G4Event* anEvent);

    G4double GetX0() const {return fx0;}
    G4double GetY0() const {return fy0;}
    G4double GetZ0() const {return fz0;}
    G4double GetSigmaX() const {return fsigmax;}
    G4double GetBetaX()  const {return fbetax;}
    G4double GetEmittanceX() const {return femittancex;}
    G4double GetBeamFocus() const {return fbeamfocus;}

    G4double GetSigmaY() const {return fsigmay;}
    G4double GetBetaY()  const {return fbetay;}
    G4double GetEmittanceY() const {return femittancey;}
    G4int GetBeamType() const {return fBeamType;}
    G4String GetMCfile() const {return fMCfile;}
    const PrimarySpectra *GetSpectraSettings() const {return fSpectra;}

    void SetSigmaX(const G4double sigma);
    void SetSigmaY(const G4double sigma);
    void SetBeamPosZ(const G4double posz);
    void SetBeamFocus(const G4double zfocus) { fbeamfocus = zfocus; }
    void SetMCWeightScale(const G4int scale) { fMCWeightScale = scale; }
    void SetBeamPosition (const G4ThreeVector pos) { fx0 = pos.x(); fy0 = pos.y(); fz0 = pos.z(); }
    void SetSkipEvents(const size_t nev) {fSkipEvents = nev; }

    G4double GetSigmaZ() const {return fsigmaz;}
    size_t GetSkipEvents() const {return fSkipEvents;}

    void AddMCSelectParticle(const G4int pdgid);
    void SavePrimaryTrack() const;

    enum tBeamType : G4int {beamGauss, beamMono, beamMC, beamMCh5, beamMCTupleG4, beamSource};

  protected:

    G4double TestHitTarget(const std::vector <double> &pp, const double *vtx);
    void InitReader();

  private:
    G4ParticleGun*         fParticleGun;
    DetectorConstruction*  fDetector;
    G4int                  fBeamType;
    G4double               fx0, fy0, fz0;
    G4double               fsigmax, fsigmay, fsigmaz;
    G4double               femittancex, femittancey;
    G4double               fbetax, fbetay;
    G4double               fbeamfocus;

    static LuxeMCGenerator   *lxgen;
    PrimarySpectra        *fSpectra;
    G4String               fMCfile;
    G4bool                 fMCList;
    G4int                  fnfixparticles;

    G4int                  fMCWeightScale;
    EventInfo             *fMCEventInfo;

    std::vector<G4int>     fselectpdgvec;
    G4bool                 fselectpdg;

    size_t                 fSkipEvents;

    G4int                  fMCPrimeParentId;
    G4int                  fMCPrimeNDescendant;
    G4double               fMCPrimeXi;
    G4int                  fMCPhysProc;

    PrimaryGeneratorMessenger* fGunMessenger;
};



#endif


