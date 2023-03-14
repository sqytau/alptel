//
/// \brief Definition of the RunAction class
//

#ifndef RunAction_h
#define RunAction_h 1

#include <map>
#include <tuple>
#include <vector>

#include "G4UserRunAction.hh"
#include "globals.hh"



class Run;
class DetectorConstruction;
class PrimaryGeneratorAction;
class HistoManager;
class RunActionMessenger;



class RunAction : public G4UserRunAction
{

public:

    RunAction(DetectorConstruction* det, PrimaryGeneratorAction* prim=0);
   ~RunAction();
   
    virtual G4Run* GenerateRun();    
    virtual void BeginOfRunAction(const G4Run*);
    virtual void   EndOfRunAction(const G4Run*);

    void AddInterceptVolume (const G4String volumename);
    void AddInterceptVolumeECut (const G4String volname, const G4double ecutt);
    void AddSensitiveVolume (const G4String volumename);

    void SetDumpGeometry(const G4bool ff) { fDumpGeometry = ff; }
    void SetSaveTrajectory(const G4bool ff) { fSaveTrajectory = ff; }
    void SetExternalRunID(const G4int rr) { fExternalRunID = rr; }

  protected:
    void SaveRunInfo();
    void SaveParameterInfo(const G4String paramid, const G4String paramval);
    void PrintRunInfo();

  private:
    DetectorConstruction*   fDetector;
    PrimaryGeneratorAction* fPrimary;
    Run*                    fRun;        
    HistoManager*           fHistoManager;

    RunActionMessenger        *fRunActionMessenger;
    std::map<G4String, std::vector<G4int> >  fInterceptVol;
    std::map<G4int, G4double>  fInterceptVolECut;
    std::map<G4String, std::tuple<G4int, G4int, G4int> > fSensitiveVol;
    
    G4bool  fDumpGeometry;
    G4bool  fSaveTrajectory;
    G4int   fExternalRunID;
    
};



#endif

