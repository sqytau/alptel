//
/// \brief Definition of the DetectorConstruction class
//

#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include <map>
#include <tuple>

#include "G4VUserDetectorConstruction.hh"
#include "G4LogicalVolume.hh"
#include "G4AssemblyVolume.hh"
#include "globals.hh"
#include "G4Cache.hh"

#include "EventAction.hh"

class G4Box;
class G4VPhysicalVolume;
class G4Material;
class G4MaterialCutsCouple;
class G4UniformMagField;
class DetectorMessenger;
class G4GlobalMagFieldMessenger;
class G4UserLimits;
class LxMagnetAssembly;
class LxDetector;
class LxTargetChamber;
class G4MagneticField;
class LxBField;
class FieldDistribution;


class DetectorConstruction : public G4VUserDetectorConstruction
{
  public:

    DetectorConstruction();
   ~DetectorConstruction();

  public:

     void SetAbsorberMaterial (G4String);
     void SetAbsorberThickness(G4double);
     void SetAbsorberSizeXY   (G4double);
     void SetAbsorberType(G4String val);
     void SetGammaBeamDumpLength(G4double val);
     void SetGammaBeamDumpMaterial(G4String val);

     void SetAbsorberZpos(G4double);

     void SetWorldMaterial(G4String);
     void SetWorldSizeZ   (G4double);
     void SetWorldSizeXY  (G4double);

     void SetMagField(G4double);

     virtual G4VPhysicalVolume* Construct();
     virtual void ConstructSDandField();

  public:

     void PrintCalorParameters();

     G4Material* GetAbsorberMaterial()  const {return fAbsorberMaterial;}
     G4double    GetAbsorberThickness() const {return fAbsorberThickness;}
     G4double    GetAbsorberSizeXY()    const {return fAbsorberSizeXY;}

     G4double    GetAbsorberZpos()    const  {return fZposAbs;}
     G4double    GetzstartAbs()       const  {return fZstartAbs;}
     G4double    GetzendAbs()         const  {return fZendAbs;}

     G4Material* GetWorldMaterial()     {return fWorldMaterial;}
     G4double    GetWorldSizeZ()      const  {return fWorldSizeZ;}
     G4double    GetWorldSizeXY()     const  {return fWorldSizeXY;}

     const G4VPhysicalVolume* GetphysiWorld() {return fPhysiWorld;}
     const G4VPhysicalVolume* GetAbsorber()   {return fPhysiAbsorber;}
     const G4MaterialCutsCouple* GetAbsorbMaterialCut()  const
                             {return fLogicAbsorber->GetMaterialCutsCouple();}

     G4double GetMagnetZend() { return fMagnetZPos + fMagnetSizeZ/2.0; }

     void AddSensorSegmentation(const G4String sname, const G4double xsize, const G4double ysize,
                                const G4int nx,const G4int ny);
     const std::tuple<G4double, G4double, G4int, G4int> &GetSensorSegmentation(const G4String sname)
                                                                        { return fSensors.at(sname); }
     void AddBFieldModel(const G4String magid, std::tuple<G4String, G4String, G4String, G4String> bdatmodel)
                    {fBFieldModelsInfo[magid].push_back(bdatmodel); };
     void DumpBFieldModel();

  private:

     G4Material*        fAbsorberMaterial;
     G4double           fAbsorberThickness;
     G4double           fAbsorberSizeXY;

     G4double           fZposAbs;
     G4double           fZstartAbs, fZendAbs;

     G4Material*        fWorldMaterial;
     G4double           fWorldSizeZ;
     G4double           fWorldSizeXY;

     G4bool             fDefaultWorld;

     G4Box*             fSolidWorld;
     G4LogicalVolume*   fLogicWorld;
     G4VPhysicalVolume* fPhysiWorld;

     G4VSolid*          fSolidAbsorber;
     G4LogicalVolume*   fLogicAbsorber;
     G4VPhysicalVolume* fPhysiAbsorber;

     G4VSolid*          fGSolidAbsorber;
     G4LogicalVolume*   fGLogicAbsorber;
     G4VPhysicalVolume* fGPhysiAbsorber;

     DetectorMessenger* fDetectorMessenger;
     G4Cache<G4GlobalMagFieldMessenger*> fFieldMessenger;

     G4double           fMagnetZPos, fMagnetSizeZ;

     G4UserLimits       *fStepLimit;

     std::map<G4String, LxMagnetAssembly*> fMagentsTypeCollection;
     std::map<G4String, std::tuple<G4double, G4double, G4int, G4int> > fSensors;
     std::map<G4String, LxDetector*>  fDetList;

     LxTargetChamber   *fTargetChamber;

     std::map<G4String, std::vector<std::tuple<G4String, G4String, G4String, G4String> > > fBFieldModelsInfo;

  private:

     void DefineMaterials();
     void ComputeCalorParameters();
     G4VPhysicalVolume* ConstructCalorimeter();

     void ConstructIP();
     void ConstructIPChamber();
     void ConstructDumpMagnet(const G4String magType = "FlashMagnet");
     void ConstructBeamDump();
     void ConstructShielding();
     void ConstructBeamPipeToDump();
     void ConstructBeamPipeToIP();
     void ConstructLuxeDetectors();
     void ConstructGammaMagnet();
     void ConstructBremsTarget();
     void ConstructGammaTarget();
     void ConstructIPMagnet(const G4String magType = "FlashMagnet");
     void ConstructCleaningMagnet();
     void ConstructMagnet(G4LogicalVolume *mfContainer, const G4ThreeVector mf);
     void ConstructComptShielding();
     void ConstructInfrastructure();
     void ConstructWalls();

     void ConstructVacuumChamber();
     void AssignRegions();

     LxBField* ComposeFieldObject(const G4String magid, const G4ThreeVector magpos);
     FieldDistribution* CreateFieldDistribution(const G4String &fmodel, const G4String &params);
     void ConfigureFieldManager(G4FieldManager *fieldMgr, G4MagneticField *mfield);
     void AddFieldToLogVolumes(G4FieldManager* fieldMgr, const std::vector<G4String> vname);
};



#endif

