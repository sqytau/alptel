
#ifndef LxDetector_h
#define LxDetector_h 1

class DetectorConstruction;
class G4AssemblyVolume;


class LxDetector
{
  public:
    LxDetector(DetectorConstruction *detc = 0): fDetector(detc) {};
    virtual ~LxDetector() {};
    virtual void Construct() = 0;
 
  protected:
    DetectorConstruction*  fDetector;
};


class WISDetectorTele: public LxDetector
{
  public:
    WISDetectorTele(DetectorConstruction *detc = 0): LxDetector(detc) {};
    virtual ~WISDetectorTele() {};
    virtual void Construct();

  protected:
    G4LogicalVolume* ConstructSensor();
    G4LogicalVolume* ConstructPCB();
    void AddSegmentation();
    void CreateMaterial();
    void ConstructShielding(G4LogicalVolume  *logicWorld);
    virtual void ConstructBottomSupport(G4LogicalVolume  *logicWorld);
    G4LogicalVolume* ConstructSr90Sourse();
};


class WISDetectorTeleFrame: public WISDetectorTele
{
  public:
    WISDetectorTeleFrame(DetectorConstruction *detc = 0): WISDetectorTele(detc) {};
    virtual ~WISDetectorTeleFrame() {};
    virtual void Construct();

  protected:
    G4LogicalVolume* ConstructROPCB();
    G4LogicalVolume* ConstructAlFrame();
//     void ConstructShielding(G4LogicalVolume  *logicWorld);
    virtual void ConstructBottomSupport(G4LogicalVolume  *logicWorld);
    G4AssemblyVolume* ConstructSupportAssembly();
};

#endif
