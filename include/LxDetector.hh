
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



class LxDetectorOPPP: public LxDetector
{
  public:
    LxDetectorOPPP(DetectorConstruction *detc = 0): LxDetector(detc) {};
    virtual ~LxDetectorOPPP() {};
    virtual void Construct();

  protected:
    void AddSegmentation();
};


class LxDetectorCompton: public LxDetector
{
  public:
    LxDetectorCompton(DetectorConstruction *detc = 0): LxDetector(detc) {};
    virtual ~LxDetectorCompton() {};
    virtual void Construct();

  protected:
    void AddSegmentation();
};

class LxDetectorComptonFluka: public LxDetector
{
  public:
    LxDetectorComptonFluka(DetectorConstruction *detc = 0): LxDetector(detc) {};
    virtual ~LxDetectorComptonFluka() {};
    virtual void Construct();

  protected:
    void AddSegmentation();
    G4AssemblyVolume* ConstructSupportAssembly(G4double &sphight);
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
    void ConstructBottomSupport(G4LogicalVolume  *logicWorld);
    void ConstructSr90Sourse(G4LogicalVolume  *logicWorld);
//     G4AssemblyVolume* ConstructSupportAssembly(G4double &sphight);
};


#endif
