
#ifndef LxECal_h
#define LxECal_h 1

class LxDetector;


class LxECal: public LxDetector
{
  public:
    LxECal(DetectorConstruction *detc = 0): LxDetector(detc), fECalLayerZ(0.0) {};
    virtual ~LxECal() {};
    virtual void Construct();

  protected:
    void CreateMaterial();
    void AddSegmentation();
    G4AssemblyVolume* ConstructSupportAssembly();
    G4AssemblyVolume* ConstructCasingAssembly(const G4double ecalz);
    void ConstructShielding();
    void ConstructDumpShielding();
    G4LogicalVolume* ConstructPCB(const G4double ecalx, const G4double ecalz);

  protected:
    G4double fECalLayerZ;
};


#endif
