
#ifndef LxHICSBeam_h
#define LxHICSBeam_h 1

class LxDetector;


class LxHICSBeam: public LxDetector
{
  public:
    LxHICSBeam(DetectorConstruction *detc = 0): LxDetector(detc), fHICSDumpAngle(0.0) {};
    virtual ~LxHICSBeam() {};
    virtual void Construct();

  protected:
    void AddSegmentation();
    void ConstructElectronShielding(G4LogicalVolume  *lgicWorld);
    void ConstructNeutronAbsorber(G4LogicalVolume  *logicWorld);
    G4AssemblyVolume* ConstructSupportAssembly(G4double &supporthight);
    G4AssemblyVolume* ConstructHICSElDetSupportAssembly(G4double &sphight);
    void CreateMaterial();
    
  private:
    G4double fHICSDumpAngle;  
};


#endif
