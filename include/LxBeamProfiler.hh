
#ifndef LxBeamProfiler_h
#define LxBeamProfiler_h 1


class LxDetector;

class LxBeamProfiler: public LxDetector
{
  public:
    LxBeamProfiler(DetectorConstruction *detc = 0): LxDetector(detc) {};
    virtual ~LxBeamProfiler() {};
    virtual void Construct();

  protected:
    void CreateMaterial();
    void AddSegmentation();
    G4AssemblyVolume* ConstructPCBAssembly();
    G4AssemblyVolume* ConstructPCBSupportAssembly();
    G4AssemblyVolume* ConstructMotorsAssembly();
    G4AssemblyVolume* ConstructSupportAssembly();
};



#endif
