
#ifndef LxBSMDetector_h
#define LxBSMDetector_h 1


class LxDetector;

class LxBSMDetector: public LxDetector
{
  public:
    LxBSMDetector(DetectorConstruction *detc = 0): LxDetector(detc) {};
    virtual ~LxBSMDetector() {};
    virtual void Construct();

  protected:
    void CreateMaterial();
    void AddSegmentation();
    G4AssemblyVolume* ConstructSupportAssembly();
};



#endif
