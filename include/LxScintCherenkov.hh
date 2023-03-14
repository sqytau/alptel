
#ifndef LxScintCherenkov_h
#define LxScintCherenkov_h 1


class LxDetector;

class LxSScreen: public LxDetector
{
  public:
    LxSScreen(DetectorConstruction *detc = 0): LxDetector(detc) {};
    virtual ~LxSScreen() {};
    virtual void Construct();

  protected:
    void ConstructSupportAssembly();
};



class LxCherenkov: public LxDetector
{
  public:
    LxCherenkov(DetectorConstruction *detc = 0): LxDetector(detc) {};
    virtual ~LxCherenkov() {};
    virtual void Construct();

  protected:
    void AddSegmentation();
//     G4AssemblyVolume* ConstructSupportAssembly(G4double &sphight);
//     G4AssemblyVolume* ConstructTableAssembly(const G4double tabley);
};



#endif
