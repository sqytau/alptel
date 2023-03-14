
#ifndef GammaCalo_h
#define GammaCalo_h 1

#include "LxDetector.hh"

class LxDetectorGammaCalo: public LxDetector
{
  public:
    LxDetectorGammaCalo(DetectorConstruction *detc = 0): LxDetector(detc) {};
    virtual ~LxDetectorGammaCalo() {};
    virtual void Construct();

  private:
    void AddSegmentation();
};

class LxDetectorGammaCalo48: public LxDetector
{
  public:
    LxDetectorGammaCalo48(DetectorConstruction *detc = 0): LxDetector(detc) {};
    virtual ~LxDetectorGammaCalo48() {};
    virtual void Construct();

  protected:
    void CreateMaterial();
  private:
    void AddSegmentation();
};

class LxDetectorGammaCalo8: public LxDetectorGammaCalo48
{
  public:
    LxDetectorGammaCalo8(DetectorConstruction *detc = 0): LxDetectorGammaCalo48(detc) {};
    virtual ~LxDetectorGammaCalo8() {};
    virtual void Construct();

 private:
    void AddSegmentation();
    
};

#endif
