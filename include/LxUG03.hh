
#ifndef LxUG03_h
#define LxUG03_h 1


class LxDetector;

class LxUG03: public LxDetector
{
  public:
    LxUG03(DetectorConstruction *detc = 0): LxDetector(detc) {};
    virtual ~LxUG03() {};
    virtual void Construct();

  protected:
    void CreateMaterial();

};


#endif
