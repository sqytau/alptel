
#ifndef STLConstruction_h
#define STLConstruction_h 1

#include <vector>
#include <string>

#include "LxDetector.hh"


class STLConstruction: public LxDetector
{
  public:
    STLConstruction(DetectorConstruction *detc = 0);
    virtual ~STLConstruction() {};
    virtual void Construct();

  protected:
    std::vector<std::string> fFileList;
};


#endif
