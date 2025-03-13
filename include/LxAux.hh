//
/// \brief LxAux class
//

#ifndef LxAux_h
#define LxAux_h 1

#include <map>
#include <tuple>

#include "G4UserRunAction.hh"
#include "globals.hh"


class G4AssemblyVolume;
class G4LogicalVolume;


class LxAux
{
public:
  LxAux() {};
  ~LxAux() {};

  static G4LogicalVolume *BuildPedestal(const G4String pname, const G4double xs, const G4double ys, const G4double zs);
  static void AddAssmblyVolumes(G4AssemblyVolume* vasm, G4AssemblyVolume* cpyasm, G4ThreeVector &translation, G4RotationMatrix *rotation);

};


#endif

