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

  static G4AssemblyVolume *BuildTable(const G4String tname, const G4double xs, const G4double ys, 
                 const G4double zs, const int nleg, const G4double ytop = -1.0, const G4double rleg = -1.0);
  static G4LogicalVolume *BuildPedestal(const G4String pname, const G4double xs, const G4double ys, const G4double zs);
  static G4AssemblyVolume *BuildHexapod(const G4String pname, G4double &hexhight);
  static G4AssemblyVolume *BuildConeHexapod(const G4String pname, G4double &hexhight);
  static void AddAssmblyVolumes(G4AssemblyVolume* vasm, G4AssemblyVolume* cpyasm, G4ThreeVector &translation, G4RotationMatrix *rotation);

};


#endif

