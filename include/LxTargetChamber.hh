
#ifndef LxTargetChamber_h
#define LxTargetChamber_h 1


class LxTargetChamber
{
  public:
    LxTargetChamber(): flogicTargetChamber(0) {};
    ~LxTargetChamber() {};

    G4LogicalVolume* GetTargetChamber() {if (!flogicTargetChamber) ConstructTargetChamber(); return flogicTargetChamber;}
    G4LogicalVolume* GetTargetVolume(const G4String volname);
    void ConstructTargetChamberSupport(G4LogicalVolume *logicWorld, const G4String tcname, const G4double zpos);

  protected:
    void ConstructTargetChamber();
    G4AssemblyVolume* ConstructTargetSupportAssembly();
//     void CreateMaterial();

  private:
    G4LogicalVolume* flogicTargetChamber;
};


#endif
