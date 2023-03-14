
#ifndef LxTrackerOPPP_h
#define LxTrackerOPPP_h 1


class LxDetector;

class LxTrackerOPPP: public LxDetector
{
  public:
    LxTrackerOPPP(DetectorConstruction *detc = 0);
    virtual ~LxTrackerOPPP() {};
    virtual void Construct();

  protected:
    void CreateMaterial();
    void AddSegmentation();
    G4AssemblyVolume* ConstructStaveAssembly();
    G4AssemblyVolume* ConstructStaveHolderAssembly();
    G4LogicalVolume*  ConstructSensorFlex();
    G4AssemblyVolume* ConstructSupportAssembly(G4double &sphight);
    G4AssemblyVolume* ConstructHexapodAssembly(G4double &hexhight);
    G4AssemblyVolume* ConstructServiceSupportAssembly(G4double &width);
    G4AssemblyVolume* ConstructCoolingPipeConnectorAssembly();
    G4AssemblyVolume* ConstructInnerServiceLinesAssembly();
    G4AssemblyVolume* ConstructOuterServiceLinesAssembly();
    G4AssemblyVolume* ConstructInnerServiceLinesTermAssembly();

    void  ConstructSideServiceLines(G4LogicalVolume  *worldVol);
    void  ConstructSideServiceLinesOutAssemblies(G4AssemblyVolume* &sideServiceAssembly,
                        G4AssemblyVolume* &sideFlatCableTermAssemblyL, G4AssemblyVolume* &sideFlatCableTermAssemblyR);

    void  ConstructSideServiceLinesInAssemblies(G4AssemblyVolume* &sideServiceAssembly,
                        G4AssemblyVolume* &sideFlatCableTermAssemblyL, G4AssemblyVolume* &sideFlatCableTermAssemblyR);

    G4AssemblyVolume* ConstructFlatCbleTerminator1Assembly();
    G4AssemblyVolume* ConstructFlatCbleTerminator2Assembly();

    void ConstructElectronicsRack();

  protected:
    G4AssemblyVolume* fCoolingPipeConnectorAssembly;
    G4AssemblyVolume* fCableTerm1Assembly;
    G4AssemblyVolume* fCableTerm2Assembly;

    // Position of the components in the assembly for detector container to join them with service lines container
    G4double fOutServFlatCblAssmblyZpos;
    G4double fInServCoolPipeAssmblyZpos;
    G4double fInServFlatCbl0AssmblyZpos;
    G4double fInServFlatCbl1AssmblyZpos;
    G4double fInServFlatCbl2AssmblyZpos;

};



#endif
