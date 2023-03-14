
#ifndef LxVacChambers_h
#define LxVacChambers_h 1


class LxVacChamberAssembly
{
  public:
    LxVacChamberAssembly() {};
    virtual ~LxVacChamberAssembly() {};
    virtual G4AssemblyVolume* GetAssembly() = 0;
};



class PipeChamberAssembly: public LxVacChamberAssembly
{
  public:
    PipeChamberAssembly(const G4String mtypename);
    virtual ~PipeChamberAssembly() { if (fChamberAssembly && !fChamberAssembly->GetImprintsCount() ) 
                                        { delete fChamberAssembly; }  };
    virtual G4AssemblyVolume* GetAssembly() { return fChamberAssembly; };

  protected: 
    void ConstructVacuumChamber();
      
  protected: 
    G4String           fChamberType;  
    G4AssemblyVolume  *fChamberAssembly;
};



class DetChamberAssembly: public LxVacChamberAssembly
{
  public:
    DetChamberAssembly(const G4String mtypename);
    virtual ~DetChamberAssembly() { if (fChamberAssembly && !fChamberAssembly->GetImprintsCount() ) 
                                        { delete fChamberAssembly; }  };
    virtual G4AssemblyVolume* GetAssembly() { return fChamberAssembly; };

  protected: 
    void ConstructVacuumChamber();
    G4AssemblyVolume*  ConstructFrameAssembly();
      
  protected: 
    G4String           fChamberType;  
    G4AssemblyVolume  *fChamberAssembly;

    G4double  fVCLength;
};


#endif
