
#ifndef LxINteratcionChambers_h
#define LxINteratcionChambers_h 1


class InteractionChamberAssembly
{
  public:
    InteractionChamberAssembly() : fIPChamberType(""), fIPChamberAssembly(0) {};
    virtual ~InteractionChamberAssembly() {};
    virtual G4AssemblyVolume* GetAssembly() { return fIPChamberAssembly; };

  protected: 
    InteractionChamberAssembly(const G4String mtypename) : fIPChamberType(mtypename), fIPChamberAssembly(0) {};
      
  protected: 
    G4String           fIPChamberType;  
    G4AssemblyVolume  *fIPChamberAssembly;
};



class SimpleIPBox: public InteractionChamberAssembly
{
  public:
    SimpleIPBox(const G4String mtypename);
    virtual ~SimpleIPBox() { if (fIPChamberAssembly && !fIPChamberAssembly->GetImprintsCount() ) 
                                        { delete fIPChamberAssembly; }  };
  protected: 
    void ConstructInteractionChamber();
};



class CircularChamber: public InteractionChamberAssembly
{
  public:
    CircularChamber(const G4String mtypename);
    virtual ~CircularChamber() { if (fIPChamberAssembly && !fIPChamberAssembly->GetImprintsCount() ) 
                                        { delete fIPChamberAssembly; }  };
  protected: 
    void ConstructInteractionChamber();
};



class LxInteractionChamber: public InteractionChamberAssembly
{
  public:
    LxInteractionChamber(const G4String mtypename);
    virtual ~LxInteractionChamber() { if (fIPChamberAssembly && !fIPChamberAssembly->GetImprintsCount() ) 
                                        { delete fIPChamberAssembly; }  };
  protected: 
    void ConstructInteractionChamber();
    void ConstructMirrors(G4LogicalVolume *lgTAUICContainer);
    G4AssemblyVolume* ConstructMirrorAssembly(const G4double r, const G4double d, const G4double h, const G4String mname);
};


#endif
