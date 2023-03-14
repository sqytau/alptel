
#ifndef LxMagnets_h
#define LxMagnets_h 1

#include <vector>

#include "G4ThreeVector.hh"

class G4MagneticField;

class LxMagnetAssembly
{
  public:
    LxMagnetAssembly() {};
    virtual ~LxMagnetAssembly() {};
    virtual G4LogicalVolume* GetFieldVolume(const G4String lvname) = 0;
    virtual G4AssemblyVolume* GetAssembly() = 0;
    virtual void CostructSupport(G4LogicalVolume* lWorld, const G4ThreeVector& pos, const G4String& mname,
                                 G4bool rotate = false) = 0;

};



class TypMBMagnetAssembly: public LxMagnetAssembly
{
  public:
    TypMBMagnetAssembly(const G4String mtypename);
    virtual ~TypMBMagnetAssembly() {if (fMagnetAssembly && !fMagnetAssembly->GetImprintsCount() ) {delete fMagnetAssembly;}};
    virtual G4LogicalVolume* GetFieldVolume(const G4String lvname);
    virtual G4AssemblyVolume* GetAssembly() { if (!fMagnetAssembly) ConstructMagnet(); return fMagnetAssembly; };
    virtual void CostructSupport(G4LogicalVolume* lWorld, const G4ThreeVector& pos, const G4String& mname,
                                 G4bool rotate = false);

  protected:
    void ConstructMagnet();

  protected:
    G4String          fMagnetType;
    G4AssemblyVolume *fMagnetAssembly;
    G4MagneticField  *fField;
};



class PDSMagnetAssembly: public LxMagnetAssembly
{
  public:
    PDSMagnetAssembly(const G4String mtypename);
    virtual ~PDSMagnetAssembly() {if (fMagnetAssembly && !fMagnetAssembly->GetImprintsCount() ) {delete fMagnetAssembly;}};
    virtual G4LogicalVolume* GetFieldVolume(const G4String lvname);
    virtual G4AssemblyVolume* GetAssembly() { ConstructMagnet(); return fMagnetAssembly; }
    virtual void CostructSupport(G4LogicalVolume* lWorld, const G4ThreeVector& pos, const G4String& mname,
                                 G4bool rotate = false);

  protected:
    void ConstructMagnet();

  protected:
    G4String          fMagnetType;
    G4AssemblyVolume *fMagnetAssembly;
};



class FlashMagnetAssembly: public LxMagnetAssembly
{
  public:
    FlashMagnetAssembly(const G4String mtypename);
    virtual ~FlashMagnetAssembly() {if (fMagnetAssembly && !fMagnetAssembly->GetImprintsCount() ) {delete fMagnetAssembly;}};
    virtual G4AssemblyVolume* GetAssembly() { ConstructMagnet(); return fMagnetAssembly; };
    virtual G4LogicalVolume* GetFieldVolume(const G4String lvname);
    virtual void CostructSupport(G4LogicalVolume* lWorld, const G4ThreeVector& pos, const G4String& mname,
                                 G4bool rotate = false);

  protected:
    void ConstructMagnet();

  protected:
    G4String          fMagnetType;
    G4AssemblyVolume *fMagnetAssembly;
};



class PDSMagnetField: public G4MagneticField
{
  public:
    PDSMagnetField(const G4ThreeVector pos);

    virtual void GetFieldValue(const G4double ppos[4], G4double *MagField) const;

  protected:
    void LoadField();
    G4double FField(const G4double x, const G4int ii) const;
    void LoadFieldData(const G4String fname);

  template <class T> class point_cmp
  {
    public:
    point_cmp(const std::vector<T> *ptr): p(ptr) {};
    bool operator() (size_t i, size_t j) const { if ((i < p->size()) && (j < p->size()) ) return (p->at(i) < p->at(j));
                                                 return false; }
    const std::vector<T> *p;
  };

  G4double FindZ(const G4ThreeVector &r1, const G4ThreeVector &r2, const G4ThreeVector &r3,
                               const double xx, const double yy) const;
  G4bool CheckCollinear(const G4ThreeVector &r1, const G4ThreeVector &r2, const G4ThreeVector &r3) const;
  G4double DataField(const G4double xx, const G4double yy) const;

  protected:
    G4ThreeVector fPosition;
    G4ThreeVector fMin, fMax;
    G4ThreeVector fBField;

    std::map<int, std::vector<double> > fInvExpParams;
    std::vector<G4ThreeVector> fFieldData;
};



class TypMBMagnetSimple: public TypMBMagnetAssembly
{
  public:
    TypMBMagnetSimple(const G4String mtypename);
    virtual ~TypMBMagnetSimple() {if (fMagnetAssembly && !fMagnetAssembly->GetImprintsCount() ) {delete fMagnetAssembly;}};
    virtual G4AssemblyVolume* GetAssembly() { ConstructMagnet(); return fMagnetAssembly; };

  protected:
    void ConstructMagnet();
};


#endif
