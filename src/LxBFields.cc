//
/// \brief Implementation of the Magnetic field models classes
//

#include <algorithm>
#include <functional>
#include <fstream>

#include "G4UnitsTable.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

#include "G4GlobalMagFieldMessenger.hh"
#include "G4AutoDelete.hh"

#include "G4UserLimits.hh"
//#include "LXSetUp.hh"

#include "LxBFields.hh"


FieldData2D::FieldData2D(const std::string fname) : FieldDistribution(), fFName(fname), fNsort(10)
{
  LoadFieldData(fFName);
  frMin = fFieldData[0];
  frMax = fFieldData[0];
  for (const auto &rr : fFieldData) {
    if (frMin.x() > rr.x()) frMin.setX(rr.x());
    if (frMin.y() > rr.y()) frMin.setY(rr.y());
    if (frMax.x() < rr.x()) frMax.setX(rr.x());
    if (frMax.y() < rr.y()) frMax.setY(rr.y());
  }
}



void FieldData2D::LoadFieldData(const std::string fname)
{
  std::fstream  fdataf;

  fdataf.open(fname, std::ios::in);
  if (!fdataf.is_open()) {
//     G4String msgstr("Failed to load magnetic field data from the file ");
//     msgstr += fname;
//     G4Exception("FieldData2D::", "LoadFieldData()", FatalException, msgstr.c_str());
  }

  double x, z, bf;
  while (!fdataf.eof()) {
      fdataf >> x >> z >> bf;
      if (fdataf.eof()) break;
      fFieldData.push_back(G4ThreeVector(x*mm, z*mm, bf*tesla));
//       std::cout << x << " " << z << "  " << bf << std::endl;
  }
  fdataf.close();
  G4cout << fFieldData.size() << " Field data points were loaded from the file " << fname << G4endl;
}



G4double FieldData2D::DataField(const G4double xx, const G4double yy) const
{
  if (!(xx > frMin.x() && xx < frMax.x() && yy > frMin.y() && yy < frMax.y()) ) return 0.0;
  
  std::vector<G4double> drv;
  std::for_each(fFieldData.begin(), fFieldData.end(), 
                [&](const G4ThreeVector &vv){drv.push_back( pow(vv.x()-xx, 2.0) + pow(vv.y()-yy, 2.0 ));} );
  std::vector<size_t> vind(drv.size());
  size_t ind = 0;
  std::generate(vind.begin(), vind.end(), [&](){return ind++;});
  point_cmp<G4double> pcmp(&drv);
  std::partial_sort (vind.begin(), vind.begin()+fNsort, vind.end(), pcmp);
   
  size_t ind3 = 2;
  while (CheckCollinear(fFieldData[vind[0]], fFieldData[vind[1]], fFieldData[vind[ind3]]) 
//          && (!CheckInside(fFieldData[vind[0]], fFieldData[vind[1]], fFieldData[vind[ind3]], xx, yy)) 
        ) {
    ind3++;
    if (ind3 > fNsort) { 
      G4String msgstr("Number of elements for partial sorting is not sufficient");
      G4Exception("PDSMagnetField::", "DataField()", FatalException, msgstr.c_str());
    }
  }

//   std::cout << "DataField: " << xx << "  " << yy << "  " << drv[vind[0]] 
//             << "  " << drv[vind[1]] << "  " << drv[vind[ind3]] << std::endl;
  return FindZ(fFieldData.at(vind[0]), fFieldData.at(vind[1]), fFieldData.at(vind[ind3]), xx, yy);
}



G4double FieldData2D::FindZ(const G4ThreeVector &r1, const G4ThreeVector &r2, const G4ThreeVector &r3, 
                               const double xx, const double yy) const
{
// Makes linear interpolation finding Z coordinate for a given (x,y) 
// as a point on the surface defined by three 3D points r1, r2, r3.
  G4ThreeVector dr1 = r2 -r1;
  G4ThreeVector dr2 = r3 -r1;

  G4ThreeVector normv = dr1.cross(dr2).unit();

  return r1.z() - (normv.x()*(xx-r1.x()) + normv.y()*(yy-r1.y())) / normv.z();
}



G4bool FieldData2D::CheckCollinear(const G4ThreeVector &r1, const G4ThreeVector &r2, const G4ThreeVector &r3) const
{
// Check if 2D projection to xy-plane of three 3D points are on one line.
  const double eps = 1.0e-15;
  G4ThreeVector dr1 = r2 -r1;
  G4ThreeVector dr2 = r3 -r1;
  dr1.setZ(0.0);  
  dr2.setZ(0.0);
  return ((dr1.cross(dr2).mag2()) < eps);
}


/////////////////////////////////////////////////////////////////////
LxBField::LxBField(const VFieldComponent *bx, const VFieldComponent *by, 
                   const VFieldComponent *bz, const G4ThreeVector pos):
  G4MagneticField(), fPos(pos)
{
  fFieldXYZ.push_back(bx);
  fFieldXYZ.push_back(by);
  fFieldXYZ.push_back(bz);
}




