#include "G4Box.hh"
#include "G4Colour.hh"
#include "G4Element.hh"
#include "G4ElementTable.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4PVParameterised.hh"
#include "G4PVPlacement.hh"
#include "G4ThreeVector.hh"
#include "G4Tubs.hh"
#include "G4RotationMatrix.hh"
#include "G4Transform3D.hh"
#include "G4TransportationManager.hh"
#include "G4SDManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4VisAttributes.hh"
//#include "NuSimCalorimeterParametrisation.hh"
//#include "NuSimCalorimeterROGeometry.hh"
//#include "NuSimCalorimeterSD.hh"
#include "NuSimDetectorConstruction.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
NuSimDetectorConstruction::NuSimDetectorConstruction()
{
#include "NuSimDetectorParameterDef.icc"
  DefineMaterials();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
NuSimDetectorConstruction::~NuSimDetectorConstruction()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void NuSimDetectorConstruction::DefineMaterials()
{
  //-------------------------------------------------------------------------
  // Materials
  //-------------------------------------------------------------------------

  G4double a, z, density;
  G4int nel;

  G4Element* H = new G4Element("Hydrogen", "H", z=1., a=  1.01*g/mole);
  G4Element* C = new G4Element("Carbon",   "C", z=6., a= 12.01*g/mole);
  G4Element* N = new G4Element("Nitrogen", "N", z=7., a= 14.01*g/mole);
  G4Element* O = new G4Element("Oxygen",   "O", z=8., a= 16.00*g/mole);

  fAir = new G4Material("Air", density= 1.29*mg/cm3, nel=2);
  fAir-> AddElement(N, 70.*perCent);
  fAir-> AddElement(O, 30.*perCent);

  fLead = new G4Material("Lead", z=82., a= 207.19*g/mole, density= 11.35*g/cm3);

  fAr = new G4Material("ArgonGas",z=18., a= 39.95*g/mole, density=1.782*mg/cm3);

  fSilicon = new G4Material("Silicon", z=14., a= 28.09*g/mole,
                                       density= 2.33*g/cm3);

  fScinti = new G4Material("Scintillator", density= 1.032*g/cm3, nel=2);
  fScinti-> AddElement(C, 9);
  fScinti-> AddElement(H, 10);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4VPhysicalVolume* NuSimDetectorConstruction::Construct()
{
  //------------------------------ experimental hall
  G4Box* experimentalHall_box =
         new G4Box("expHall_b", fexpHall_x, fexpHall_y, fexpHall_z);
  G4LogicalVolume* experimentalHall_log =
         new G4LogicalVolume(experimentalHall_box, fAir,"expHall_L", 0,0,0);
  G4VPhysicalVolume * experimentalHall_phys =
         new G4PVPlacement(0, G4ThreeVector(), experimentalHall_log,
                           "expHall_P", 0, false,0);
  G4VisAttributes* experimentalHallVisAtt =
         new G4VisAttributes(G4Colour(1.,1.,1.));
  experimentalHallVisAtt-> SetForceWireframe(true);
  experimentalHall_log-> SetVisAttributes(experimentalHallVisAtt);

/*
  //------------------------------ calorimeter
  G4VSolid* calorimeter_tubs
    = new G4Tubs("calorimeter_tubs", fcaloTubs_rmin, fcaloTubs_rmax,
                  fcaloTubs_dz, fcaloTubs_sphi, fcaloTubs_dphi);
  G4LogicalVolume* calorimeter_log
    = new G4LogicalVolume(calorimeter_tubs, fScinti, "caloT_L",0,0,0);
  // G4VPhysicalVolume * calorimeter_phys =
      new G4PVPlacement(0,G4ThreeVector(), calorimeter_log, "caloM_P",
                        experimentalHall_log, false,0);
  G4VisAttributes* calorimeter_logVisATT
    = new G4VisAttributes(G4Colour(1.0,1.0,0.0));
  calorimeter_logVisATT->SetForceWireframe(true);
  calorimeter_log->SetVisAttributes(calorimeter_logVisATT);

  //------------------------------- Lead layers
  // As an example for Parameterised volume
  // dummy values for G4Tubs -- modified by parameterised volume
  G4VSolid* caloLayer_tubs
    = new G4Tubs("caloLayer_tubs", fcaloRing_rmin, fcaloRing_rmax,
                  fcaloRing_dz, fcaloRing_sphi, fcaloRing_dphi);
  G4LogicalVolume* caloLayer_log
    = new G4LogicalVolume(caloLayer_tubs, fLead, "caloR_L",0,0,0);
  G4VPVParameterisation* calorimeterParam
    = new NuSimCalorimeterParametrisation;
  // dummy value : kXAxis -- modified by parameterised volume
  // G4VPhysicalVolume * caloLayer_phys =
      new G4PVParameterised("caloLayer_phys",caloLayer_log,calorimeter_log,
                           kXAxis, fnocaloLayers, calorimeterParam);
  G4VisAttributes* caloLayer_logVisAtt
    = new G4VisAttributes(G4Colour(0.7,1.0,0.0));
  caloLayer_logVisAtt->SetForceWireframe(true);
  caloLayer_log->SetVisAttributes(caloLayer_logVisAtt);

  //------------------------------------------------------------------
  // Sensitive Detector
  //------------------------------------------------------------------

  G4SDManager* SDman = G4SDManager::GetSDMpointer();

  G4String calorimeterSDname = "/mydet/calorimeter";
  NuSimCalorimeterSD* calorimeterSD = new NuSimCalorimeterSD(calorimeterSDname);
  G4String ROgeometryName = "CalorimeterROGeom";
  G4VReadOutGeometry* calRO = new NuSimCalorimeterROGeometry(ROgeometryName);
  calRO->BuildROGeometry();
  calorimeterSD->SetROgeometry(calRO);
  SDman->AddNewDetector(calorimeterSD);
  calorimeter_log->SetSensitiveDetector(calorimeterSD);
*/

  //------------------------------------------------------------------
  // Digitizer modules
  //------------------------------------------------------------------

  return experimentalHall_phys;
}
