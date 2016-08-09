//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
/// \file eventgenerator/HepMC/HepMCEx01/src/NuSimPrimaryGeneratorMessenger.cc
/// \brief Implementation of the NuSimPrimaryGeneratorMessenger class
//
// $Id: NuSimPrimaryGeneratorMessenger.cc 77801 2013-11-28 13:33:20Z gcosmo $
//

#include "G4UIcommand.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIdirectory.hh"
#include "G4UIparameter.hh"
#include "NuSimPrimaryGeneratorMessenger.hh"
#include "NuSimPrimaryGeneratorAction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
NuSimPrimaryGeneratorMessenger::NuSimPrimaryGeneratorMessenger
                                (NuSimPrimaryGeneratorAction* genaction)
  : primaryAction(genaction)
{
  mydetdir = new G4UIdirectory("/mydet/");
  mydetdir-> SetGuidance("NuSim detector control commands.");

  dir= new G4UIdirectory("/generator/");
  dir-> SetGuidance("Control commands for primary generator");

  select= new G4UIcmdWithAString("/generator/select", this);
  select-> SetGuidance("select generator type");
  select-> SetParameterName("generator_type", false, false);
  select-> SetCandidates("particleGun pythia hepmcAscii");
  select-> SetDefaultValue("particleGun");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
NuSimPrimaryGeneratorMessenger::~NuSimPrimaryGeneratorMessenger()
{
  delete select;
  delete mydetdir;
  delete dir;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void NuSimPrimaryGeneratorMessenger::SetNewValue(G4UIcommand* command,
                                              G4String newValues)
{
  if ( command == select ) {
    primaryAction-> SetGenerator(newValues);
    G4cout << "current generator type: "
            << primaryAction-> GetGeneratorName() << G4endl;
  }
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4String NuSimPrimaryGeneratorMessenger::GetCurrentValue(G4UIcommand* command)
{
  G4String cv, st;
  if ( command == select ) {
    cv = primaryAction-> GetGeneratorName();
  }

 return cv;
}
