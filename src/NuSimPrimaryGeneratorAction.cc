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
/// \file eventgenerator/HepMC/HepMCEx01/src/NuSimPrimaryGeneratorAction.cc
/// \brief Implementation of the NuSimPrimaryGeneratorAction class
//
// $Id: NuSimPrimaryGeneratorAction.cc 77801 2013-11-28 13:33:20Z gcosmo $
//

#include "NuSimPrimaryGeneratorAction.hh"
#include "NuSimPrimaryGeneratorMessenger.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "HepMCG4AsciiReader.hh"
#include "HepMCG4PythiaInterface.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
NuSimPrimaryGeneratorAction::NuSimPrimaryGeneratorAction()
{
  // default generator is particle gun.
  currentGenerator = particleGun= new G4ParticleGun();
  currentGeneratorName = "particleGun";
  hepmcAscii = new HepMCG4AsciiReader();
#ifdef G4LIB_USE_PYTHIA
  pythiaGen = new HepMCG4PythiaInterface();
#else
  pythiaGen = 0;
#endif
  gentypeMap["particleGun"] = particleGun;
  gentypeMap["hepmcAscii"] = hepmcAscii;
  gentypeMap["pythia"] = pythiaGen;

  messenger= new NuSimPrimaryGeneratorMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
NuSimPrimaryGeneratorAction::~NuSimPrimaryGeneratorAction()
{
  delete messenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void NuSimPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  if(currentGenerator)
    currentGenerator-> GeneratePrimaryVertex(anEvent);
  else
    G4Exception("NuSimPrimaryGeneratorAction::GeneratePrimaries",
                "PrimaryGeneratorAction001", FatalException,
                "generator is not instanciated." );
}
