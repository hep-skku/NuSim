#include "FTFP_BERT.hh"
#include "G4RunManager.hh"
#include "G4UImanager.hh"

#include "NuSimDetectorConstruction.hh"
#include "NuSimPrimaryGeneratorAction.hh"
/*
#include "NuSimEventAction.hh"
#include "NuSimRunAction.hh"
#include "NuSimStackingAction.hh"
#include "NuSimSteppingVerbose.hh"
#include "NuSimSteppingAction.hh"
#include "NuSimTrackingAction.hh"
*/

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

#ifdef G4UI_USE
#include "G4UIExecutive.hh"
#endif

int main(int argc,char** argv)
{
  //G4VSteppingVerbose* verbosity = new NuSimSteppingVerbose;
  //G4VSteppingVerbose::SetInstance(verbosity);

  G4RunManager* runManager = new G4RunManager;

  G4VUserDetectorConstruction* detector = new NuSimDetectorConstruction;
  runManager->SetUserInitialization(detector);

  G4VUserPhysicsList* physics = new FTFP_BERT;
  runManager->SetUserInitialization(physics);

  runManager->Initialize();

  G4VUserPrimaryGeneratorAction* gen_action = new NuSimPrimaryGeneratorAction;
  runManager->SetUserAction(gen_action);
/*
  //
  G4UserRunAction* run_action = new NuSimRunAction;
  runManager->SetUserAction(run_action);
  //
  G4UserEventAction* event_action = new NuSimEventAction;
  runManager->SetUserAction(event_action);
  //
  G4UserStackingAction* stacking_action = new NuSimStackingAction;
  runManager->SetUserAction(stacking_action);
  //
  G4UserTrackingAction* tracking_action = new NuSimTrackingAction;
  runManager->SetUserAction(tracking_action);
  //
  G4UserSteppingAction* stepping_action = new NuSimSteppingAction;
  runManager->SetUserAction(stepping_action);
*/

#ifdef G4VIS_USE
  // Visualization, if you choose to have it!
  G4VisManager* visManager = new G4VisExecutive;
  visManager->Initialize();
#endif

  G4UImanager* UImanager = G4UImanager::GetUIpointer();
  if (argc!=1) {  // batch mode
#ifdef G4VIS_USE
    visManager->SetVerboseLevel("quiet");
#endif
    G4String command = "/control/execute ";
    G4String fileName = argv[1];
    UImanager->ApplyCommand(command+fileName);
  }
  else { // interactive mode : define UI session
#ifdef G4UI_USE
    G4UIExecutive* ui = new G4UIExecutive(argc, argv);
    ui->SessionStart();
    delete ui;
#endif
  }

#ifdef G4VIS_USE
  delete visManager;
#endif
  delete runManager;
  //delete verbosity;

  return 0;
}

