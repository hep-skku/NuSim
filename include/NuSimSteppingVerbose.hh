#ifndef NuSimSteppingVerbose_h
#define NuSimSteppingVerbose_h 1

#include "G4SteppingVerbose.hh"

class NuSimSteppingVerbose;


class NuSimSteppingVerbose : public G4SteppingVerbose {
public:
  NuSimSteppingVerbose();
  ~NuSimSteppingVerbose();

  void StepInfo();
  void TrackingStarted();

};

#endif
