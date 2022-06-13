//
/// \file SteppingAction.hh
/// \brief Definition of the SteppingAction class

#ifndef SteppingAction_h
#define SteppingAction_h 1

#include "G4UserSteppingAction.hh"

// class Birks;
class G4EmSaturation;

class SteppingAction : public G4UserSteppingAction {
public:
  SteppingAction();
  virtual ~SteppingAction();

  virtual void UserSteppingAction(const G4Step *aStep);

private:
  //	Birks *birks;
  G4EmSaturation *emSat;
};

#endif
