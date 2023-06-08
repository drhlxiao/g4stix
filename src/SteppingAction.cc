/***************************************************************
 * \brief Implementation of the SteppingAction class
 * Author  : Hualin Xiao
 * Date    : Feb., 2015
 * Version : 1.10
 ***************************************************************/
#include "SteppingAction.hh"

#include "AnalysisManager.hh"
#include "EventAction.hh"
#include "G4Step.hh"
#include "G4SystemOfUnits.hh"
#include "G4Track.hh"
#include "TString.h"

SteppingAction::SteppingAction() : G4UserSteppingAction() {}

SteppingAction::~SteppingAction() {}

void SteppingAction::UserSteppingAction(const G4Step *aStep) {
  AnalysisManager *analysisManager = AnalysisManager::GetInstance();

  analysisManager->ProcessStep(aStep);
}
