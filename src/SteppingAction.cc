/***************************************************************
 * \brief Implementation of the SteppingAction class
 * Author  : Hualin Xiao
 * Date    : Feb., 2015
 * Version : 1.10
 ***************************************************************/
#include "G4SystemOfUnits.hh"
#include "TString.h"

#include "EventAction.hh"
#include "G4Step.hh"
#include "G4Track.hh"
//#include "Birks.hh"
#include "AnalysisManager.hh"
#include "SteppingAction.hh"

SteppingAction::SteppingAction() : G4UserSteppingAction() {}

SteppingAction::~SteppingAction() {}

void SteppingAction::UserSteppingAction(const G4Step *aStep) {

  AnalysisManager *analysisManager = AnalysisManager::GetInstance();

  analysisManager->processStep(aStep);
}
