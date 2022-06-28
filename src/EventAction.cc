/***************************************************************
 * descriptions of EventAction
 * Author  : Hualin Xiao
 * Date    : Jan., 2015
 * Version : 1.10
 ***************************************************************/

#include "G4Event.hh"
#include "G4GenericMessenger.hh"
#include "G4RunManager.hh"

#include <iomanip>

#include "AnalysisManager.hh"
#include "EventAction.hh"

EventAction::EventAction() {}

EventAction::~EventAction() {}

void EventAction::BeginOfEventAction(const G4Event *event) {
  AnalysisManager *analysisManager = AnalysisManager::GetInstance();

  analysisManager->InitEvent(event);
  fEventID = event->GetEventID();
  if (fEventID % 100000 == 0) {
    G4cout << "\n---> Begin of event: " << fEventID << G4endl;
  }
}

void EventAction::EndOfEventAction(const G4Event *event) {
  AnalysisManager *analysisManager = AnalysisManager::GetInstance();
  analysisManager->ProcessEvent(event);
  if (fEventID % 100000 == 0) {
    G4cout << "---> End of event: " << fEventID << G4endl;
  }
}
