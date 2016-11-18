// -*- C++ -*-
//
// Package:    PulseStudies/PrintRecHitSampleFourCorrections
// Class:      PrintRecHitSampleFourCorrections
// 
/**\class PrintRecHitSampleFourCorrections PrintRecHitSampleFourCorrections.cc PulseStudies/PrintRecHitSampleFourCorrections/plugins/PrintRecHitSampleFourCorrections.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/


// system include files
#include <memory>
#include <iostream>
#include <algorithm>
#include <vector>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Provenance/interface/RunLumiEventNumber.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"

#include "CalibCalorimetry/EcalLaserCorrection/interface/EcalLaserDbService.h"
#include "CalibCalorimetry/EcalLaserCorrection/interface/EcalLaserDbRecord.h"
#include "CondFormats/EcalObjects/interface/EcalADCToGeVConstant.h"
#include "CondFormats/DataRecord/interface/EcalADCToGeVConstantRcd.h"
#include "CondFormats/EcalObjects/interface/EcalIntercalibConstants.h"
#include "CondFormats/EcalObjects/interface/EcalIntercalibConstantsMC.h"
#include "CondFormats/DataRecord/interface/EcalIntercalibConstantsRcd.h"
#include "CondFormats/DataRecord/interface/EcalIntercalibConstantsMCRcd.h"

#include "PulseStudies/PulseTree/interface/RecHitSampleFourCorrector.h"

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class PrintRecHitSampleFourCorrections : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit PrintRecHitSampleFourCorrections(const edm::ParameterSet&);
      ~PrintRecHitSampleFourCorrections();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

  RecHitSampleFourCorrector corrector_;

      // ----------member data ---------------------------

  edm::EDGetTokenT<EcalRecHitCollection> tok_ebrechits_;
  edm::EDGetTokenT<EcalRecHitCollection> tok_eerechits_;

  edm::ESHandle<EcalLaserDbService> laser;
  edm::ESHandle<EcalIntercalibConstants> ical;
  edm::ESHandle<EcalADCToGeVConstant> agc;

  edm::RunNumber_t _run;
  edm::LuminosityBlockNumber_t _lumi;
  edm::EventNumber_t _evt;

  void PrintRecHit(const EcalRecHit *rh, edm::Timestamp runTime_);

};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
PrintRecHitSampleFourCorrections::PrintRecHitSampleFourCorrections(const edm::ParameterSet& iConfig)

{

   //now do what ever initialization is needed
   usesResource("TFileService");
   edm::Service<TFileService> fs;

   tok_ebrechits_ = consumes<EcalRecHitCollection>(edm::InputTag("reducedEgamma","reducedEBRecHits"));
   tok_eerechits_ = consumes<EcalRecHitCollection>(edm::InputTag("reducedEgamma","reducedEERecHits"));

}


PrintRecHitSampleFourCorrections::~PrintRecHitSampleFourCorrections()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
PrintRecHitSampleFourCorrections::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   iSetup.get<EcalLaserDbRecord>().get(laser);
   iSetup.get<EcalIntercalibConstantsRcd>().get(ical);
   iSetup.get<EcalADCToGeVConstantRcd>().get(agc);

   Handle<EcalRecHitCollection> ebrechithandle;
   iEvent.getByToken(tok_ebrechits_,ebrechithandle);
   auto ebrechits = ebrechithandle.product();
   Handle<EcalRecHitCollection> eerechithandle;
   iEvent.getByToken(tok_eerechits_,eerechithandle);
   auto eerechits = eerechithandle.product();

   auto thistime = iEvent.eventAuxiliary().time();

   _run = iEvent.eventAuxiliary().run();
   _lumi = iEvent.eventAuxiliary().luminosityBlock();
   _evt = iEvent.eventAuxiliary().event();

   for (auto rh = ebrechits->begin(); rh!=ebrechits->end(); rh++) {
     const EcalRecHit* _rh = &(*rh);
     PrintRecHit(_rh,thistime);
   }
   for (auto rh = eerechits->begin(); rh!=eerechits->end(); rh++) {
     const EcalRecHit* _rh = &(*rh);
     PrintRecHit(_rh,thistime);
   }

}

void PrintRecHitSampleFourCorrections::PrintRecHit(const EcalRecHit *rh, edm::Timestamp runTime_){
  auto lc = laser->getLaserCorrection(rh->id(), runTime_);
  auto icalMap = ical->getMap();
  auto it = icalMap.find(rh->id());
  auto ic = (it!=icalMap.end()) ? (*it) : 0;
  auto _agc = (rh->id().subdetId()==EcalBarrel) ? float(agc->getEBValue()) : float(agc->getEEValue());
  std::cout << _run << " " << _lumi << " " << _evt << " " << rh->id().rawId() << " " << rh->energy() << " " << corrector_.RecHitCorrectedEnergy(rh,lc,ic,_agc) << std::endl;
}

// ------------ method called once each job just before starting event loop  ------------
void 
PrintRecHitSampleFourCorrections::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
PrintRecHitSampleFourCorrections::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
PrintRecHitSampleFourCorrections::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(PrintRecHitSampleFourCorrections);
