// -*- C++ -*-
//
// Package:    PulseStudies/PulseTree
// Class:      PulseTree
// 
/**\class PulseTree PulseTree.cc PulseStudies/PulseTree/plugins/PulseTree.cc

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
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/EcalDigi/interface/EcalDigiCollections.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "CondFormats/EcalObjects/interface/EcalMGPAGainRatio.h"
#include "CondFormats/EcalObjects/interface/EcalGainRatios.h"
#include "CondFormats/DataRecord/interface/EcalGainRatiosRcd.h"

#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/EcalAlgo/interface/EcalBarrelGeometry.h"
#include "Geometry/EcalAlgo/interface/EcalEndcapGeometry.h"

#include "TTree.h"

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class PulseTree : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit PulseTree(const edm::ParameterSet&);
      ~PulseTree();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

  void FillDigi(EcalDataFrame digi, const EcalUncalibratedRecHitCollection *rechits, const EcalUncalibratedRecHitCollection *w_rechits);
  bool FilterBx(unsigned bx);
  void WriteAverageOutput();

  float MultiFitParametricCorrection(float amplitude_multifit_intime_uncal, float eta);

      // ----------member data ---------------------------

  edm::EDGetTokenT<EBDigiCollection> tok_ebdigis_;
  edm::EDGetTokenT<EEDigiCollection> tok_eedigis_;
  edm::EDGetTokenT<EcalUncalibratedRecHitCollection> tok_ebrechits_;
  edm::EDGetTokenT<EcalUncalibratedRecHitCollection> tok_w_ebrechits_;
  edm::EDGetTokenT<EcalUncalibratedRecHitCollection> tok_eerechits_;
  edm::EDGetTokenT<EcalUncalibratedRecHitCollection> tok_w_eerechits_;
  bool do_EB;
  bool do_EE;
  bool do_average;
  bool split_lumis;
  bool do_filter_bx;
  UInt_t n_pedestal_samples;
  bool save_rechits;

  TTree *outTree;
  UInt_t t_run;
  UShort_t t_lumi;
  UShort_t t_bx;
  UInt_t t_id;
  float t_gainratios[2];
  float t_eta;
  float t_phi;
  float t_pulse[10];
  float t_gain[10];
  float t_pedestal;
  float t_pedestal_rms;
  UShort_t t_gainmask;
  UInt_t t_nevt;
  float t_multifit[10];
  float t_weights;
  float t_corrected_multifit;

  std::vector<std::pair<UInt_t,UShort_t> > summed_index; // key: (id,gain[0])
  std::vector<std::vector<double> > summed_pulses;
  std::vector<std::vector<double> > summed_gains;
  std::vector<std::pair<double,double> > summed_pedestal;
  std::vector<long> summed_count;
  float min_amplitude_in_average;

  float min_rechit_amplitude_weights;

  UShort_t old_lumi;
  std::vector<unsigned> bx_to_keep;
  bool bx_invert_selection;
  std::vector<UShort_t> seen_lumis;

  edm::ESHandle<CaloGeometry> geometry;
  edm::ESHandle<EcalGainRatios> gRatio;

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
PulseTree::PulseTree(const edm::ParameterSet& iConfig)

{

   //now do what ever initialization is needed
   usesResource("TFileService");
   edm::Service<TFileService> fs;

   n_pedestal_samples = iConfig.getUntrackedParameter<unsigned int>("nPedestalSamples",3);
   do_average = iConfig.getUntrackedParameter<bool>("doAverage",false);
   split_lumis = iConfig.getUntrackedParameter<bool>("splitByLumi",false) && do_average;
   min_amplitude_in_average = iConfig.getUntrackedParameter<double>("minAmplitudeForAverage",-9e9);
   old_lumi = 0;

   min_rechit_amplitude_weights = iConfig.getUntrackedParameter<double>("minRecHitAmplitudeWeights",-9e9);

   do_EB = iConfig.getUntrackedParameter<bool>("processEB",true);
   do_EE = iConfig.getUntrackedParameter<bool>("processEE",true);
   if (do_EB) tok_ebdigis_ = consumes<EBDigiCollection>(iConfig.getParameter<edm::InputTag>("EBDigiCollection"));
   if (do_EE) tok_eedigis_ = consumes<EEDigiCollection>(iConfig.getParameter<edm::InputTag>("EEDigiCollection"));

   save_rechits = iConfig.getUntrackedParameter<bool>("saveRecHits",false) && (!do_average);
   if (save_rechits) {
     tok_ebrechits_ = consumes<EcalUncalibratedRecHitCollection>(edm::InputTag("ecalMultiFitUncalibRecHit","EcalUncalibRecHitsEB"));
     tok_w_ebrechits_ = consumes<EcalUncalibratedRecHitCollection>(edm::InputTag("ecalUncalibRecHit","EcalUncalibRecHitsEB"));
     tok_eerechits_ = consumes<EcalUncalibratedRecHitCollection>(edm::InputTag("ecalMultiFitUncalibRecHit","EcalUncalibRecHitsEE"));
     tok_w_eerechits_ = consumes<EcalUncalibratedRecHitCollection>(edm::InputTag("ecalUncalibRecHit","EcalUncalibRecHitsEE"));
   }

   outTree = fs->make<TTree>("pulses","pulses");
   outTree->Branch("run",&t_run,"run/i");
   if (!do_average || split_lumis) outTree->Branch("lumi",&t_lumi,"lumi/s");
   if (!do_average) outTree->Branch("bx",&t_bx,"bx/s");
   outTree->Branch("id",&t_id,"id/i");
   outTree->Branch("gainratios",t_gainratios,"gainratios[2]/F");
   outTree->Branch("eta",&t_eta,"eta/F");
   outTree->Branch("phi",&t_phi,"phi/F");
   outTree->Branch("pulse",t_pulse,"pulse[10]/F");
   outTree->Branch("gain",t_gain,"gain[10]/F");
   outTree->Branch("pedestal",&t_pedestal,"pedestal/F");
   if (do_average) outTree->Branch("pedestal_rms",&t_pedestal_rms,"pedestal_rms/F");
   if (!do_average) outTree->Branch("gainmask",&t_gainmask,"gainmask/s");
   if (do_average) outTree->Branch("nevt",&t_nevt,"nevt/i");
   if (save_rechits) {
     outTree->Branch("ampl_multifit",t_multifit,"ampl_multifit[10]/F");
     outTree->Branch("ampl_weights",&t_weights,"ampl_weights/F");
     outTree->Branch("ampl_corrected_multifit",&t_corrected_multifit,"ampl_corrected_multifit/F");
   }

   bx_to_keep = iConfig.getUntrackedParameter<std::vector<unsigned> >("filterBx",std::vector<unsigned>());
   bx_invert_selection = iConfig.getUntrackedParameter<bool>("invertBxSelection",false);
   do_filter_bx = (bx_to_keep.size()!=0);

}


PulseTree::~PulseTree()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
PulseTree::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   iSetup.get<CaloGeometryRecord>().get(geometry);
   iSetup.get<EcalGainRatiosRcd>().get(gRatio);

   t_run = iEvent.eventAuxiliary().run();
   t_lumi = iEvent.eventAuxiliary().luminosityBlock();
   t_bx = iEvent.eventAuxiliary().bunchCrossing();
   if (do_filter_bx && !FilterBx(t_bx)) return;

   if (split_lumis && (t_lumi!=old_lumi) && old_lumi!=0) {
     if (std::find(seen_lumis.begin(),seen_lumis.end(),t_lumi)!=seen_lumis.end()){
       std::cout << "Warning: lumisection " << t_lumi << " already found previously, will be filled two times" << std::endl;
     }
     WriteAverageOutput();
     summed_index.clear();
     summed_pulses.clear();
     summed_gains.clear();
     summed_pedestal.clear();
     summed_count.clear();
     old_lumi=t_lumi;
   }

   Handle<EBDigiCollection> ebdigihandle;
   Handle<EEDigiCollection> eedigihandle;
   Handle<EcalUncalibratedRecHitCollection> ebrechithandle;
   Handle<EcalUncalibratedRecHitCollection> w_ebrechithandle;
   const EcalUncalibratedRecHitCollection *ebrechits = NULL;
   const EcalUncalibratedRecHitCollection *w_ebrechits = NULL;
   Handle<EcalUncalibratedRecHitCollection> eerechithandle;
   Handle<EcalUncalibratedRecHitCollection> w_eerechithandle;
   const EcalUncalibratedRecHitCollection *eerechits = NULL;
   const EcalUncalibratedRecHitCollection *w_eerechits = NULL;
   if (save_rechits){
     iEvent.getByToken(tok_ebrechits_,ebrechithandle);
     ebrechits = ebrechithandle.product();
     iEvent.getByToken(tok_w_ebrechits_,w_ebrechithandle);
     w_ebrechits = w_ebrechithandle.product();
     iEvent.getByToken(tok_eerechits_,eerechithandle);
     eerechits = eerechithandle.product();
     iEvent.getByToken(tok_w_eerechits_,w_eerechithandle);
     w_eerechits = w_eerechithandle.product();
   }
   if (do_EB){
     iEvent.getByToken(tok_ebdigis_,ebdigihandle);
     auto ebdigis = ebdigihandle.product();
     for (uint i=0; i<ebdigis->size(); i++) FillDigi((*ebdigis)[i],ebrechits,w_ebrechits);
   }
   if (do_EE){
     iEvent.getByToken(tok_eedigis_,eedigihandle);
     auto eedigis = eedigihandle.product();
     for (uint i=0; i<eedigis->size(); i++) FillDigi((*eedigis)[i],eerechits,w_eerechits);
   }

}

void
PulseTree::FillDigi(EcalDataFrame digi, const EcalUncalibratedRecHitCollection *rechits, const EcalUncalibratedRecHitCollection *w_rechits){

  t_id = UInt_t(digi.id());
  t_gainmask = 0;

  for (int j=0; j<10; j++){ 
    t_pulse[j] = float(digi[j]&0xFFF);
    t_gain[j] = float((digi[j]>>12)&0x3);
    t_gainmask |= 1<<(int(t_gain[j]));
  }

  t_pedestal = 0;
  for (uint j=0; j<n_pedestal_samples; j++) t_pedestal += t_pulse[j];
  t_pedestal/=n_pedestal_samples;
  
  t_pedestal_rms = 0;
  t_nevt = 0;

  auto _gr = gRatio.product();
  t_gainratios[0]=(*_gr)[t_id].gain12Over6();
  t_gainratios[1]=t_gainratios[0]*(*_gr)[t_id].gain6Over1();

  if (!do_average) {
    auto detid = DetId(t_id);
    if (save_rechits){
      for (int j=0; j<10; j++) t_multifit[j] = 0;
      t_weights = 0;
      t_corrected_multifit = 0;
      auto subGeom =  geometry->getSubdetectorGeometry(detid);
      auto cellGeom = subGeom->getGeometry(detid);
      t_eta = cellGeom->getPosition().eta();
      t_phi = cellGeom->getPosition().phi();
      auto it = rechits->find(detid);
      if (it==rechits->end()) std::cout << "Warning: rechit (multifit) not found" << std::endl;
      else for (int j=0; j<10; j++) t_multifit[j] = (j==5) ? it->amplitude() : it->outOfTimeAmplitude(j);
      auto it2 = w_rechits->find(detid);
      if (it2==w_rechits->end()) std::cout << "Warning: rechit (weights) not found" << std::endl;
      else t_weights = it2->amplitude();
      t_corrected_multifit = it->amplitude()*MultiFitParametricCorrection(it->amplitude(),t_eta);
      if (t_weights<min_rechit_amplitude_weights) return;
    }
    outTree->Fill();
  }
  else {
    if ((-t_pedestal + *std::max_element(t_pulse,t_pulse+10)) < min_amplitude_in_average) return; // do not average empty pulses
    auto thispair = std::pair<UInt_t,UShort_t>(t_id,int(t_gain[0]));
    auto it = find(summed_index.begin(),summed_index.end(),thispair);
    if (it==summed_index.end()){
      summed_index.push_back(thispair);
      summed_pulses.push_back(std::vector<double>(10,0));
      summed_gains.push_back(std::vector<double>(10,0));
      summed_pedestal.push_back(std::pair<double,double>(0,0));
      summed_count.push_back(0);
      it = find(summed_index.begin(),summed_index.end(),thispair);
    }
    uint idx = it-summed_index.begin();
    for (int j=0; j<10; j++){
      summed_pulses[idx][j]+=t_pulse[j];
      summed_gains[idx][j]+=t_gain[j];
    }
    summed_pedestal[idx] = std::pair<double,double>(summed_pedestal[idx].first+t_pedestal,summed_pedestal[idx].second+std::pow(t_pedestal,2));
    summed_count[idx]+=1;
  }

}

bool
PulseTree::FilterBx(unsigned bx){
  for (auto x : bx_to_keep) if (bx==x) return (!bx_invert_selection);
  return bx_invert_selection;
}

void
PulseTree::WriteAverageOutput(){

  t_gainmask = 0;
  t_bx = 0;
  for (int j=0; j<10; j++) t_multifit[j] = 0;
  t_weights = 0;
  t_corrected_multifit = 0;

  for (uint idx = 0; idx<summed_count.size(); idx++){
    float den = summed_count[idx];
    for (int j=0; j<10; j++){
      t_pulse[j] = summed_pulses[idx][j]/den;
      t_gain[j] = summed_gains[idx][j]/den;
    }
    t_pedestal = summed_pedestal[idx].first/den;
    t_pedestal_rms = std::sqrt(summed_pedestal[idx].second/den - std::pow(summed_pedestal[idx].first/den,2));
    t_id = summed_index[idx].first;
    auto _gr = gRatio.product();
    t_gainratios[0]=(*_gr)[t_id].gain12Over6();
    t_gainratios[1]=t_gainratios[0]*(*_gr)[t_id].gain6Over1();
    auto detid = DetId(t_id);
    auto subGeom =  geometry->getSubdetectorGeometry(detid);
    auto cellGeom = subGeom->getGeometry(detid);
    t_eta = cellGeom->getPosition().eta();
    t_phi = cellGeom->getPosition().phi();
    t_nevt = summed_count[idx];
    outTree->Fill();
  }

}

float
PulseTree::MultiFitParametricCorrection(float amplitude_multifit_intime_uncal, float eta){

  if (fabs(eta)<1.479){
    if (amplitude_multifit_intime_uncal-2800.000000<200) return 1.000000;
    else if (amplitude_multifit_intime_uncal-3200.000000<200) return 0.999503;
    else if (amplitude_multifit_intime_uncal-3600.000000<200) return 0.998589;
    else if (amplitude_multifit_intime_uncal-4000.000000<200) return 0.996788;
    else if (amplitude_multifit_intime_uncal-4400.000000<200) return 1.000612;
    else if (amplitude_multifit_intime_uncal-4800.000000<200) return 1.010879;
    else if (amplitude_multifit_intime_uncal-5200.000000<200) return 1.025890;
    else if (amplitude_multifit_intime_uncal-5600.000000<200) return 1.035207;
    else if (amplitude_multifit_intime_uncal-6000.000000<200) return 1.047249;
    else if (amplitude_multifit_intime_uncal-6400.000000<200) return 1.060303;
    else if (amplitude_multifit_intime_uncal-6800.000000<200) return 1.070912;
    else if (amplitude_multifit_intime_uncal-7200.000000<200) return 1.079619;
    else if (amplitude_multifit_intime_uncal-7600.000000<200) return 1.086805;
    else if (amplitude_multifit_intime_uncal-8000.000000<200) return 1.090260;
    else if (amplitude_multifit_intime_uncal-8400.000000<200) return 1.066447;
    else if (amplitude_multifit_intime_uncal-8800.000000<200) return 1.041410;
    else if (amplitude_multifit_intime_uncal-9200.000000<200) return 1.024488;
    else if (amplitude_multifit_intime_uncal-9600.000000<200) return 1.014497;
    else if (amplitude_multifit_intime_uncal-10000.000000<200) return 1.016732;
    else if (amplitude_multifit_intime_uncal-10400.000000<200) return 1.029175;
    else if (amplitude_multifit_intime_uncal-10800.000000<200) return 1.032594;
    else if (amplitude_multifit_intime_uncal-11200.000000<200) return 1.038157;
    else if (amplitude_multifit_intime_uncal-11600.000000<200) return 1.045066;
    else return 1.039888;
      }
  else {
      if (amplitude_multifit_intime_uncal-2800.000000<200) return 1.000000;
      else if (amplitude_multifit_intime_uncal-3200.000000<200) return 0.999389;
      else if (amplitude_multifit_intime_uncal-3600.000000<200) return 0.998585;
      else if (amplitude_multifit_intime_uncal-4000.000000<200) return 1.006239;
      else if (amplitude_multifit_intime_uncal-4400.000000<200) return 1.008576;
      else if (amplitude_multifit_intime_uncal-4800.000000<200) return 1.009653;
      else if (amplitude_multifit_intime_uncal-5200.000000<200) return 1.010044;
      else if (amplitude_multifit_intime_uncal-5600.000000<200) return 1.009461;
      else if (amplitude_multifit_intime_uncal-6000.000000<200) return 1.009145;
      else if (amplitude_multifit_intime_uncal-6400.000000<200) return 1.009076;
      else if (amplitude_multifit_intime_uncal-6800.000000<200) return 1.009964;
      else if (amplitude_multifit_intime_uncal-7200.000000<200) return 1.008127;
      else if (amplitude_multifit_intime_uncal-7600.000000<200) return 1.008139;
      else if (amplitude_multifit_intime_uncal-8000.000000<200) return 1.009407;
      else if (amplitude_multifit_intime_uncal-8400.000000<200) return 1.009913;
      else if (amplitude_multifit_intime_uncal-8800.000000<200) return 1.008228;
      else if (amplitude_multifit_intime_uncal-9200.000000<200) return 1.007341;
      else if (amplitude_multifit_intime_uncal-9600.000000<200) return 1.006882;
      else if (amplitude_multifit_intime_uncal-10000.000000<200) return 1.007504;
      else if (amplitude_multifit_intime_uncal-10400.000000<200) return 1.006405;
      else if (amplitude_multifit_intime_uncal-10800.000000<200) return 1.006475;
      else if (amplitude_multifit_intime_uncal-11200.000000<200) return 1.004366;
      else if (amplitude_multifit_intime_uncal-11600.000000<200) return 1.004993;
      else return 1.004440;
  }

}

// ------------ method called once each job just before starting event loop  ------------
void 
PulseTree::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
PulseTree::endJob() 
{

  if (do_average) WriteAverageOutput();

}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
PulseTree::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(PulseTree);
