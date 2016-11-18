#ifndef __RecHitSampleFourCorrector_h__
#define __RecHitSampleFourCorrector_h__

#include "DataFormats/EcalRecHit/interface/EcalUncalibratedRecHit.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHit.h"

class RecHitSampleFourCorrector {

 public:
  RecHitSampleFourCorrector();
  ~RecHitSampleFourCorrector();

  float RecHitCorrectedEnergy(const EcalRecHit *rh, float lc, float ic, float agc);
  float MultiFitParametricCorrection(float amplitude_multifit_intime_uncal, float chi2, uint32_t recoflag);

 private:
  double CorrectionFunction1(double x, double chi2);
  double CorrectionFunction2(double x, double chi2);
  double CorrectionFunction3(double x, double chi2);

};

#endif
