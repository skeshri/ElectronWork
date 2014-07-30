#ifndef VARIABLELIMITS_HH
#define VARIABLELIMITS_HH

#include "TString.h"

namespace VarLims{

  const int nVarLimits = 9;
  
  struct VariableLimits {
    TString name;
    float   max;
    VariableLimits(TString nameIn, float maxIn):
      name(nameIn), max(maxIn){};
  };

  // For the veto working point, do not restrict anything
  VariableLimits * limitsNoRestrictions[nVarLimits] = {
    new VariableLimits("full5x5_sigmaIetaIeta"   ,1e30),
    new VariableLimits("dEtaIn"                  ,1e30),
    new VariableLimits("dPhiIn"                  ,1e30),
    new VariableLimits("hOverE"                  ,1e30),
    new VariableLimits("relIsoWithDBeta"         ,1e30),
    new VariableLimits("ooEmooP"                 ,1e30),
    new VariableLimits("d0"                      ,1e30),
    new VariableLimits("dz"                      ,1e30),
    new VariableLimits("expectedMissingInnerHits",1e30)
  };

  // For Loose, Medium, Tight working points:
  //   restrict d0 (at a very loose value)
  //   restrict missing hits, to Run 1-motivated value
  //      (to make sure ID discriminates well against 
  //       hard photons faking electrons, which may not be
  //       plentiful in the sample on which we are optimizing).
  //   restrict H/E, dPhiIn: all are very loose common sense
  //       restrictions.
  VariableLimits * limitsWPAnyV1[nVarLimits] = {
    new VariableLimits("full5x5_sigmaIetaIeta"   ,1e30),
    new VariableLimits("dEtaIn"                  ,0.04), //eff EB(EE)100.0(99.9)%
    new VariableLimits("dPhiIn"                  ,0.15), //eff EB(EE) 99.5(99.6)%
    new VariableLimits("hOverE"                  ,0.15), //eff EB(EE) 99.4(99.3)%
    new VariableLimits("relIsoWithDBeta"         ,1e30),
    new VariableLimits("ooEmooP"                 ,1e30),
    new VariableLimits("d0"                      ,0.10), //eff EB(EE) 99.9(99.2)% 
    new VariableLimits("dz"                      ,1e30),
    new VariableLimits("expectedMissingInnerHits",1e30)  //eff EB(EE) 99.8(99.4)%
  };

}

#endif
