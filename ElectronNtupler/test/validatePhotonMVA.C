{

  gROOT->ProcessLine(".L computePhotonIDEfficiencyMVA.C+");

  // Compute photon ID efficiency using the MVA cut corresponding to 90%
  // signal efficiency for pt ranges 40-60 and 60+ for barrel and endcap,
  // as done on slide 12 in
  //    https://indico.cern.ch/event/369185/contribution/2/material/slides/0.pdf

  computePhotonIDEfficiencyMVA(0,0);
  computePhotonIDEfficiencyMVA(0,1);
  computePhotonIDEfficiencyMVA(1,0);
  computePhotonIDEfficiencyMVA(1,1);


}
