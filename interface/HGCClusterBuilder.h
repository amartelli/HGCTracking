#ifndef RecoparticleFlow_HGCTracking_HGCClusterBuilder_h
#define RecoparticleFlow_HGCTracking_HGCClusterBuilder_h


#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include <cmath>
#include <cstdlib>
#include <iostream>

/// Class that does the trajectory building starting from a seed


class HGCClusterBuilder {
    public:
  HGCClusterBuilder();
  ~HGCClusterBuilder() {};

  void initState(float energy, float size, int layer, bool is2DClSeed);
  inline bool is2DClusterSeed() {return is2DClSeed; };


  inline const float getSeedLayer(){return seedLayer;};
  inline const float getSeedEnergy(){return seedEnergy;};
  inline const float getSeedSize(){return seedSize;};

  void evolveState(float energy, float size, int layer); 
  void updateEvolution();
  std::pair<float, float> evaluateESChi2(float sumEL, float sumSL, int L, int oroginL, int originSize);  
  
  inline float getTotalStateEnergyL(int layer){return energyL[layer-1];};

 private:
  bool is2DClSeed;

  float seedEnergy;
  float seedSize;
  float seedLayer;

  std::vector<int> nValidMeas;

  /* //to fit evolution */
  float Esum;
  float Ssum;
  float Lsum;

  float Emean; 
  float Smean; 
  float Lmean; 

  float aE; 
  float bE;
  float aS;
  float bS;
  
  std::array<float, 52> energyL;
  std::array<float, 52> sizeL;
};

#endif
