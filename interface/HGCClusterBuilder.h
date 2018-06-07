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



//monitor evolution
/*
struct evolutionState{ 
float energyTot; 
float sizeTot; 
//LocalPoint avgPos; 
int layer; 
int nUnits; 
} evoState; 
*/

class HGCClusterBuilder {
    public:
  HGCClusterBuilder();
  ~HGCClusterBuilder() {};

  void initState(float energy, float size, int layer);
  
  inline void updateLayer(int layer) {layer_m = layer;};
  void evolveState(float energy, float size, int layer, unsigned int nUnits);
  void updateEvolution(int layer);

  float evaluateESChi2(float sumEL, float sumSL, int L, int oroginL, int originSize);

  
  inline float getEvoETot(int layer){return energyTot_m[layer-1];};
  inline float getEvoSTot(int layer ){return sizeTot_m[layer-1];};
  inline float getEvoL(){return layer_m;};
  inline float getEvoSeedL(){return seedL_m;};
  inline float getEvoUnits(int layer){return nUnits_m[layer-1];};

  float getEvoCumulE(int layer );
  
  /*
  inline float getEvoETot(){return energyTot_m;};
  inline float getEvoSTot(){return sizeTot_m;}; 
  inline float getEvoL(){return layer_m;};                       
  inline float getEvoSeedL(){return seedL_m;};                   
  inline float getEvoUnits(){return nUnits_m;};  
  inline float getEvoCumulE(){return energyFromSeed_m;};         
  */

 private:

  //to fit evolution
  float Lsum;
  float Esum;
  float Ssum;
  int Nlayers;
  float Emean;
  float Smean;
  float Lmean;

  float aE;
  float bE;
  float aS;
  float bS;


  
  std::array<float, 52> energyTot_m;            
  std::array<float, 52> sizeTot_m;              
  //LocalPoint avgPos;        
  std::array<int, 52> layerOK_m;
  int layer_m;                  
  std::array<unsigned int, 52> nUnits_m;                 
  int seedL_m;
  std::array<float,52> energyFromSeed_m;            



  /*
  float energyTot_m;     
  float sizeTot_m;       
  //LocalPoint avgPos;       
  int layer_m;           
  int seedL_m;               
  unsigned int nUnits_m; 
  float energyFromSeed_m;    
  */
};

#endif
