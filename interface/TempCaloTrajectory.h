#ifndef RecoparticleFlow_HGCTracking_TempCaloTrajectory_h
#define RecoparticleFlow_HGCTracking_TempCaloTrajectory_h


#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "TrackingTools/PatternTools/interface/TempTrajectory.h"
#include "RecoParticleFlow/HGCTracking/interface/HGCClusterBuilder.h"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <utility>
/// Class that does the trajectory building starting from a seed

//class TempTrajectory;


class TempCaloTrajectory{
 public:
  TempCaloTrajectory();
  TempCaloTrajectory(PropagationDirection direction, unsigned char nhseed);
  TempCaloTrajectory(const TempCaloTrajectory& rh) = default;

  TempCaloTrajectory & operator=(TempCaloTrajectory const & rh) = default;
  
  TempCaloTrajectory(TempCaloTrajectory && rh)  noexcept: 
  tmpTraj(std::move(rh.tmpTraj)),
    nValidMeas(rh.nValidMeas),
    Esum(rh.Esum), 
    Ssum(rh.Ssum), 
    Lsum(rh.Lsum), 
    Emean(rh.Emean), 
    Smean(rh.Smean), 
    Lmean(rh.Lmean), 
    aE(rh.aE), 
    bE(rh.bE), 
    aS(rh.aS), 
    bS(rh.bS), 
    energyL(rh.energyL), 
    sizeL(rh.sizeL) 
      { 
	rh.nValidMeas.clear();
	rh.Esum = 0.;
	rh.Ssum = 0.;
	rh.Lsum = 0.;
	rh.Emean = 0.;
	rh.Smean = 0.;
	rh.Lmean = 0.;
	rh.aE = 0.;
	rh.bE = 0.;
	rh.aS = 0.;
	rh.bS = 0.;
	rh.energyL.fill(0);
	rh.sizeL.fill(0);
      }

  TempCaloTrajectory & operator=(TempCaloTrajectory && rh) noexcept { 
    using std::swap;
    swap(tmpTraj, rh.tmpTraj); 
    swap(nValidMeas, rh.nValidMeas); 
    swap(Esum, rh.Esum); 
    swap(Ssum, rh.Ssum); 
    swap(Lsum, rh.Lsum); 
    swap(Emean, rh.Emean); 
    swap(Smean, rh.Smean); 
    swap(Lmean, rh.Lmean); 
    swap(aE, rh.aE); swap(bE, rh.bE); swap(aS, rh.aS); swap(bS, rh.bS); 
    swap(energyL, rh.energyL); 
    swap(sizeL, rh.sizeL); 
    return *this; 
   } 


  TempCaloTrajectory(TempTrajectory, TempCaloTrajectory);
  ~TempCaloTrajectory(){};
  
  
  inline TempTrajectory* getTempTraj(){ /*std::cout << " called getTempTraj() " << std::endl;*/ return  &tmpTraj; };

  std::pair<float, float> getPredState(HGCClusterBuilder* evoState, int layer) const;
  void updateState(float energy, float size, int layer);
  void updateEvolution();

  inline bool isValid() const {return tmpTraj.isValid(); };
  inline int foundHits() const {return tmpTraj.foundHits(); };
  inline int lostHits() const {return tmpTraj.lostHits(); };
  inline float chiSquared() const {return tmpTraj.chiSquared();};
  inline const TempTrajectory::DataContainer & measurements() const {return tmpTraj.measurements(); };
  inline PropagationDirection direction() const {return tmpTraj.direction(); };
  inline void invalidate() {tmpTraj.invalidate(); return;};
  inline void pop() {tmpTraj.pop(); return;};
  inline void addMeasurement(TrajectoryMeasurement trjM) {tmpTraj.push(trjM); return;};
  inline void addMeasurement(TrajectoryMeasurement trjM, float chi2) {tmpTraj.push(trjM, chi2); return;};


  inline float getCurrentE(int layer){return energyL[layer-1];};
  inline float getCurrentS(int layer){return sizeL[layer-1];}; 
  

 private:
  TempTrajectory tmpTraj;

  //to fit evolution
  std::vector<int> nValidMeas;

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
