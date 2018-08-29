#include "RecoParticleFlow/HGCTracking/interface/TempCaloTrajectory.h"
#include "FWCore/Utilities/interface/Exception.h"

TempCaloTrajectory::TempCaloTrajectory():
  nValidMeas(0),
  Esum(0), Ssum(0), Lsum(0), 
  Emean(0), Smean(0), Lmean(0), 
  aE(0), bE(0), aS(0), bS(0){ 
  energyL.fill(0.);
  sizeL.fill(0);
}


TempCaloTrajectory::TempCaloTrajectory(PropagationDirection direction, unsigned char nhseed):
  tmpTraj(TempTrajectory(direction, nhseed)), 
  nValidMeas(0), 
  Esum(0), Ssum(0), Lsum(0),
  Emean(0), Smean(0), Lmean(0), 
  aE(0), bE(0), aS(0), bS(0){ 
  energyL.fill(0.);
  sizeL.fill(0);
}


TempCaloTrajectory::TempCaloTrajectory(TempTrajectory traj, TempCaloTrajectory caloTraj){
  using std::swap;
  swap(nValidMeas, caloTraj.nValidMeas);
  swap(Esum, caloTraj.Esum);
  swap(Ssum, caloTraj.Ssum);
  swap(Lsum, caloTraj.Lsum);
  swap(Emean, caloTraj.Emean);
  swap(Smean, caloTraj.Smean);
  swap(Lmean, caloTraj.Lmean);
  swap(aE, caloTraj.aE); swap(bE, caloTraj.bE); swap(aS, caloTraj.aS); swap(bS, caloTraj.bS);
  swap(energyL, caloTraj.energyL);
  swap(sizeL, caloTraj.sizeL);
  swap(tmpTraj, traj);
}


std::pair<float, float> TempCaloTrajectory::getPredState(HGCClusterBuilder* evoState, int layer) const{

  std::cout << " in TempCaloTrajectory::getPredState layer = " << layer << " nValidMeas.size() = " << nValidMeas.size()
	    << " energyL[layer-1] = " << energyL[layer-1]
	    << " evoState->getSeedLayer() = " << evoState->getSeedLayer() 
	    << " evoState->getSeedEnergy() = " << evoState->getSeedEnergy() << std::endl;
  
  if(nValidMeas.size() == 0) return std::pair<float, float>(0, 0);
  else if(nValidMeas.size() < 3) return std::pair<float, float>(energyL[layer-2], sizeL[layer-2]);
  //  else if(nValidMeas.size() < 3) return std::pair<float, float>(evoState->getSeedEnergy(), evoState->getSeedSize());


  float expE = aE + bE * layer;
  float expS = aS + bS * layer;

  if(expE <= 0 ) expE = energyL[layer-2];
  if(expS <= 0 ) expS = sizeL[layer-2];

  std::cout << " >>> expE = " << expE  << " expS = " << expS << std::endl;

  return  std::pair<float, float>(expE, expS);
}

void TempCaloTrajectory::updateState(float energy, float size, int layer){
  if(energy == 0) return;

  //do not increase if filling same layer with different starting point
  if(std::find(nValidMeas.begin(), nValidMeas.end(), layer) == nValidMeas.end()){
    nValidMeas.push_back(layer);
    Lsum += layer;
  }

  Esum += energy;
  Ssum += size;

  energyL[layer-1] += energy;
  sizeL[layer-1] += size;
  
  if(nValidMeas.size() > 3){
    Esum = 0;
    Ssum = 0;
    Lsum = 0;
    for(unsigned int it = nValidMeas.size()-3; it < nValidMeas.size(); ++it){
      int lVal = nValidMeas.at(it);
      Esum += energyL[lVal-1];
      Ssum += sizeL[lVal-1];
      Lsum += lVal;
    }
  }
  //std::cout << " >>> TempCaloTrajectory::updateState => energyL[layer-1] =  " << energyL[layer-1] << std::endl;
  //for(int ij=seedL_m-1; ij<layer; ++ij) energyFromSeed_m[ij] += energy;
  
  // std::cout << " now move to updateEvolution " << std::endl;
  updateEvolution();
}


void TempCaloTrajectory::updateEvolution(){

  //std::cout << " in TempCaloTrajectory::updateEvolution nValidMeas.size() = " << nValidMeas.size() << std::endl;

  if(nValidMeas.size() < 3) return;

  float OneoverN = 1./nValidMeas.size();

  // Emean = Esum*OneoverN;
  // Smean = Ssum*OneoverN;
  // Lmean = Lsum*OneoverN;


  //just interpolate the last 3 measurements
  // to allow change in the trend (derivative => increase max decrease)
  Emean = Esum/3.;
  Smean = Ssum/3.;
  Lmean = Lsum/3.;

  float xYE = 0;
  float xYS = 0;
  float xxL = 0;

  for(unsigned int it = nValidMeas.size()-3; it < nValidMeas.size(); ++it){
    int lVal = nValidMeas.at(it);
    xYE += (energyL[lVal-1] - Emean) * (lVal - Lmean);
    xYS += (sizeL[lVal-1] - Smean) * (lVal - Lmean);
    xxL += (lVal - Lmean)*(lVal - Lmean);

    // std::cout << " it = " << it << " energyL[lVal-1] = " << energyL[lVal-1] << " sizeL[lVal-1] = " << sizeL[lVal-1] << " lVal = " << lVal << std::endl;
    // std::cout << " xYE = " << xYE << " xYS = " << xYS << " xxL = " << xxL << std::endl;
  }

  /*
  xYE = xYE*OneoverN;
  xYS = xYS*OneoverN;
  xxL = xxL*OneoverN;
  */

  xYE = xYE/3.;
  xYS = xYS/3.;
  xxL = xxL/3.;


  //std::cout << " xYE = " << xYE << " xYS = " << xYS << " xxL = " << xxL << std::endl;

  bE = xYE/xxL;
  aE = Emean - bE * Lmean;
  bS = xYS/xxL;
  aS = Smean - bS * Lmean;


  //std::cout << " bE = " << bE << " aE = " << aE << " bS = " << bS << " aS = " << aS << std::endl;

}




