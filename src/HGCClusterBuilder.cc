#include "RecoParticleFlow/HGCTracking/interface/HGCClusterBuilder.h"
#include "FWCore/Utilities/interface/Exception.h"

/*
HGCClusterBuilder::HGCClusterBuilder():
  energyTot_m(-1.), sizeTot_m(-1), layer_m(-1), 
  seedL_m(-1), nUnits_m(-1), energyFromSeed_m(-1)
{ }
*/
HGCClusterBuilder::HGCClusterBuilder():
  Lsum(0), Esum(0), Ssum(0), Nlayers(0), 
  Emean(0), Smean(0), Lmean(0), 
  aE(0), bE(0), aS(0), bS(0){ 
  energyTot_m.fill(0.);
  sizeTot_m.fill(0);
  layerOK_m.fill(0);
  layer_m = 0;
  seedL_m = 0; 
  nUnits_m.fill(0); 
  energyFromSeed_m.fill(0);
}


void HGCClusterBuilder::initState(float energy, float size, int L){
  
  energyTot_m[L-1] = energy;
  sizeTot_m[L-1] = size; 
  layerOK_m[L-1] = 1;
  layer_m = L;
  seedL_m = L;
  nUnits_m[L-1] = 1;
  energyFromSeed_m[L-1] = energy;

  Nlayers += 1;
  Esum += energy;
  Lsum += L;
  Ssum += size;
  /*
  energyTot_m = energy;
  sizeTot_m = size; 
  layer_m = L;
  seedL_m = L;
  nUnits_m = 1;
  energyFromSeed_m = energy;
  */

  //std::cout << " >>  initState energyTot_m = " << energyTot_m  << std::endl;
}


float HGCClusterBuilder::getEvoCumulE(int layer ){
  float totE = 0;
  for(int ij=seedL_m-1; ij<layer; ++ij) totE += energyTot_m[ij];
  return totE;
}


void HGCClusterBuilder::evolveState(float energy, float size, int layer, unsigned int nUnits){
  //  std::cout << " evolveState startE = " << energyTot_m << std::endl;
  if(energy == 0) return;

  //do not increase if filling same layer with different starting point
  if(layerOK_m[layer-1] == 0){
    Nlayers += 1;
    Lsum += layer;
  }
  Esum += energy;
  Ssum += size;

  energyTot_m[layer-1] += energy;
  sizeTot_m[layer-1] += size;
  layerOK_m[layer-1] = 1;
  layer_m = layer;
  nUnits_m[layer-1] += nUnits;
  for(int ij=seedL_m-1; ij<layer; ++ij) energyFromSeed_m[ij] += energy;

  updateEvolution(layer);

  //  std::cout << " evolveState endE = " << energyTot_m << std::endl;
}


void HGCClusterBuilder::updateEvolution(int layer){

  //  std::cout << " updateEvolution layers filled = " << layer - seedL_m << std::endl;
  if(Nlayers < 2) return;

  float OneoverN = 1./Nlayers;

  Emean = Esum*OneoverN;
  Smean = Ssum*OneoverN;
  Lmean = Lsum*OneoverN;

  float xYE = 0;  
  float xYS = 0;
  float xxL = 0;

  for(int ij=seedL_m; ij<=layer; ++ij){
    if(layerOK_m[ij-1] == 0) continue;
    //    std::cout << " layer ok = " << ij << std::endl;
    xYE += (energyTot_m[ij-1] - Emean) * (ij - Lmean); 
    xYS += (sizeTot_m[ij-1] - Smean) * (ij - Lmean); 
    xxL += (ij - Lmean)*(ij - Lmean);
  }

  xYE = xYE*OneoverN;
  xYS = xYS*OneoverN;  
  xxL = xxL*OneoverN;

  bE = xYE/xxL;
  aE = Emean - bE * Lmean;

  bS = xYS/xxL;
  aS = Smean - bS * Lmean;

  // std::cout << " Emean = " << Emean << " Smean = " << Smean << " Lmean = " << Lmean << std::endl;
  // std::cout << " aE = " << aE << " bE = " << bE << " aS = " << aS << " bS = " << bS << std::endl;

};



float HGCClusterBuilder::evaluateESChi2(float sumEL, float sumSL, int layer, int originL, int originSize){

  if(Nlayers < 2) {
    //    std::cout << " first layers -- " << std::endl;
    float expE = (sumEL - energyTot_m[seedL_m-1] ) / energyTot_m[seedL_m-1];
    float expS = (sumSL - sizeTot_m[seedL_m-1] ) / sizeTot_m[seedL_m-1];
    return expE*expE + expS*expS;
  }

  float expE = aE + bE * layer;
  float expS = aS + bS * layer;

  if((expE < 0 || expS < 0) && energyTot_m[originL-1] > 0){
    float expE = (sumEL - energyTot_m[originL-1] ) / energyTot_m[originL-1];
    float expS = (sumSL - sizeTot_m[originL-1] ) / sizeTot_m[originL-1];
    //    std::cout << " exp neg before ok " << std::endl;
    return expE*expE + expS*expS;
  }
  else if((expE < 0 || expS < 0)){
    float totE = getEvoCumulE(layer) / Nlayers;
    float expE = (sumEL - totE ) / totE;
    //    std::cout << " exp neg before NO " << std::endl;
    return expE*expE;
  }

  float rescaleFactor =(nUnits_m[originL-1] > 0) ? originSize/sizeTot_m[originL-1] : 1;
  //start from invalid hits
  if(originSize == -1) rescaleFactor = 1;

  //  std::cout << " >>> expected energy = " << expE << " exp Size = " << expS << "current L = " << originL << " scaleF = " << rescaleFactor << std::endl;

  float Echi2 = (sumEL - rescaleFactor * expE)/expE;
  float Schi2 = (sumSL - rescaleFactor * expS)/expS;


  //  std::cout << " >>> Echi2 = " << Echi2 << " Schi2 = " << Schi2 << " quad sum = " << Echi2*Echi2 +Schi2*Schi2 << std::endl;

  /*
  if(energyTot_m[layer-1] == 0. || sizeTot_m[layer-1] == 0) {
    throw cms::Exception("Configuration") << " Empty starting state for HGCClusterBuilder energyTot_m = " 
					  << energyTot_m[layer-1] << " sizeTot_m = " << sizeTot_m[layer-1] << " \n ";
  }
  */

  // float Echi2 = 1.;
  //  float Echi2 = (nUnits_m[prevLayer-1] != 0) ? (expectedE - sumEL)/nUnits_m[prevLayer-1] : (expectedE - sumEL)/sumSL;
  //  float Schi2 = (sizeTot_m[layer-1] - sumSL) / sizeTot_m[layer-1];
  //float totalChi = (layer < 35 || layer - seedL_m < 20)? Echi2+Schi2 : std::abs(Echi2)+std::abs(Schi2);
  //  float totalChi = (layer < 35 || layer - seedL_m < 20)? Echi2+Schi2 : std::abs(Echi2)+std::abs(Schi2);
  //return totalChi;
  return Echi2*Echi2 +Schi2*Schi2;

}
