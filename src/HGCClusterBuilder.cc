#include "RecoParticleFlow/HGCTracking/interface/HGCClusterBuilder.h"
#include "FWCore/Utilities/interface/Exception.h"

/*
HGCClusterBuilder::HGCClusterBuilder():
  energyTot_m(-1.), sizeTot_m(-1), layer_m(-1), 
  seedL_m(-1), nUnits_m(-1), energyFromSeed_m(-1)
{ }
*/

HGCClusterBuilder::HGCClusterBuilder():
  is2DClSeed(0), seedEnergy(0), seedSize(0), seedLayer(0),
  nValidMeas(0),
  Esum(0), Ssum(0), Lsum(0),
  Emean(0), Smean(0), Lmean(0),
  aE(0), bE(0), aS(0), bS(0){
  energyL.fill(0.);
  sizeL.fill(0);
}


void HGCClusterBuilder::initState(float energy, float size, int L, bool kSeed){ 
  is2DClSeed = kSeed;

  seedEnergy = energy;
  seedLayer = L;
  seedSize = size;

  if(kSeed == false) return;

  // to cross-check if gets double counted
  energyL[L-1] = energy;
  sizeL[L-1] = size; 
  nValidMeas.push_back(L);

  Esum += energy;
  Ssum += size;
  Lsum += L;

  nValidMeas.push_back(L);
}

/*
float HGCClusterBuilder::getEvoCumulE(int layer ){
  float totE = 0;
  for(int ij=seedL_m-1; ij<layer; ++ij) totE += energyTot_m[ij];
  return totE;
}
*/


void HGCClusterBuilder::evolveState(float energy, float size, int layer){
  if(energy == 0) return;

  if(layer == seedLayer){
    //do not increase if filling same layer with different starting point
    if(std::find(nValidMeas.begin(), nValidMeas.end(), layer) == nValidMeas.end()){
      nValidMeas.push_back(layer);
      Lsum += layer;
    }
    
    Esum += energy;
    Ssum += size;
    
    energyL[layer-1] += energy;
    sizeL[layer-1] += size;
  }

  //std::cout << " Esum = " << Esum << " energyL[layer-1] = " << energyL[layer-1] << " nValidMeas.size() = " << nValidMeas.size()  << std::endl;
  //std::cout << " >>> now updateEvolution " << std::endl;
  updateEvolution();
}


void HGCClusterBuilder::updateEvolution(){
  if(nValidMeas.size() < 3) return;

  float OneoverN = 1./nValidMeas.size();

  Emean = Esum*OneoverN;
  Smean = Ssum*OneoverN;
  Lmean = Lsum*OneoverN;
  float xYE = 0;
  float xYS = 0;
  float xxL = 0;

  for(unsigned int it = nValidMeas.size()-3; it < nValidMeas.size(); ++it){
    int lVal = nValidMeas.at(it);
    xYE += (energyL[lVal] - Emean) * (lVal - Lmean);
    xYS += (sizeL[lVal] - Smean) * (lVal - Lmean);
    xxL += (lVal - Lmean)*(lVal - Lmean);
  }

  xYE = xYE*OneoverN;
  xYS = xYS*OneoverN;
  xxL = xxL*OneoverN;

  bE = xYE/xxL;
  aE = Emean - bE * Lmean;
  bS = xYS/xxL;
  aS = Smean - bS * Lmean;
}

std::pair<float, float> HGCClusterBuilder::evaluateESChi2(float sumEL, float sumSL, int layer, int originL, int originSize){

  std::cout << " in HGCClusterBuilder::evaluateESChi2 layer = " << layer << " nValidMeas.size() = " << nValidMeas.size()
            << "  energyL[layer-2] = " << energyL[layer-2] << " energyL[layer-1] = " <<energyL[layer-1]
            << " evoState->getSeedEnergy() = " << seedEnergy << std::endl;


  if(nValidMeas.size() == 0) return std::pair<float, float>(0, 0);
  else if(nValidMeas.size() < 3) return std::pair<float, float>(energyL[layer-2], sizeL[layer-2]);
  //  else if(nValidMeas.size() < 3) return std::pair<float, float>(evoState->getSeedEnergy(), evoState->getSeedSize());
  float expE = aE + bE * layer;
  float expS = aS + bS * layer;

  std::cout << " >>> expE = " << expE  << " expS = " << expS << std::endl;


  return  std::pair<float, float>(expE, expS);


  /*
  if(Nlayers < 2) {
    //    std::cout << " first layers -- " << std::endl;
    float expE = (sumEL - energyL[seedLayer-1] ) / energyL[seedL_m-1];
    float expS = (sumSL - sizeL[seedLayer-1] ) / sizeL[seedL_m-1];
    return expE*expE + expS*expS;
  }

  float expE = aE + bE * layer;
  float expS = aS + bS * layer;

  if((expE < 0 || expS < 0) && energyL[originL-1] > 0){
    float expE = (sumEL - energyL[originL-1] ) / energyTot_m[originL-1];
    float expS = (sumSL - sizeL[originL-1] ) / sizeTot_m[originL-1];
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

 
  // if(energyTot_m[layer-1] == 0. || sizeTot_m[layer-1] == 0) {
  //   throw cms::Exception("Configuration") << " Empty starting state for HGCClusterBuilder energyTot_m = " 
  // 					  << energyTot_m[layer-1] << " sizeTot_m = " << sizeTot_m[layer-1] << " \n ";
  // }
 

  // float Echi2 = 1.;
  //  float Echi2 = (nUnits_m[prevLayer-1] != 0) ? (expectedE - sumEL)/nUnits_m[prevLayer-1] : (expectedE - sumEL)/sumSL;
  //  float Schi2 = (sizeTot_m[layer-1] - sumSL) / sizeTot_m[layer-1];
  //float totalChi = (layer < 35 || layer - seedL_m < 20)? Echi2+Schi2 : std::abs(Echi2)+std::abs(Schi2);
  //  float totalChi = (layer < 35 || layer - seedL_m < 20)? Echi2+Schi2 : std::abs(Echi2)+std::abs(Schi2);
  //return totalChi;
  return Echi2*Echi2 +Schi2*Schi2;
  */
}

