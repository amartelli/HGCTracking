#include "RecoParticleFlow/HGCTracking/interface/TrajectoryCleanerByBranchFilter.h"
#include "TrackingTools/PatternTools/interface/TrajectoryMeasurement.h"
#include "DataFormats/HGCRecHit/interface/HGCRecHitCollections.h"

//logic:
// share if >= 3 among first 5 hits are shared
// if track has less than 5 hits reuire that fraction of shared energy > 1.5 equalsharing and that max 1 hit is on shared
// need option for tracks sharing fraction of energy in the bulk non at the beginning
// => try with mergin if they share > 80%

//template<typename Traj>
std::map<unsigned int, std::vector<TrajectoryMeasurement > > TrajectoryCleanerByBranchFilter::mergeTraj(const std::vector<Trajectory> &trajs)
{

  bool printInfo = false;
  bool printInfo2 = false;
  bool printInfo3 = false;
  int totLoop = 1;

  //not used - try clever criteria...
  unsigned int shareFirstN = 3;
  unsigned int checkFirstN = 5;

  //map to record the shared tracks - throughout common measurements - 
  //beging with each track is shared with itself
  std::map<unsigned int, std::set<unsigned int> > firstShareWithSecond;
  for(unsigned int ij=0; ij<trajs.size(); ++ij){
    std::set<unsigned int> dummy;
    dummy.insert(ij);
    firstShareWithSecond[ij] = dummy;
  }
  std::map<unsigned int, std::vector<TrajectoryMeasurement> > allHits;
  std::map<unsigned int, std::vector<TrajectoryMeasurement> > allHits_Init;
  std::map<unsigned int, std::vector<float > > allHitsEnergyFr_Init;
  std::map<unsigned int, std::vector<float > > allHitsEnergyFr;

  std::set<reco::CaloClusterPtr> outCl;
  std::set<const HGCRecHit*> outRh;   


  unsigned int trajSize = trajs.size();

  std::vector<std::pair<int, int> >nTrkSize;

  if(printInfo)  std::cout << " total trk size = " << trajSize << std::endl;
  for (unsigned int i = 0; i < trajSize; ++i) {
    const Trajectory &t1 = trajs[i];
    auto first = t1.measurements();
    
    if(printInfo) std::cout << " track " << i << " size meas = " << first.size() << std::endl;
    
    float trackSumEnergy = 0.;
    //initialize possible non-shared tracks
    // map[track_index] = vector<measurements>
    for(auto tm1 : first){
   
      float energy = 0.;
      if (typeid(*tm1.recHit()) == typeid(HGCTrackingRecHitFromCluster)) {
	reco::CaloClusterPtr cl = (dynamic_cast<const HGCTrackingRecHitFromCluster&>(*tm1.recHit())).objRef();
	energy = cl->energy();
      }
      else if (typeid(*tm1.recHit()) == typeid(HGCTrackingClusteringRecHit)) {
	//should implement case for this
      }
      else if (typeid(*tm1.recHit()) == typeid(HGCTrackingRecHitFromHit)) {
	const HGCRecHitRef rh = (dynamic_cast<const HGCTrackingRecHitFromHit&>(*tm1.recHit()) ).objRef();
	energy = rh->energy();
      }
      else{
	//case of invalid hits: ok during trajectory building but useless when filling collection
	continue;
	//std::cout << "Invalid valid hit" << typeid(*tm1.recHit()).name() << " is valid = " << (*tm1.recHit()).isValid() << std::endl;
      }
      
      allHits_Init[i].push_back(tm1);
      
      //how much each measurement contribute to the energy of the track
      if(printInfo) std::cout << " filling map with energy = " << energy << " and trackSumEnergy = " << trackSumEnergy << std::endl;
      allHitsEnergyFr_Init[i].push_back(energy);
      trackSumEnergy += energy;
    }
    
    nTrkSize.push_back(std::pair<int, int>(i, allHits_Init[i].size()));    
    
    if(trackSumEnergy > 0.) for(unsigned int iTM=0; iTM < allHitsEnergyFr_Init[i].size(); ++iTM) {
	allHitsEnergyFr_Init[i][iTM] /= trackSumEnergy;
	if(printInfo) std::cout << " fraction = " << iTM << " mappa = " << allHitsEnergyFr_Init[i].at(iTM) << std::endl;
      }
  }
  
  std::cout << " map allHits_Init.size = " << allHits_Init.size() << " traj size = " << trajSize << std::endl; 
  if(allHits_Init.size() != trajSize) std::cout << " size problem " << std::endl;
  
  
  //longer first => will copy less objects later
  auto sortTrajs = [&](const std::pair<int, int>& a, const std::pair<int, int>& b) 
    {
      return (a.second > b.second);
    };
  
  std::sort(nTrkSize.begin(), nTrkSize.end(), sortTrajs);
  
  for(auto itMap=0; itMap<int(allHits_Init.size()); ++itMap ){
    allHits[itMap] = allHits_Init[nTrkSize.at(itMap).first];
    allHitsEnergyFr[itMap] = allHitsEnergyFr_Init[nTrkSize.at(itMap).first];
  }

  std::cout << " allHits.size = " << allHits.size() << std::endl;
  for(auto ip:allHits) std::cout << "trj.size = " << ip.second.size() << std::endl;
  
  std::cout << " fine per ora " << std::endl;
  
  //iterations might help? not clear for the moment...
  for(int nLoop=0; nLoop<totLoop; ++nLoop){

    if(printInfo2) std::cout << " pre allHits.size = " << allHits.size() << std::endl;

    //check for sharing 
    //if shared: copy non-common measurements in longer track
    //flag the sharing of the 2 tracks
    //cancel all measurements from shorter track (they are already all assigned to the longer)
    //clean map
    int count1 = -1;
    for (std::map<unsigned int, std::vector<TrajectoryMeasurement> >::iterator first = allHits.begin(); first != allHits.end(); ++first) {
      ++count1;
      int count2 = 0;
      for (std::map<unsigned int, std::vector<TrajectoryMeasurement> >::iterator second = std::next(first,1); second != allHits.end() && allHits.size() > 1; ++second){
	++count2;
	

	if(printInfo) std::cout << " first = " << count1 << " second = " << count2 
				<< " => first->first = " << first->first<< " size map = " << allHits.size() << " trk Size = " << first->second.size() << std::endl;
	
	if(printInfo)    std::cout << "\n first = " << count1 << " => first->first = " << first->first<< " size map = " << allHits.size() 
				   << " trk Size = " << first->second.size() << std::endl;
	if(printInfo)      std::cout << " second = " << count2 << " => second->first = " << second->first << " trk Size = " << second->second.size() << std::endl;                                    
	
	
	//both tracks share with same tracks => already copied
	if(*(firstShareWithSecond[first->first].begin()) == *(firstShareWithSecond[second->first].begin())) {
	  if(printInfo)	std::cout << " " << first->first << " already shared with " << second->first << " through " << *(firstShareWithSecond[first->first].begin()) << std::endl;
	  continue;
	}

	float totalEshared = 0.;
	bool foundShared = false;
	bool foundSharedNhits = false;
	unsigned int totalShared = 0;

	bool foundSharedEnergy = false;
	unsigned int firstShared = 0;

	unsigned int secondTrackSize = second->second.size();
	unsigned int firstTrackSize = first->second.size();
	std::vector<bool> firstTrksharedHits( firstTrackSize, false);
	std::vector<bool> sharedHits( secondTrackSize, false);

	//invent something clever...
	//bool sharedShortTracks = false;
	//unsigned int minTracktotalEshared = 0;
	//unsigned int minTrackLength = std::min(secondTrackSize, firstTrackSize);
      
	unsigned int itTm1 = 0;
	for(auto tm1 : first->second){
	  
	  unsigned int itTm2 = 0;
	  for(auto tm2 : second->second){
	  
	  if (tm1.recHit()->sharesInput(&*tm2.recHit(), TrackingRecHit::all)) {
	    totalEshared += allHitsEnergyFr[second->first].at(itTm2);
	    ++totalShared;

	    //experimenting criteria...
	    /*
	    if(first->first == 0 || first->first == 71 || second->first == 71 || second->first == 200){
	      std::cout << "allHitsEnergyFr[first->first].size() = " << allHitsEnergyFr[first->first].size() << " firstTrackSize = " << firstTrackSize 
			<< " allHitsEnergyFr[second->first].size() = " << allHitsEnergyFr[second->first].size() << " secondTrackSize = " << secondTrackSize << std::endl;
	    }

	    minTracktotalEshared += (minTrackLength == firstTrackSize) ? allHitsEnergyFr[first->first].at(itTm1) : allHitsEnergyFr[second->first].at(itTm2);
	    */
	    /*
	    if(first->first == 0 || first->first == 71 || first->first == 72 ){
	      std::cout << " minTrackLength = " << minTrackLength << " totalShared = " << totalShared
			<< " totalEshared = " << totalEshared << " trsh minTrack = " << 1.2/(float)minTrackLength 
			<< " trsh 1st = " << 1.2/(float)firstTrackSize << " trsh 2nd = " << 1.2/(float)secondTrackSize 
			<< " minTracktotalEshared = " << minTracktotalEshared << std::endl;

	    }
	    */

	    
	    //flag as shared if common measurements account for > 90% of short track
	    //in other words: if shorter track (b) is common to a longer one (a) for more that its 90% of energy 
	    //consider b as part of a
	    //just flag the sharing once (foundSharedNhits)
	    if(totalEshared > 0.9 && foundSharedNhits == false){
	      firstShareWithSecond[second->first].insert(first->first);
	      foundSharedNhits = true;
	    }

	    if(printInfo) std::cout << ">>> track " << first->first << " shares with track " << second->first << " hit i " << itTm1 << " hit j " << itTm2 << std::endl;
	    foundShared = true;
	    //save info on shared measurements to avoid replica in the copy
	    sharedHits[itTm2] = true;
	    firstTrksharedHits[itTm1] = true;
	  }// found sharing
	  ++itTm2;
	}
	++itTm1;
      }//loop 2nd and 1st hits

      if(printInfo)       std::cout << " >>> foundShared = " << foundShared << " foundSharedNhits = " << foundSharedNhits << " totalEshared = " << totalEshared << std::endl;
      if(printInfo)      std::cout << " end loop trk i = " << first->first  << " and j = " <<  second->first << std::endl;


      //foundShared => the tracks share at least 1 meas
      //foundSharedNhits => the sharing criteria satisfy the condistion to merge the tracks
      //BUT if foundShared && !foundSharedNhits => end up with duplicated measurements for the moment!
      if(foundShared == true && foundSharedNhits == true){

	//all measurements already common => go directly to cleaning
	if(totalShared == secondTrackSize) {
	  std::cout << " secondTrackSize == totalShared "<< std::endl;
	  continue;
	}

	//identify the original track that "connects the sharing" and copy new measurements 
	//map order by key => begin with lower track_index == longer track
	unsigned int addToTrk = std::min(*(firstShareWithSecond[first->first].begin()), *(firstShareWithSecond[second->first].begin()));

	if(printInfo) {
	std::cout << " to be added to track i0 = " << *(firstShareWithSecond[first->first].begin()) << " j0 = " << *(firstShareWithSecond[second->first].begin())
		  << " min of 2 = " << addToTrk << std::endl;
	std::cout << " >>> found sharing, hits from trk " << second->first << " shared with track " << std::endl;
	for(auto pp : firstShareWithSecond[second->first]) std::cout << pp << std::endl;	
	}

	//min of the 2 is from trk 1 
	if(addToTrk == *(firstShareWithSecond[first->first].begin())){
	  unsigned int itTm2 = 0;
	  if(printInfo) std::cout << " 1st already shared -- add 2nd allHits[addToTrk].size() = " << allHits[addToTrk].size() << std::endl;
	  for(auto tm2 : second->second){
	    if(sharedHits[itTm2] == false){
	      allHits[addToTrk].push_back((tm2));
	      if(printInfo) 	    std::cout << " add hit " << itTm2 << std::endl;
	    }
	    ++itTm2;
	  }
	  if(printInfo) 	  std::cout << " post Add allHits[addToTrk].size() = " << allHits[addToTrk].size() << std::endl;
	}
	//min of the 2 is from trk 2 
	else if(addToTrk == *(firstShareWithSecond[second->first].begin())){
	  if(printInfo) std::cout << " >>> addToTrk !=  >>> add firstShareWithSecond: first = " << first->first << " second = " << addToTrk << std::endl; 
	  firstShareWithSecond[first->first].insert(addToTrk);

	  if(printInfo) 	  std::cout << " 2nd already shared -- add 1std allHits[addToTrk].size() = " << allHits[addToTrk].size() << std::endl;
	  // need also to copy eventual hits non shared
	  itTm1 = 0;
	  //for(auto tm1 : allHits[first->first]){
	  for(auto tm1 : first->second){
	    if(firstTrksharedHits[itTm1] == false) {
	      allHits[addToTrk].push_back((tm1));
	      if(printInfo) std::cout << " aggiungo also hit " << itTm1 << " from track " << first->first << std::endl;
	    }
	    ++itTm1;
	  }
	  if(printInfo) 	  std::cout << " post Add allHits[addToTrk].size() = " << allHits[addToTrk].size() << std::endl;
	}
      }//foundshared
      else if(foundShared == true){ // the tracks are sharing something but < 90% energy...
	// need to think something for this 
	// would like to get rid of duplication
	// merging everything that shares 1 measurement is not a solution...
	if(printInfo || printInfo3){
	  std::cout << " ############## shared hits but not enought => totalEshared = " << totalEshared << std::endl;
	  unsigned int itTm2 = 0;
	  for(auto tm2 : second->second){
	    //auto tm2 = *tm2P;
	    std::cout << " hit shared " << itTm2 << " = " << sharedHits[itTm2] << std::endl;
	    ++itTm2;
	  }
	  if(printInfo3)	std::cout << " ############## shared hits but not enought " << std::endl;
	}
	//need to think how to assign the shared hits
      }

      }//end loop over second
    }//loop over first


    if(printInfo)  std::cout << " start cleaning the map => remove shared tracks" << std::endl;

    for(auto ij: firstShareWithSecond){
      if(ij.second.size() == 1 && *(ij.second.begin()) == ij.first) {
	if(printInfo) std::cout << " non shared track " << ij.first << std::endl;
	continue;
      }
      else{
	if(*(ij.second.begin()) < ij.first) {
	  if(printInfo)	std::cout << " track " << ij.first;
	  for(auto kl: ij.second) {
	    if(printInfo)	  std::cout << " shared with track " << kl << std::endl;
	  }
	  
	  allHits[ij.first].clear(); 
	  allHits.erase(allHits.find(ij.first));
	  
	  if(printInfo)	std::cout << " >>> empty track = " << ij.first << std::endl;
	}
      }
    }
    
    //possible duplicated measurements exists, example tracks a,b,c 
    //a shares few measurements with b but not enough
    //a shared with c => copy c in a
    //b shared with c => copy b in a (c refers to a) 
    //if b_x is shared with a, but not with c => b_x will be replicated
    if(printInfo)   std::cout << " >>> now get rid of duplicates " << std::endl;


    allHitsEnergyFr.clear();
    for( auto& mapIt : allHits){
      if(printInfo)     std::cout << " trj = " << mapIt.first << std::endl;
      outCl.clear();
      outRh.clear();

      float trackSumEnergy = 0.;

      //if(mapIt.first == 0) continue;
      if(printInfo2) {
	std::cout << " mapIt.second.size() = " << mapIt.second.size() << std::endl;
      }
      unsigned int rhCount = 0;
      for(std::vector<TrajectoryMeasurement>::iterator rhIt = mapIt.second.begin(); rhIt != mapIt.second.end(); ++rhIt){
	if(printInfo)       std::cout << " rhCount = " << rhCount << " rhValid = " << (*rhIt->recHit()).isValid() <<std::endl;
	
	if (typeid(*rhIt->recHit()) == typeid(HGCTrackingRecHitFromCluster) ) {
	  if(printInfo) 	std::cout << " type Cl " << std::endl;
	  unsigned int prevClSize = outCl.size();
	  reco::CaloClusterPtr cl = (dynamic_cast<const HGCTrackingRecHitFromCluster&>(*rhIt->recHit())).objRef();
	  if(printInfo) 	std::cout << " cast ok insert " << std::endl;
	  outCl.insert(cl);
	  if(outCl.size() == prevClSize){
	    if(printInfo) 	  std::cout << " remove duplicated Cl " << std::endl;
	    mapIt.second.erase(rhIt);
	    --rhIt;
	  }
	  else{
	    float energy = cl->energy();
	    allHitsEnergyFr[mapIt.first].push_back(energy);
	    trackSumEnergy += energy;
	    if(printInfo) 	  std::cout << " save energy " << std::endl;
	  }
	} 
	else if (typeid(*rhIt->recHit()) == typeid(HGCTrackingClusteringRecHit)) {
	  //should implement case for this
	  //std::cout << " BAU BAU " << std::endl;
	} 
	else if (typeid(*rhIt->recHit()) == typeid(HGCTrackingRecHitFromHit)) {
	  // case HGCTrackingRecHitFromHit
	  if(printInfo) 	std::cout << " type Rh " << std::endl;
	  unsigned int prevRhSize = outRh.size();
	  const HGCRecHitRef rh = (dynamic_cast<const HGCTrackingRecHitFromHit&>(*rhIt->recHit()) ).objRef();
	  if(printInfo) 	std::cout << " cast ok insert " << std::endl;
	  outRh.insert(&*(rh));
	  if(outRh.size() == prevRhSize) {
	    if(printInfo) 	  std::cout << " remove duplicated Rh " << std::endl;
	    mapIt.second.erase(rhIt);
	    --rhIt;
	  }
	  else{
	    float energy = rh->energy();
	    allHitsEnergyFr[mapIt.first].push_back(energy);
	    trackSumEnergy += energy;
	    if(printInfo) 	  std::cout << " save energy " << std::endl;
	  }
	}
	else{
	  //	throw cms::Exception("Invalid valid hit", typeid(*meas[i].recHit()).name());
	  std::cout << " Invalid valid hit " << typeid(*rhIt->recHit()).name() << std::endl;
	  std::cout << " ??? " << std::endl;
	}
	++rhCount;
      }//measurements
      
    if(printInfo)           std::cout << " mapIt.second.size() = " << mapIt.second.size() << " and allHitsEnergyFr[mapIt.first].size() = " << allHitsEnergyFr[mapIt.first].size() << std::endl;
    
    if(trackSumEnergy > 0. && totLoop > 1) for(unsigned int iTM=0; iTM < allHitsEnergyFr[mapIt.first].size(); ++iTM) {
	allHitsEnergyFr[mapIt.first][iTM] /= trackSumEnergy;
	if(printInfo2 ) std::cout << " fraction = " << iTM << " mappa = " << allHitsEnergyFr[mapIt.first].at(iTM) << " size = " << allHitsEnergyFr[mapIt.first].size()  << std::endl;
      }
    }

    if(printInfo2) std::cout << " post allHits.size = " << allHits.size() << " allHitsEnergyFr.size() = " << allHitsEnergyFr.size() << std::endl;
    
    firstShareWithSecond.clear();
    if(totLoop > 1) for(auto ij : allHits){
	std::set<unsigned int> dummy;
	dummy.insert(ij.first);
	firstShareWithSecond[ij.first] = dummy;
      }
  
  }// loop over loops

  return allHits;
}
