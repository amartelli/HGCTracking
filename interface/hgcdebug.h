#ifndef RecoHGCal_HGCTracking_hgcdebug_h
#define RecoHGCal_HGCTracking_hgcdebug_h

#include "SimDataFormats/CaloAnalysis/interface/CaloParticle.h"

namespace hgctracking { 
    extern int g_debuglevel;
}

typedef std::unordered_multimap<uint32_t,std::pair<const CaloParticle *,float>> CaloTruthRevMap;

#endif
