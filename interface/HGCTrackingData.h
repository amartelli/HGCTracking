#ifndef RecoHGCal_HGCTracking_interface_HGCTrackingData_h
#define RecoHGCal_HGCTracking_interface_HGCTrackingData_h

#include "RecoHGCal/HGCTracking/interface/HGCTrackingDiskData.h"
#include "RecoHGCal/HGCTracking/interface/HGCTracker.h"
#include <map>

/// Class that holds the rechits of HGC Layer, and does the lookup
//  Almost like a MeasurementTrackerEvent in the tracker
//
//  note: the class may be expensive both to copy and to move


class HGCTrackingData {
    public:
        typedef HGCTrackingDiskData Disk;
        typedef Disk::TColl TColl;

        HGCTrackingData(const HGCTracker &tracker, const HGCTrackingBasicCPE *cpe) : tracker_(&tracker), cpe_(cpe) {}
        void addData(const edm::Handle<TColl> &data, int subdet, float thrSoN=1.) ;
        void addClusters(const edm::Handle<reco::CaloClusterCollection> &data) ;
        const Disk & diskData(const HGCDiskGeomDet *disk) const { 
            auto it = data_.find(disk);
            if (it == data_.end()) throwBadDisk_(disk);
            return it->second; 
        }

    private:
        const HGCTracker * tracker_;
        const HGCTrackingBasicCPE* cpe_;
        std::map<const HGCDiskGeomDet *, Disk> data_; // FIXME better type later

        void throwBadDisk_(const HGCDiskGeomDet *disk) const ;
};


#endif
