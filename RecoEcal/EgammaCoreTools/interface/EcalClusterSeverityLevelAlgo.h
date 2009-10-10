#ifndef RecoEcal_EgammaCoreTools_EcalClusterSeverityLevelAlgo_hh
#define RecoEcal_EgammaCoreTools_EcalClusterSeverityLevelAlgo_hh

#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/CaloRecHit/interface/CaloCluster.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "RecoLocalCalo/EcalRecAlgos/interface/EcalSeverityLevelAlgo.h"
#include "CondFormats/EcalObjects/interface/EcalChannelStatus.h"

class CaloTopology;
class EBDetId;
class EcalClusterSeverityLevelAlgo {
 public:

  // the severity is the fraction of cluster energy 
  // taken by good channels  
  // (e.g. not noisy, dead and recovered etc.)
  static float goodFraction( const reco::CaloCluster & , const EcalRecHitCollection &,  const EcalChannelStatus & );
  // fraction of SC energy around closest problematic
  static float fractionAroundClosestProblematic( const reco::CaloCluster & , const EcalRecHitCollection &,  const EcalChannelStatus &, const CaloTopology* topology );
  // retrieve closest problematic channel wrt seed crystal using as distance sqrt(ieta^2+ieta^2+iphi^2+iphi^2). Return a null detId in case not found within a search region of 11 (ieta) x 51 (iphi)  
  static DetId closestProblematic( const reco::CaloCluster & , const EcalRecHitCollection &,  const EcalChannelStatus &, const CaloTopology* topology );
  // retrieve the distance in ieta,iphi (number of crystals) of the closest problematic channel wrt seed crystal (defined as above)
  // return -1,-1 if no crystal is found within a search region of 11 (eta) x 51 (phi)
  static std::pair<int,int> etaphiDistanceClosestProblematic( const reco::CaloCluster & , const EcalRecHitCollection &,  const EcalChannelStatus &, const CaloTopology* topology );

 private:

  //return the distance in eta units between two EBDetId
  static int distanceEta(const EBDetId& a,const EBDetId& b); 
  //return the distance in phi units between two EBDetId
  static int distancePhi(const EBDetId& a,const EBDetId& b); 
  //return the distance in eta units between two EBDetId
  static int distanceX(const EEDetId& a,const EEDetId& b); 
  //return the distance in phi units between two EBDetId
  static int distanceY(const EEDetId& a,const EEDetId& b); 
  //return the chStatus as taken from the DB object
  static uint16_t retrieveDBStatus( const DetId , const EcalChannelStatus &chStatus );

};

#endif
