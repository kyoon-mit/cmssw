#include <iostream>

// Author: Arabella Martelli, Felice Pantaleo, Marco Rovere
// arabella.martelli@cern.ch, felice.pantaleo@cern.ch, marco.rovere@cern.ch
// Date: 06/2019
#include <algorithm>
#include <set>
#include <vector>

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "SeedingRegionByTracks.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"
#include "DataFormats/GeometrySurface/interface/BoundDisk.h"

using namespace ticl;

SeedingRegionByTracks::SeedingRegionByTracks(const edm::ParameterSet &conf, edm::ConsumesCollector &sumes)
    : SeedingRegionAlgoBase(conf, sumes),
      tracks_token_(sumes.consumes<reco::TrackCollection>(conf.getParameter<edm::InputTag>("tracks"))),
      cutTk_(conf.getParameter<std::string>("cutTk")),
      detector_(conf.getParameter<std::string>("detector")),
      propName_(conf.getParameter<std::string>("propagator")),
      bfield_token_(sumes.esConsumes<MagneticField, IdealMagneticFieldRecord, edm::Transition::BeginRun>()),
      propagator_token_(sumes.esConsumes<Propagator, TrackingComponentsRecord, edm::Transition::BeginRun>(
          edm::ESInputTag("", propName_))) {
  std::string detectorName_ = (detector_ == "HFNose") ? "HGCalHFNoseSensitive" : "HGCalEESensitive";
  hdc_token_ = sumes.esConsumes<HGCalDDDConstants, IdealGeometryRecord, edm::Transition::BeginRun>(
      edm::ESInputTag("", detectorName_));
}

SeedingRegionByTracks::~SeedingRegionByTracks() {}

void SeedingRegionByTracks::initialize(const edm::EventSetup &es) {
  edm::ESHandle<HGCalDDDConstants> hdc = es.getHandle(hdc_token_);
  hgcons_ = hdc.product();

  buildFirstLayers();

  bfield_ = es.getHandle(bfield_token_);
  propagator_ = es.getHandle(propagator_token_);
}

void SeedingRegionByTracks::makeRegions(const edm::Event &ev,
                                        const edm::EventSetup &es,
                                        std::vector<TICLSeedingRegion> &result) {
  edm::Handle<reco::TrackCollection> tracks_h;
  ev.getByToken(tracks_token_, tracks_h);
  edm::ProductID trkId = tracks_h.id();
  auto bFieldProd = bfield_.product();
  const Propagator &prop = (*propagator_);

  int nTracks = tracks_h->size();

  for (int i = 0; i < nTracks; ++i) {
    const reco::Track &tk = (*tracks_h)[i];
    if (!cutTk_((tk))) {
      continue;
    }

    FreeTrajectoryState fts = trajectoryStateTransform::outerFreeState((tk), bFieldProd);
    int iSide = int(tk.eta() > 0);
    TrajectoryStateOnSurface tsos;
    
    if (detector_ != "HFNose") {
      tsos = prop.propagate(fts, firstDisk_[iSide]->surface());
    } else {
      std::cout << "TEPX --> HFNose propagator begins w/ " << propName_ << std::endl;

      if (propName_ == "AnalyticalPropagator") {
        // make sure to check theMaxRelativeChangeInBz variable in AnalyticalPropagator.cc
	// todo: implement steps
      }

      tsos = prop.propagate(fts, firstDisk_[iSide]->surface());

      std::cout << "tsos valid (bool): " << tsos.isValid() << std::endl;
      std::cout << "---- ends ----" << std::endl;
    }

    /*
    auto p0 = GlobalPoint(GlobalPoint::Polar(0.0603764,0,0));
    auto p1 = GlobalPoint(GlobalPoint::Polar(0.0603764,0,100));
    auto p2 = GlobalPoint(GlobalPoint::Polar(0.0603764,0,200));
    auto p3 = GlobalPoint(GlobalPoint::Polar(0.0603764,0,300));
    auto p4 = GlobalPoint(GlobalPoint::Polar(0.0603764,0,400));
    auto p5 = GlobalPoint(GlobalPoint::Polar(0.0603764,0,500));
    auto p6 = GlobalPoint(GlobalPoint::Polar(0.0603764,0,600));
    auto p7 = GlobalPoint(GlobalPoint::Polar(0.0603764,0,700));
    auto p8 = GlobalPoint(GlobalPoint::Polar(0.0603764,0,800));
    auto p9 = GlobalPoint(GlobalPoint::Polar(0.0603764,0,900));
    auto p10 = GlobalPoint(GlobalPoint::Polar(0.0603764,0,1000));
    auto p11 = GlobalPoint(GlobalPoint::Polar(0.0603764,0,1100));
    auto p12 = GlobalPoint(GlobalPoint::Polar(0.0603764,0,1200));
    auto p13 = GlobalPoint(GlobalPoint::Polar(0.0603764,0,-100));
    auto p14 = GlobalPoint(GlobalPoint::Polar(0.0603764,0,-200));
    auto p15 = GlobalPoint(GlobalPoint::Polar(0.0603764,0,-300));
    auto p16 = GlobalPoint(GlobalPoint::Polar(0.0603764,0,-400));
    auto p17 = GlobalPoint(GlobalPoint::Polar(0.0603764,0,-500));
    auto p18 = GlobalPoint(GlobalPoint::Polar(0.0603764,0,-600));
    auto p19 = GlobalPoint(GlobalPoint::Polar(0.0603764,0,-700));
    auto p20 = GlobalPoint(GlobalPoint::Polar(0.0603764,0,-800));
    auto p21 = GlobalPoint(GlobalPoint::Polar(0.0603764,0,-900));
    auto p22 = GlobalPoint(GlobalPoint::Polar(0.0603764,0,-1000));
    auto p23 = GlobalPoint(GlobalPoint::Polar(0.0603764,0,-1100));
    auto p24 = GlobalPoint(GlobalPoint::Polar(0.0603764,0,-1200));

    std::cout << "(theta,phi) = (Polar(0.0603764,0)" <<  std::endl;
    std::cout << "check eta, z: " << p10.eta() << " " << p10.z() << std::endl;
    std::cout << "p0: " << p0 << std::endl;

    std::cout << "r = 0: " << bFieldProd->inTesla(p0).z() << std::endl;
    std::cout << "r = 100: " << bFieldProd->inTesla(p1).z() << std::endl;
    std::cout << "r = 200: " << bFieldProd->inTesla(p2).z() << std::endl;
    std::cout << "r = 300: " << bFieldProd->inTesla(p3).z() << std::endl;
    std::cout << "r = 400: " << bFieldProd->inTesla(p4).z() << std::endl;
    std::cout << "r = 500: " << bFieldProd->inTesla(p5).z() << std::endl;
    std::cout << "r = 600: " << bFieldProd->inTesla(p6).z() << std::endl;
    std::cout << "r = 700: " << bFieldProd->inTesla(p7).z() << std::endl;
    std::cout << "r = 800: " << bFieldProd->inTesla(p8).z() << std::endl;
    std::cout << "r = 900: " << bFieldProd->inTesla(p9).z() << std::endl;
    std::cout << "r = 1000: " << bFieldProd->inTesla(p10).z() << std::endl;
    std::cout << "r = 1100: " << bFieldProd->inTesla(p11).z() << std::endl;
    std::cout << "r = 1200: " << bFieldProd->inTesla(p12).z() << std::endl;
    std::cout << "r = -100: " << bFieldProd->inTesla(p13).z() << std::endl;
    std::cout << "r = -200: " << bFieldProd->inTesla(p14).z() << std::endl;
    std::cout << "r = -300: " << bFieldProd->inTesla(p15).z() << std::endl;
    std::cout << "r = -400: " << bFieldProd->inTesla(p16).z() << std::endl;
    std::cout << "r = -500: " << bFieldProd->inTesla(p17).z() << std::endl;
    std::cout << "r = -600: " << bFieldProd->inTesla(p18).z() << std::endl;
    std::cout << "r = -700: " << bFieldProd->inTesla(p19).z() << std::endl;
    std::cout << "r = -800: " << bFieldProd->inTesla(p20).z() << std::endl;
    std::cout << "r = -900: " << bFieldProd->inTesla(p21).z() << std::endl;
    std::cout << "r = -1000: " << bFieldProd->inTesla(p22).z() << std::endl;
    std::cout << "r = -1100: " << bFieldProd->inTesla(p23).z() << std::endl;
    std::cout << "r = -1200: " << bFieldProd->inTesla(p24).z() << std::endl;
    */

    if (tsos.isValid()) {
      result.emplace_back(tsos.globalPosition(), tsos.globalMomentum(), iSide, i, trkId);
    }
  }

  // sorting seeding region by descending momentum
  std::sort(result.begin(), result.end(), [](const TICLSeedingRegion &a, const TICLSeedingRegion &b) {
    return a.directionAtOrigin.perp2() > b.directionAtOrigin.perp2();
  });
}

void SeedingRegionByTracks::fillPSetDescription(edm::ParameterSetDescription &desc) {
  desc.add<edm::InputTag>("tracks", edm::InputTag("generalTracks"));
  desc.add<std::string>("cutTk",
                        "1.48 < abs(eta) < 3.0 && pt > 1. && quality(\"highPurity\") && "
                        "hitPattern().numberOfLostHits(\"MISSING_OUTER_HITS\") < 5");
  desc.add<std::string>("propagator", "PropagatorWithMaterial");
  desc.add<std::string>("detector", "HGCAL");
  SeedingRegionAlgoBase::fillPSetDescription(desc);
}

void SeedingRegionByTracks::buildFirstLayers() {
  float zVal = hgcons_->waferZ(1, true);
  std::pair<double, double> rMinMax = hgcons_->rangeR(zVal, true);

  for (int iSide = 0; iSide < 2; ++iSide) {
    float zSide = (iSide == 0) ? (-1. * zVal) : zVal;
    firstDisk_[iSide] =
        std::make_unique<GeomDet>(Disk::build(Disk::PositionType(0, 0, zSide),
                                              Disk::RotationType(),
                                              SimpleDiskBounds(rMinMax.first, rMinMax.second, zSide - 0.5, zSide + 0.5))
                                      .get());
  }
}
