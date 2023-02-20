// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

///
/// \file   qaKFParticleLc.cxx
/// \author Annalena Kalteyer <annalena.sophie.kalteyer@cern.ch>, GSI Darmstadt
/// \brief  Task to test the performance of the KFParticle package on the Lc to pKpi decay
///

#include "Tools/KFparticle/qaKFParticle.h"
#include <CCDB/BasicCCDBManager.h>
#include <string>
#include <TDatabasePDG.h>
#include <TPDGCode.h>
#include "TableHelper.h"

/// includes O2
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/DCA.h"
#include "ReconstructionDataFormats/Track.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "DetectorsBase/GeometryManager.h"
#include "DetectorsBase/Propagator.h"

/// includes O2Physics
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/Core/trackUtilities.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/Core/RecoDecay.h"
#include "Common/Core/TrackSelectorPID.h"
#include "Tools/KFparticle/KFUtilities.h"

/// includes KFParticle
#include "KFParticle.h"
#include "KFPTrack.h"
#include "KFPVertex.h"
#include "KFParticleBase.h"
#include "KFVertex.h"

#ifndef HomogeneousField

#define HomogeneousField

#endif

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::dataformats;

struct qaKFParticle {

  /// general steering settings
  Configurable<bool> isRun3{"isRun3", true, "Is Run3 dataset"};
  Configurable<std::string> ccdbUrl{"ccdburl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> ccdbPathLut{"ccdbPathLut", "GLO/Param/MatLUT", "Path for LUT parametrization"};
  Configurable<std::string> ccdbPathGeo{"ccdbPathGeo", "GLO/Config/GeometryAligned", "Path of the geometry file"};
  Configurable<std::string> ccdbPathGrp{"ccdbPathGrp", "GLO/GRP/GRP", "Path of the grp file (Run 2)"};
  Configurable<std::string> ccdbPathGrpMag{"ccdbPathGrpMag", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object (Run 3)"};
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  o2::base::MatLayerCylSet* lut;
  o2::base::Propagator::MatCorrType matCorr = o2::base::Propagator::MatCorrType::USEMatCorrLUT;
  int runNumber;
  double magneticField = 0.;

  /// Histogram Configurables
  ConfigurableAxis binsPt{"binsPt", {VARIABLE_WIDTH, 0.0, 1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 11., 12., 13., 14., 15., 16., 17., 24., 36., 50.0}, ""};

  /// option to select good events
  Configurable<bool> eventSelection{"eventSelection", true, "select good events"}; // currently only sel8 is defined for run3
  /// options to select only specific tracks
  Configurable<int> trackSelection{"trackSelection", 1, "Track selection: 0 -> No Cut, 1 -> kGlobalTrack, 2 -> kGlobalTrackWoPtEta, 3 -> kGlobalTrackWoDCA, 4 -> kQualityTracks, 5 -> kInAcceptanceTracks"};

  /// Particle Identification
  // TPC PID
  Configurable<double> ptPidTpcMinPi{"ptPidTpcMinPi", 0.15, "Lower bound of track pT for TPC PID"};
  Configurable<double> ptPidTpcMaxPi{"ptPidTpcMaxPi", 5., "Upper bound of track pT for TPC PID"};
  Configurable<double> nSigmaTpcMaxPi{"nSigmaTpcMaxPi", 3., "Nsigma cut on TPC only"};
  Configurable<double> nSigmaTpcCombinedMaxPi{"nSigmaTpcCombinedMaxPi", 5., "Nsigma cut on TPC combined with TOF"};
  Configurable<double> ptPidTpcMinKa{"ptPidTpcMinKa", 0.15, "Lower bound of track pT for TPC PID"};
  Configurable<double> ptPidTpcMaxKa{"ptPidTpcMaxKa", 5., "Upper bound of track pT for TPC PID"};
  Configurable<double> nSigmaTpcMaxKa{"nSigmaTpcMaxKa", 3., "Nsigma cut on TPC only"};
  Configurable<double> nSigmaTpcCombinedMaxKa{"nSigmaTpcCombinedMaxKa", 5., "Nsigma cut on TPC combined with TOF"};
  Configurable<double> ptPidTpcMinPr{"ptPidTpcMinPr", 0.15, "Lower bound of track pT for TPC PID"};
  Configurable<double> ptPidTpcMaxPr{"ptPidTpcMaxPr", 5., "Upper bound of track pT for TPC PID"};
  Configurable<double> nSigmaTpcMaxPr{"nSigmaTpcMaxPr", 3., "Nsigma cut on TPC only"};
  Configurable<double> nSigmaTpcCombinedMaxPr{"nSigmaTpcCombinedMaxPr", 5., "Nsigma cut on TPC combined with TOF"};
  // TOF PID
  Configurable<double> ptPidTofMinPi{"ptPidTofMinPi", 0.15, "Lower bound of track pT for TOF PID"};
  Configurable<double> ptPidTofMaxPi{"ptPidTofMaxPi", 5., "Upper bound of track pT for TOF PID"};
  Configurable<double> nSigmaTofMaxPi{"nSigmaTofMaxPi", 3., "Nsigma cut on TOF only"};
  Configurable<double> nSigmaTofCombinedMaxPi{"nSigmaTofCombinedMaxPi", 5., "Nsigma cut on TOF combined with TPC"};
  Configurable<double> ptPidTofMinKa{"ptPidTofMinKa", 0.15, "Lower bound of track pT for TOF PID"};
  Configurable<double> ptPidTofMaxKa{"ptPidTofMaxKa", 5., "Upper bound of track pT for TOF PID"};
  Configurable<double> nSigmaTofMaxKa{"nSigmaTofMaxKa", 3., "Nsigma cut on TOF only"};
  Configurable<double> nSigmaTofCombinedMaxKa{"nSigmaTofCombinedMaxKa", 5., "Nsigma cut on TOF combined with TPC"};
  Configurable<double> ptPidTofMinPr{"ptPidTofMinPr", 0.15, "Lower bound of track pT for TOF PID"};
  Configurable<double> ptPidTofMaxPr{"ptPidTofMaxPr", 5., "Upper bound of track pT for TOF PID"};
  Configurable<double> nSigmaTofMaxPr{"nSigmaTofMaxPr", 3., "Nsigma cut on TOF only"};
  Configurable<double> nSigmaTofCombinedMaxPr{"nSigmaTofCombinedMaxPr", 5., "Nsigma cut on TOF combined with TPC"};
  /// singe track selections
  Configurable<float> d_pTMin{"d_pTMin", 0.3, "minimum momentum for tracks"};
  Configurable<float> d_etaRange{"d_etaRange", 0.8, "eta Range for tracks"};
  Configurable<float> d_dcaXYTrackPV{"d_dcaXYTrackPV", 2., "DCA XY of the daughter tracks to the PV"};
  Configurable<float> d_dcaZTrackPV{"d_dcaZTrackPV", 10., "DCA Z of the daughter tracks to the PV"};
  /// D0 selections
  Configurable<bool> applySelectionDoWithTopoConst{"applySelectionDoWithTopoConst", true, "Apply selections on the D0 after constraining it to the PV"};
  Configurable<float> d_pTMinD0{"d_pTMinD0", 0., "minimum momentum for D0 candidates"};
  Configurable<float> d_pTMaxD0{"d_pTMaxD0", 36., "maximum momentum for D0 candidates"};
  Configurable<float> d_massMinD0{"d_massMinD0", 1.65, "minimum mass for D0"};
  Configurable<float> d_massMaxD0{"d_massMaxD0", 2.08, "minimum mass for D0"};
  Configurable<float> d_cosPA{"d_cosPA", -1., "minimum cosine Pointing angle for D0"};
  Configurable<float> d_cosPAXY{"d_cosPAXY", -1., "minimum cosine Pointing angle for D0"};
  Configurable<float> d_decayLength{"d_decayLength", 0., "minimum decay length for D0"};
  Configurable<float> d_normdecayLength{"d_normdecayLength", 100., "minimum normalised decay length for D0"};
  Configurable<float> d_chi2topoD0{"d_chi2topoD0", 1000., "maximum chi2 topological of D0 to PV"};
  Configurable<float> d_dist3DSVDau{"d_dist3DSVDau", 1000., "maximum geometrical distance 3D daughter tracks at the SV"};
  Configurable<float> d_cosThetaStar{"d_cosThetaStarPi", 1000., "maximum cosine theta star"};
  Configurable<float> d_distPiToSV{"d_distPiToSV", 1000., "maximum distance Pi to SV"};
  Configurable<float> d_distKaToSV{"d_distKaToSV", 1000., "maximum distance Ka to SV"};
  Configurable<float> d_d0Pid0Ka{"d_d0Pid0Ka", 1000., "Product impact parameters"};
  /// Option to write D0 variables in a tree
  Configurable<double> d_DwnSmplFact{"d_DwnSmplFact", 1., "Downsampling factor for tree"};
  Configurable<bool> writeTree{"writeTree", false, "write daughter variables in a tree"};
  Configurable<bool> writeQAHistograms{"writeQAHistograms", false, "write all QA histograms"};

  // Define which track selection should be used:
  // 0 -> No track selection is applied
  // 1 kGlobalTrack = kQualityTracks | kPrimaryTracks | kInAcceptanceTracks
  //        kQualityTracks = kTrackType | kTPCNCls | kTPCCrossedRows | kTPCCrossedRowsOverNCls | kTPCChi2NDF |
  //                         kTPCRefit | kITSNCls | kITSChi2NDF | kITSRefit | kITSHits
  //        kPrimaryTracks = kGoldenChi2 | kDCAxy | kDCAz
  //        kInAcceptanceTracks = kPtRange | kEtaRange
  // 2 kGlobalTrackWoPtEta = kQualityTracks | kPrimaryTracks
  // 3 kGlobalTrackWoDCA = kQualityTracks | kInAcceptanceTracks
  // 4 kQualityTracks
  // 5 kInAcceptanceTracks
  Filter trackFilter = (trackSelection.node() == 0) ||
                       ((trackSelection.node() == 1) && requireGlobalTrackInFilter()) ||
                       ((trackSelection.node() == 2) && requireGlobalTrackWoPtEtaInFilter()) ||
                       ((trackSelection.node() == 3) && requireGlobalTrackWoDCAInFilter()) ||
                       ((trackSelection.node() == 4) && requireQualityTracksInFilter()) ||
                       ((trackSelection.node() == 5) && requireTrackCutInFilter(TrackSelectionFlags::kInAcceptanceTracks));

  using CollisionTableData = soa::Join<aod::Collisions, aod::EvSels, aod::Mults>;
  using BigTracks = soa::Join<aod::Tracks, aod::TracksCov, aod::TracksExtra>;
  using BigTracksExtended = soa::Join<BigTracks, aod::TracksDCA>;
  using BigTracksPID = soa::Join<BigTracksExtended, aod::pidTPCFullEl, aod::pidTPCFullMu, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr, aod::pidTOFFullEl, aod::pidTOFFullMu, aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr>;
  using TrackTableData = soa::Join<BigTracksPID, aod::TrackSelection>;

  /// Table to be produced
  Produces<o2::aod::TreeKF> rowKF;

  void initMagneticFieldCCDB(o2::aod::BCsWithTimestamps::iterator const& bc, int& mRunNumber,
                             o2::framework::Service<o2::ccdb::BasicCCDBManager> const& ccdb, std::string ccdbPathGrp, o2::base::MatLayerCylSet* lut,
                             bool isRun3)
  {

    if (mRunNumber != bc.runNumber()) {

      LOGF(info, "====== initCCDB function called (isRun3==%d)", isRun3);
      if (!isRun3) { // Run 2 GRP object
        o2::parameters::GRPObject* grpo = ccdb->getForTimeStamp<o2::parameters::GRPObject>(ccdbPathGrp, bc.timestamp());
        if (grpo == nullptr) {
          LOGF(fatal, "Run 2 GRP object (type o2::parameters::GRPObject) is not available in CCDB for run=%d at timestamp=%llu", bc.runNumber(), bc.timestamp());
        }
        o2::base::Propagator::initFieldFromGRP(grpo);
        o2::base::Propagator::Instance()->setMatLUT(lut);
        LOGF(info, "Setting magnetic field to %d kG for run %d from its GRP CCDB object (type o2::parameters::GRPObject)", grpo->getNominalL3Field(), bc.runNumber());
      } else { // Run 3 GRP object
        o2::parameters::GRPMagField* grpo = ccdb->getForTimeStamp<o2::parameters::GRPMagField>(ccdbPathGrp, bc.timestamp());
        if (grpo == nullptr) {
          LOGF(fatal, "Run 3 GRP object (type o2::parameters::GRPMagField) is not available in CCDB for run=%d at timestamp=%llu", bc.runNumber(), bc.timestamp());
        }
        o2::base::Propagator::initFieldFromGRP(grpo);
        o2::base::Propagator::Instance()->setMatLUT(lut);
        LOGF(info, "Setting magnetic field to current %f A for run %d from its GRP CCDB object (type o2::parameters::GRPMagField)", grpo->getL3Current(), bc.runNumber());
      }
      mRunNumber = bc.runNumber();
    }
  } /// end initMagneticFieldCCDB

  void init(InitContext const&)
  {
    ccdb->setURL(ccdbUrl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    lut = o2::base::MatLayerCylSet::rectifyPtrFromFile(ccdb->get<o2::base::MatLayerCylSet>(ccdbPathLut));
    if (!o2::base::GeometryManager::isGeometryLoaded()) {
      ccdb->get<TGeoManager>(ccdbPathGeo);
    }
    runNumber = 0;
  } /// End init

  /// Function to select collisions
  template <typename T>
  bool isSelectedCollision(const T& collision)
  {
    /// Trigger selection
    if (eventSelection && !(isRun3 ? collision.sel8() : collision.sel7())) { // currently only sel8 is defined for run3
      return false;
    }
    /// Reject collisions with negative covariance matrix elemts on the digonal
    if (collision.covXX() < 0. || collision.covYY() < 0. || collision.covZZ() < 0.) {
      return false;
    }
    return true;
  }

  /// Function for single track selection
  template <typename T>
  bool isSelectedTracks(const T& track1, const T& track2, const T& track3)
  {
    if ((track1.p() < d_pTMin) ||  (track2.p() < d_pTMin) || (track3.p() < d_pTMin)) {
      return false;
    }
    /// Eta range
    if ((abs(track1.eta()) > d_etaRange) || (abs(track2.eta()) > d_etaRange) || (abs(track3.eta()) > d_etaRange)) {
      return false;
    }
    /// DCA XY of the daughter tracks to the primaty vertex
    if ((track1.dcaXY() > d_dcaXYTrackPV) || (track2.dcaXY() > d_dcaXYTrackPV) || (track3.dcaXY() > d_dcaXYTrackPV)) {
      return false;
    }
    /// DCA Z of the daughter tracks to the primaty vertex
    if ((track1.dcaZ() > d_dcaZTrackPV) || (track2.dcaZ() > d_dcaZTrackPV) || (track3.dcaZ() > d_dcaZTrackPV)) {
      return false;
    }
    /// reject if the tracks with sum of charge != 1
    if (abs(track1.sign() + track2.sign() + track3.sign()) != 1) {
      return false;
    }
    return true;
  }

  template <typename T>
  bool isSelectedDaughters(const T& KFPion, const T& KFKaon, const T& KFDZero, const T& KFPV)
  {
    /// distance 3D daughter tracks at the secondary vertex
    if (KFPion.GetDistanceFromParticle(KFKaon) > d_dist3DSVDau) {
      return false;
    }
    /// distance Pi to SV
    if (KFPion.GetDistanceFromVertex(KFDZero) > d_distPiToSV) {
      return false;
    }
    /// distance Ka to SV
    if (KFKaon.GetDistanceFromVertex(KFDZero) > d_distKaToSV) {
      return false;
    }
    float d0pid0ka = KFPion.GetDistanceFromVertexXY(KFPV) * KFKaon.GetDistanceFromVertexXY(KFPV);
    if (d0pid0ka > d_d0Pid0Ka) {
      return false;
    }
    return true;
  }

  template <typename T>
  bool isSelectedDoGeo(const T& KFDZero, const T& KFPV, float cosThetaStar)
  {
    /// Pt selection
    if (KFDZero.GetPt() < d_pTMinD0 || KFDZero.GetPt() > d_pTMaxD0) {
      return false;
    }
    /// Mass window selection
    if (KFDZero.GetMass() < d_massMinD0 || KFDZero.GetMass() > d_massMaxD0) {
      return false;
    }
    /// cosine pointing angle selection
    if (cpaFromKF(KFDZero, KFPV) < d_cosPA) {
      return false;
    }
    /// cosine pointing XY angle selection
    if (cpaXYFromKF(KFDZero, KFPV) < d_cosPAXY) {
      return false;
    }
    /// cosine theta star
    if (cosThetaStar > d_cosThetaStar) {
      return false;
    }
    return true;
  }

  template <typename T>
  bool isSelectedDoTopo(const T& KFDZero_PV, const T& KFPion, const T& KFKaon, const T& KFDZero_DecayVtx, const T& KFPV)
  {
    /// Pt selection
    if (KFDZero_PV.GetPt() < d_pTMinD0 || KFDZero_PV.GetPt() > d_pTMaxD0) {
      return false;
    }
    /// Mass window selection
    if (KFDZero_PV.GetMass() == 0.) {
      return false;
    }
    /// cosine pointing angle selection
    if (cpaFromKF(KFDZero_DecayVtx, KFPV) < d_cosPA) {
      return false;
    }
    /// decay length selection
    if (KFDZero_PV.GetDecayLength() < d_decayLength) {
      return false;
    }
    /// decay length error selection
    float normdecayLength = KFDZero_PV.GetDecayLength() / KFDZero_PV.GetErrDecayLength();
    if (normdecayLength < d_normdecayLength) {
      return false;
    }
    /// chi2 topological of DZero to PV
    float chi2topo = KFDZero_PV.GetChi2() / KFDZero_PV.GetNDF();
    if (chi2topo > d_chi2topoD0) {
      return false;
    }
    return true;
  }

  template <typename T1, typename T2, typename T3>
  void writeVarTree(const T1& kfpTrackPi, const T1& kfpTrackKa, const T2& KFPion, const T2& KFKaon, const T2& KFDZero_PV, const T2& KFDZero, const T2& KFPV, const T2& KFDZero_DecayVtx, float TPCnSigmaPi, float TOFnSigmaPi, float TPCnSigmaKa, float TOFnSigmaKa, float cosThetaStar, const T3& track1, const int source)
  {

    float d0pid0ka = KFPion.GetDistanceFromVertexXY(KFPV) * KFKaon.GetDistanceFromVertexXY(KFPV);
    float chi2geo = KFDZero.GetChi2() / KFDZero.GetNDF();
    float normdecayLength = KFDZero_PV.GetDecayLength() / KFDZero_PV.GetErrDecayLength();
    float chi2topo = KFDZero_PV.GetChi2() / KFDZero_PV.GetNDF();
    const double pseudoRndm = track1.pt() * 1000. - (int64_t)(track1.pt() * 1000);
    if (pseudoRndm < d_DwnSmplFact) {
      if (writeTree) {
        /// Filling the D0 tree
        rowKF(runNumber,
              KFPion.GetPt(),
              KFKaon.GetPt(),
              KFPion.GetDistanceFromVertexXY(KFPV),
              KFKaon.GetDistanceFromVertexXY(KFPV),
              KFPion.GetDistanceFromVertex(KFPV),
              KFKaon.GetDistanceFromVertex(KFPV),
              TPCnSigmaPi,
              TPCnSigmaKa,
              TOFnSigmaPi,
              TOFnSigmaKa,
              KFPion.GetDistanceFromVertexXY(KFDZero),
              KFKaon.GetDistanceFromVertexXY(KFDZero),
              cosThetaStar,
              KFPion.GetDistanceFromParticle(KFKaon),
              d0pid0ka,
              KFDZero.GetPt(),
              KFDZero.GetMass(),
              cpaFromKF(KFDZero, KFPV),
              cpaXYFromKF(KFDZero, KFPV),
              KFDZero.GetDistanceFromVertex(KFPV),
              KFDZero.GetDistanceFromVertexXY(KFPV),
              chi2geo,
              KFDZero_PV.GetPt(),
              KFDZero_PV.GetMass(),
              KFDZero_PV.GetDecayLength(),
              KFDZero_PV.GetDecayLengthXY(),
              cpaFromKF(KFDZero_DecayVtx, KFPV),
              KFDZero_PV.GetLifeTime(),
              normdecayLength,
              KFDZero_PV.GetDistanceFromVertex(KFPV),
              KFDZero_PV.GetDistanceFromVertexXY(KFPV),
              chi2topo,
              source);
      }
    }
  }

  /// Process function for data
  void processData(CollisionTableData::iterator const& collision, soa::Filtered<TrackTableData> const& tracks, aod::BCsWithTimestamps const&)
  {

    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    if (runNumber != bc.runNumber()) {
      initMagneticFieldCCDB(bc, runNumber, ccdb, isRun3 ? ccdbPathGrpMag : ccdbPathGrp, lut, isRun3);
      magneticField = o2::base::Propagator::Instance()->getNominalBz();
      /// Set magnetic field for KF vertexing
      #ifdef HomogeneousField
      KFParticle::SetField(magneticField);
      #endif
    }
    /// Apply event selection
    if (!isSelectedCollision(collision)) {
      return;
    }
    /// set KF primary vertex
    KFPVertex kfpVertex = createKFPVertexFromCollision(collision);
    KFParticle KFPV(kfpVertex);

    TrackSelectorPID selectorPion(kPiPlus);
    selectorPion.setRangePtTPC(ptPidTpcMinPi, ptPidTpcMaxPi);
    selectorPion.setRangeNSigmaTPC(-nSigmaTpcMaxPi, nSigmaTpcMaxPi);
    selectorPion.setRangeNSigmaTPCCondTOF(-nSigmaTpcCombinedMaxPi, nSigmaTpcCombinedMaxPi);
    selectorPion.setRangePtTOF(ptPidTofMinPi, ptPidTofMaxPi);
    selectorPion.setRangeNSigmaTOF(-nSigmaTofMaxPi, nSigmaTofMaxPi);
    selectorPion.setRangeNSigmaTOFCondTPC(-nSigmaTofCombinedMaxPi, nSigmaTofCombinedMaxPi);

    TrackSelectorPID selectorKaon(kKPlus);
    selectorKaon.setRangePtTPC(ptPidTpcMinKa, ptPidTpcMaxKa);
    selectorKaon.setRangeNSigmaTPC(-nSigmaTpcMaxKa, nSigmaTpcMaxKa);
    selectorKaon.setRangeNSigmaTPCCondTOF(-nSigmaTpcCombinedMaxKa, nSigmaTpcCombinedMaxKa);
    selectorKaon.setRangePtTOF(ptPidTofMinKa, ptPidTofMaxKa);
    selectorKaon.setRangeNSigmaTOF(-nSigmaTofMaxKa, nSigmaTofMaxKa);
    selectorKaon.setRangeNSigmaTOFCondTPC(-nSigmaTofCombinedMaxKa, nSigmaTofCombinedMaxKa);

    TrackSelectorPID selectorProton(kProton);
    selectorProton.setRangePtTPC(ptPidTpcMinPr, ptPidTpcMaxPr);
    selectorProton.setRangeNSigmaTPC(-nSigmaTpcMaxPr, nSigmaTpcMaxPr);
    selectorProton.setRangeNSigmaTPCCondTOF(-nSigmaTpcCombinedMaxPr, nSigmaTpcCombinedMaxPr);
    selectorProton.setRangePtTOF(ptPidTofMinPr, ptPidTofMaxPr);
    selectorProton.setRangeNSigmaTOF(-nSigmaTofMaxPr, nSigmaTofMaxPr);
    selectorProton.setRangeNSigmaTOFCondTPC(-nSigmaTofCombinedMaxPr, nSigmaTofCombinedMaxPr);
    

    for (auto& [track1, track2, track3] : combinations(soa::CombinationsStrictlyUpperIndexPolicy(tracks, tracks, tracks))) {

      /// Apply single track selection
      if (!isSelectedTracks(track1, track2, track3)) {
        continue;
      }

      KFPTrack kfpTrackPi;
      KFPTrack kfpTrackKa;
      KFPTrack kfpTrackPr;

      bool CandLc = false;
      bool CandLcbar = false;

      float TPCnSigmaPi = 0;
      float TPCnSigmaKa = 0;
      float TPCnSigmaPr = 0;
      float TOFnSigmaPi = 0;
      float TOFnSigmaKa = 0;
      float TOFnSigmaPr = 0;
      int source = 0;

      int pidKaonTr1 = selectorKaon.getStatusTrackPIDAll(track1);
      int pidKaonTr2 = selectorKaon.getStatusTrackPIDAll(track2);
      int pidKaonTr3 = selectorKaon.getStatusTrackPIDAll(track3);
      int pidPionTr1 = selectorPion.getStatusTrackPIDAll(track1);
      int pidPionTr2 = selectorPion.getStatusTrackPIDAll(track2);
      int pidPionTr3 = selectorPion.getStatusTrackPIDAll(track3);
      int pidProtonTr1 = selectorProton.getStatusTrackPIDAll(track1);
      int pidProtonTr2 = selectorProton.getStatusTrackPIDAll(track2);
      int pidProtonTr3 = selectorProton.getStatusTrackPIDAll(track3);

      /// Select Lc and Lcbar candidates
      if (pidProtonTr1 == TrackSelectorPID::Status::PIDAccepted && pidPionTr2 == TrackSelectorPID::Status::PIDAccepted && pidKaonTr3 == TrackSelectorPID::Status::PIDAccepted) {
        if (track1.sign() == 1 && track2.sign() == 1 && track3.sign() == -1) {
          CandLc = true;
          source = 1;
          kfpTrackPosPi = createKFPTrackFromTrack(track1);
          kfpTrackNegKa = createKFPTrackFromTrack(track2);
          TPCnSigmaPosPi = track1.tpcNSigmaPi();
          TPCnSigmaNegKa = track2.tpcNSigmaKa();
          TOFnSigmaPosPi = track1.tofNSigmaPi();
          TOFnSigmaNegKa = track2.tofNSigmaKa();
        } else if (track1.sign() == -1 && track2.sign() == 1) {
          CandD0bar = true;
          source = 2;
          kfpTrackNegPi = createKFPTrackFromTrack(track1);
          kfpTrackPosKa = createKFPTrackFromTrack(track2);
          TPCnSigmaNegPi = track1.tpcNSigmaPi();
          TPCnSigmaPosKa = track2.tpcNSigmaKa();
          TOFnSigmaNegPi = track1.tofNSigmaPi();
          TOFnSigmaPosKa = track2.tofNSigmaKa();
        } else {
          continue;
        }
      }
      if (pidKaonTr1 == TrackSelectorPID::Status::PIDAccepted && pidPionTr2 == TrackSelectorPID::Status::PIDAccepted) {
        if (track1.sign() == 1 && track2.sign() == -1) {
          CandD0bar = true;
          source = 2;
          kfpTrackNegPi = createKFPTrackFromTrack(track2);
          kfpTrackPosKa = createKFPTrackFromTrack(track1);
          TPCnSigmaNegPi = track2.tpcNSigmaPi();
          TPCnSigmaPosKa = track1.tpcNSigmaKa();
          TOFnSigmaNegPi = track2.tofNSigmaPi();
          TOFnSigmaPosKa = track1.tofNSigmaKa();
        } else if (track1.sign() == -1 && track2.sign() == 1) {
          CandD0 = true;
          source = 1;
          kfpTrackPosPi = createKFPTrackFromTrack(track2);
          kfpTrackNegKa = createKFPTrackFromTrack(track1);
          TPCnSigmaPosPi = track2.tpcNSigmaPi();
          TPCnSigmaNegKa = track1.tpcNSigmaKa();
          TOFnSigmaPosPi = track2.tofNSigmaPi();
          TOFnSigmaNegKa = track1.tofNSigmaKa();
        } else {
          continue;
        }
      }
      if (!CandD0 && !CandD0bar) {
        continue;
      }
      if (CandD0 && CandD0bar) {
        source = 3;
      }

      KFParticle KFPosPion(kfpTrackPosPi, 211);
      KFParticle KFNegPion(kfpTrackNegPi, 211);
      KFParticle KFPosKaon(kfpTrackPosKa, 321);
      KFParticle KFNegKaon(kfpTrackNegKa, 321);

      int NDaughters = 2;
      float cosThetaStar = 0;

      if (CandD0) {
        KFParticle KFDZero;
        const KFParticle* D0Daughters[2] = {&KFPosPion, &KFNegKaon};
        KFDZero.SetConstructMethod(2);
        KFDZero.Construct(D0Daughters, NDaughters);
        /// Apply daughter selection
        if (!isSelectedDaughters(KFPosPion, KFNegKaon, KFDZero, KFPV)) {
          continue;
        }
        /// Apply selection on geometrically reconstructed D0
        cosThetaStar = cosThetaStarFromKF(1, 421, 211, 321, KFPosPion, KFNegKaon);
        if (!isSelectedDoGeo(KFDZero, KFPV, cosThetaStar)) {
          continue;
        }
        /// Apply a topological constraint of the D0 to the PV.
        /// Parameters will be given at the primary vertex.
        KFParticle KFDZero_PV = KFDZero;
        KFDZero_PV.SetProductionVertex(KFPV);
        /// Transport the D0 after the topological constraint back to the decay vertex
        KFParticle KFDZero_DecayVtx = KFDZero_PV;
        KFDZero_DecayVtx.TransportToDecayVertex();
        if (applySelectionDoWithTopoConst) {
          /// Apply selection on D0 after constraint to the PV
          if (!isSelectedDoTopo(KFDZero_PV, KFPosPion, KFNegKaon, KFDZero_DecayVtx, KFPV)) {
            continue;
          }
        }

        writeVarTree(kfpTrackPosPi, kfpTrackNegKa, KFPosPion, KFNegKaon, KFDZero_PV, KFDZero, KFPV, KFDZero_DecayVtx, TPCnSigmaPosPi, TOFnSigmaPosPi, TPCnSigmaNegKa, TOFnSigmaNegKa, cosThetaStar, track1, source);
      }
      if (CandD0bar) {
        KFParticle KFDZeroBar;
        const KFParticle* D0BarDaughters[2] = {&KFNegPion, &KFPosKaon};
        KFDZeroBar.SetConstructMethod(2);
        KFDZeroBar.Construct(D0BarDaughters, NDaughters);
        /// Apply daughter selection
        if (!isSelectedDaughters(KFNegPion, KFPosKaon, KFDZeroBar, KFPV)) {
          continue;
        }
        /// Apply selection on geometrically reconstructed D0
        cosThetaStar = cosThetaStarFromKF(0, 421, 321, 211, KFPosKaon, KFNegPion);
        if (!isSelectedDoGeo(KFDZeroBar, KFPV, cosThetaStar)) {
          continue;
        }
        /// Apply a topological constraint of the D0 to the PV.
        /// Parameters will be given at the primary vertex.
        KFParticle KFDZeroBar_PV = KFDZeroBar;
        KFDZeroBar_PV.SetProductionVertex(KFPV);
        /// Transport the D0 after the topological constraint back to the decay vertex
        KFParticle KFDZeroBar_DecayVtx = KFDZeroBar_PV;
        KFDZeroBar_DecayVtx.TransportToDecayVertex();
        if (applySelectionDoWithTopoConst) {
          /// Apply selection on D0 after constraint to the PV
          if (!isSelectedDoTopo(KFDZeroBar_PV, KFNegPion, KFPosKaon, KFDZeroBar_DecayVtx, KFPV)) {
            continue;
          }
        }
        writeVarTree(kfpTrackNegPi, kfpTrackPosKa, KFNegPion, KFPosKaon, KFDZeroBar_PV, KFDZeroBar, KFPV, KFDZeroBar_DecayVtx, TPCnSigmaNegPi, TOFnSigmaNegPi, TPCnSigmaPosKa, TOFnSigmaPosKa, cosThetaStar, track1, source);
      }
    }
  }
  PROCESS_SWITCH(qaKFParticle, processData, "process data", true);

  /// Process function for MC
  using CollisionTableMC = soa::Join<CollisionTableData, aod::McCollisionLabels>;
  using CollisionTableDataMult = soa::Join<aod::Collisions, aod::Mults, aod::McCollisionLabels>;
  using TrackTableMC = soa::Join<TrackTableData, aod::McTrackLabels>;
  Preslice<aod::McCollisionLabels> perMcCollision = aod::mccollisionlabel::mcCollisionId;
  void processMC(CollisionTableMC::iterator const& collision, CollisionTableMC const& collisions, soa::Filtered<TrackTableMC> const& tracks, aod::McParticles const& mcParticles, aod::McCollisions const& mcCollisions, aod::BCsWithTimestamps const&)
  {
    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    if (runNumber != bc.runNumber()) {
      initMagneticFieldCCDB(bc, runNumber, ccdb, isRun3 ? ccdbPathGrpMag : ccdbPathGrp, lut, isRun3);
      magneticField = o2::base::Propagator::Instance()->getNominalBz();
      /// Set magnetic field for KF vertexing
      #ifdef HomogeneousField
      KFParticle::SetField(magneticField);
      #endif
    }
    /// Remove Collisions without a MC Collision
    if (!collision.has_mcCollision()) {
      return;
    }
    /// Apply event selection
    if (!isSelectedCollision(collision)) {
      return;
    }
    /// set KF primary vertex
    KFPVertex kfpVertex = createKFPVertexFromCollision(collision);
    KFParticle KFPV(kfpVertex);

    KFPVertex kfpVertexDefault = createKFPVertexFromCollision(collision);
    KFParticle KFPVDefault(kfpVertexDefault);

    TrackSelectorPID selectorPion(kPiPlus);
    selectorPion.setRangePtTPC(ptPidTpcMinPi, ptPidTpcMaxPi);
    selectorPion.setRangeNSigmaTPC(-nSigmaTpcMaxPi, nSigmaTpcMaxPi);
    selectorPion.setRangeNSigmaTPCCondTOF(-nSigmaTpcCombinedMaxPi, nSigmaTpcCombinedMaxPi);
    selectorPion.setRangePtTOF(ptPidTofMinPi, ptPidTofMaxPi);
    selectorPion.setRangeNSigmaTOF(-nSigmaTofMaxPi, nSigmaTofMaxPi);
    selectorPion.setRangeNSigmaTOFCondTPC(-nSigmaTofCombinedMaxPi, nSigmaTofCombinedMaxPi);

    TrackSelectorPID selectorKaon(kKPlus);
    selectorPion.setRangePtTPC(ptPidTpcMinKa, ptPidTpcMaxKa);
    selectorPion.setRangeNSigmaTPC(-nSigmaTpcMaxKa, nSigmaTpcMaxKa);
    selectorPion.setRangeNSigmaTPCCondTOF(-nSigmaTpcCombinedMaxKa, nSigmaTpcCombinedMaxKa);
    selectorPion.setRangePtTOF(ptPidTofMinKa, ptPidTofMaxKa);
    selectorPion.setRangeNSigmaTOF(-nSigmaTofMaxKa, nSigmaTofMaxKa);
    selectorPion.setRangeNSigmaTOFCondTPC(-nSigmaTofCombinedMaxKa, nSigmaTofCombinedMaxKa);

    for (auto& [track1, track2] : combinations(soa::CombinationsStrictlyUpperIndexPolicy(tracks, tracks))) {

      /// Apply single track selection
      if (!isSelectedTracks(track1, track2)) {
        continue;
      }

      /// Check whether the track was assigned to the true MC PV
      auto particle1 = track1.mcParticle();
      auto particle2 = track2.mcParticle();
      auto collMC = particle1.mcCollision();
      auto mcCollID_recoColl = track1.collision_as<CollisionTableMC>().mcCollisionId();
      auto mcCollID_particle = particle1.mcCollisionId();
      bool indexMatchOK = (mcCollID_recoColl == mcCollID_particle);
      if (!indexMatchOK) {
        const auto matchedCollisions = collisions.sliceBy(perMcCollision, collMC.globalIndex());
        int i = 0;
        std::array<float, 5> dcaZ{100, 100, 100, 100, 100};
        float min = 100;
        for (auto matchedCollision : matchedCollisions) {
          dcaZ[i] = abs(matchedCollision.posZ() - collMC.posZ());
          if (i == 0) {
            min = dcaZ[i];
          }
          if (i > 0) {
            if (dcaZ[i] < dcaZ[i - 1]) {
              min = dcaZ[i];
            }
          }

          i = i + 1;
        }
        if (min > 10.) {
        }
        int j = 0;
        for (auto matchedCollision : matchedCollisions) {
          if (i == 1) {
            kfpVertex = createKFPVertexFromCollision(matchedCollision);
            KFParticle KFPVNew(kfpVertex);
            KFPV = KFPVNew;
          }
          if (i > 1) {
            if (abs(matchedCollision.posZ() - collMC.posZ()) == min) {
              kfpVertex = createKFPVertexFromCollision(matchedCollision);
              KFParticle KFPVNew(kfpVertex);
              KFPV = KFPVNew;
            }
          }
          j = j + 1;
        }
      }

      /// Remove fake tracks
      if (!track1.has_mcParticle() || !track2.has_mcParticle()) {
        continue;
      }
      /// Remove unmatched tracks
      if (track1.collisionId() <= 0 || track1.collisionId() <= 0) {
        continue;
      }

      int8_t sign = 0;
      int8_t flag = RecoDecay::OriginType::None;
      int sourceD0 = 0;
      int sourceD0Bar = 0;

      auto indexRec = RecoDecay::getMatchedMCRec(mcParticles, std::array{track1, track2}, 421, array{211, -321}, true, &sign);
      if (indexRec > -1) {
        auto particle = mcParticles.rawIteratorAt(indexRec);
        flag = RecoDecay::getCharmHadronOrigin(mcParticles, particle);
        if (flag == RecoDecay::OriginType::Prompt) {
          sourceD0 |= kPrompt;
          sourceD0Bar |= kPrompt;
        } else if (flag == RecoDecay::OriginType::NonPrompt) {
          sourceD0 |= kNonPrompt;
          sourceD0Bar |= kNonPrompt;
        }
        if (flag != RecoDecay::OriginType::Prompt && flag != RecoDecay::OriginType::NonPrompt) {
          continue;
        }
      } else {
        continue;
      }

      KFPTrack kfpTrackPosPi;
      KFPTrack kfpTrackNegPi;
      KFPTrack kfpTrackPosKa;
      KFPTrack kfpTrackNegKa;

      bool CandD0 = false;
      bool CandD0bar = false;
      float TPCnSigmaPosPi = 0;
      float TPCnSigmaNegPi = 0;
      float TPCnSigmaPosKa = 0;
      float TPCnSigmaNegKa = 0;
      float TOFnSigmaPosPi = 0;
      float TOFnSigmaNegPi = 0;
      float TOFnSigmaPosKa = 0;
      float TOFnSigmaNegKa = 0;

      int pidKaonTr1 = selectorKaon.getStatusTrackPIDAll(track1);
      int pidKaonTr2 = selectorKaon.getStatusTrackPIDAll(track2);
      int pidPionTr1 = selectorPion.getStatusTrackPIDAll(track1);
      int pidPionTr2 = selectorPion.getStatusTrackPIDAll(track2);


      int pdgMother = mcParticles.rawIteratorAt(indexRec - mcParticles.offset()).pdgCode();
      /// Select D0 and D0bar candidates
      if (pidPionTr1 == TrackSelectorPID::Status::PIDAccepted && pidKaonTr2 == TrackSelectorPID::Status::PIDAccepted) {
        if (track1.sign() == 1 && track2.sign() == -1) {
          CandD0 = true;
          particle1 = track1.mcParticle();
          particle2 = track2.mcParticle();
          if (pdgMother == -421) {
            sourceD0 |= kReflection;
          }
          kfpTrackPosPi = createKFPTrackFromTrack(track1);
          kfpTrackNegKa = createKFPTrackFromTrack(track2);
          TPCnSigmaPosPi = track1.tpcNSigmaPi();
          TPCnSigmaNegKa = track2.tpcNSigmaKa();
          TOFnSigmaPosPi = track1.tofNSigmaPi();
          TOFnSigmaNegKa = track2.tofNSigmaKa();
        } else if (track1.sign() == -1 && track2.sign() == 1) {
          CandD0bar = true;
          if (pdgMother == 421) {
            sourceD0Bar |= kReflection;
          }
          kfpTrackNegPi = createKFPTrackFromTrack(track1);
          kfpTrackPosKa = createKFPTrackFromTrack(track2);
          TPCnSigmaNegPi = track1.tpcNSigmaPi();
          TPCnSigmaPosKa = track2.tpcNSigmaKa();
          TOFnSigmaNegPi = track1.tofNSigmaPi();
          TOFnSigmaPosKa = track2.tofNSigmaKa();
        } else {
          continue;
        }
      }
      if (pidKaonTr1 == TrackSelectorPID::Status::PIDAccepted && pidPionTr2 == TrackSelectorPID::Status::PIDAccepted) {
        if (track1.sign() == 1 && track2.sign() == -1) {
          CandD0bar = true;
          if (pdgMother == 421) {
            sourceD0Bar |= kReflection;
          }
          kfpTrackNegPi = createKFPTrackFromTrack(track2);
          kfpTrackPosKa = createKFPTrackFromTrack(track1);
          TPCnSigmaNegPi = track2.tpcNSigmaPi();
          TPCnSigmaPosKa = track1.tpcNSigmaKa();
          TOFnSigmaNegPi = track2.tofNSigmaPi();
          TOFnSigmaPosKa = track1.tofNSigmaKa();
        } else if (track1.sign() == -1 && track2.sign() == 1) {
          CandD0 = true;
          if (pdgMother == -421) {
            sourceD0 |= kReflection;
          }
          kfpTrackPosPi = createKFPTrackFromTrack(track2);
          kfpTrackNegKa = createKFPTrackFromTrack(track1);
          TPCnSigmaPosPi = track2.tpcNSigmaPi();
          TPCnSigmaNegKa = track1.tpcNSigmaKa();
          TOFnSigmaPosPi = track2.tofNSigmaPi();
          TOFnSigmaNegKa = track1.tofNSigmaKa();
        } else {
          continue;
        }
      }
      if (!CandD0 && !CandD0bar) {
        continue;
      }

      KFParticle KFPosPion(kfpTrackPosPi, 211);
      KFParticle KFNegPion(kfpTrackNegPi, 211);
      KFParticle KFPosKaon(kfpTrackPosKa, 321);
      KFParticle KFNegKaon(kfpTrackNegKa, 321);

      int NDaughters = 2;
      float cosThetaStar = 0;

      if (CandD0) {
        KFParticle KFDZero;
        const KFParticle* D0Daughters[2] = {&KFPosPion, &KFNegKaon};
        KFDZero.SetConstructMethod(2);
        KFDZero.Construct(D0Daughters, NDaughters);
        /// Apply daughter selection
        if (!isSelectedDaughters(KFPosPion, KFNegKaon, KFDZero, KFPV)) {
          continue;
        }
        /// Apply selection on geometrically reconstructed D0
        cosThetaStar = cosThetaStarFromKF(1, 421, 211, 321, KFPosPion, KFNegKaon);
        if (!isSelectedDoGeo(KFDZero, KFPV, cosThetaStar)) {
          continue;
        }
        /// Apply a topological constraint of the D0 to the PV.
        /// Parameters will be given at the primary vertex.
        KFParticle KFDZero_PV = KFDZero;
        KFDZero_PV.SetProductionVertex(KFPV);
        /// Transport the D0 after the topological constraint back to the decay vertex
        KFParticle KFDZero_DecayVtx = KFDZero_PV;
        KFDZero_DecayVtx.TransportToDecayVertex();
        if (applySelectionDoWithTopoConst) {
          /// Apply selection on D0 after constraint to the PV
          if (!isSelectedDoTopo(KFDZero_PV, KFPosPion, KFNegKaon, KFDZero_DecayVtx, KFPV)) {
            continue;
          }
        }
        writeVarTree(kfpTrackPosPi, kfpTrackNegKa, KFPosPion, KFNegKaon, KFDZero_PV, KFDZero, KFPV, KFDZero_DecayVtx, TPCnSigmaPosPi, TOFnSigmaPosPi, TPCnSigmaNegKa, TOFnSigmaNegKa, cosThetaStar, track1, sourceD0);
      }
      if (CandD0bar) {
        KFParticle KFDZeroBar;
        const KFParticle* D0BarDaughters[2] = {&KFNegPion, &KFPosKaon};
        KFDZeroBar.SetConstructMethod(2);
        KFDZeroBar.Construct(D0BarDaughters, NDaughters);
        /// Apply daughter selection
        if (!isSelectedDaughters(KFNegPion, KFPosKaon, KFDZeroBar, KFPV)) {
          continue;
        }
        /// Apply selection on geometrically reconstructed D0
        cosThetaStar = cosThetaStarFromKF(0, 421, 321, 211, KFPosKaon, KFNegPion);
        if (!isSelectedDoGeo(KFDZeroBar, KFPV, cosThetaStar)) {
          continue;
        }
        /// Apply a topological constraint of the D0 to the PV.
        /// Parameters will be given at the primary vertex.
        KFParticle KFDZeroBar_PV = KFDZeroBar;
        KFDZeroBar_PV.SetProductionVertex(KFPV);
        /// Transport the D0 after the topological constraint back to the decay vertex
        KFParticle KFDZeroBar_DecayVtx = KFDZeroBar_PV;
        KFDZeroBar_DecayVtx.TransportToDecayVertex();
        if (applySelectionDoWithTopoConst) {
          /// Apply selection on D0 after constraint to the PV
          if (!isSelectedDoTopo(KFDZeroBar_PV, KFNegPion, KFPosKaon, KFDZeroBar_DecayVtx, KFPV)) {
            continue;
          }
        }
        writeVarTree(kfpTrackPosPi, kfpTrackNegKa, KFPosPion, KFNegKaon, KFDZeroBar_PV, KFDZeroBar, KFPV, KFDZeroBar_DecayVtx, TPCnSigmaNegPi, TOFnSigmaNegPi, TPCnSigmaPosKa, TOFnSigmaPosKa, cosThetaStar, track1, sourceD0Bar);
      }
    }
  }
  PROCESS_SWITCH(qaKFParticle, processMC, "process mc", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<qaKFParticle>(cfgc)};
}
