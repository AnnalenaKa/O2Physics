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
/// \file   pidTPCFull.cxx
/// \author Annalena Kalteyer annalena.sophie.kalteyer@cern.ch
/// \author Christian Sonnabend christian.sonnabend@cern.ch
/// \author Nicolò Jacazio nicolo.jacazio@cern.ch
/// \brief  Task to produce PID tables for TPC split for each particle.
///         Only the tables for the mass hypotheses requested are filled, the others are sent empty.
///         QA histograms for the TPC PID can be produced by adding `--add-qa 1` to the workflow
///

// ROOT includes
#include "TFile.h"

// O2 includes
#include "Framework/AnalysisTask.h"
#include "ReconstructionDataFormats/Track.h"
#include "ReconstructionDataFormats/TrackParametrization.h"
#include "ReconstructionDataFormats/DCA.h"
#include <CCDB/BasicCCDBManager.h>
#include "DataFormatsParameters/GRPObject.h"
#include "DetectorsBase/GeometryManager.h"
#include "DetectorsBase/Propagator.h"
#include "Common/Core/PID/PIDResponse.h"
#include "Common/Core/PID/TPCPIDResponse.h"
#include "Common/Core/trackUtilities.h"
#include "Framework/AnalysisDataModel.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "TableHelper.h"
#include "Common/TableProducer/PID/pidTPCML.h"
#include "DPG/Tasks/qaPIDTPC.h"

#include<iostream>
using namespace std;

using namespace o2;
using namespace o2::framework;
using namespace o2::pid;
using namespace o2::pid::tpc;
using namespace o2::framework::expressions;
using namespace o2::track;

namespace o2
{
namespace analysis
{
namespace trackpropagation
{
constexpr long run3grp_timestamp = (1619781650000 + 1619781529000) / 2;
} // namespace trackpropagation
} // namespace analysis
namespace aod
{
namespace tpc
{
DECLARE_SOA_COLUMN(TPCMomentumDiff, tpcMomentumDiff, float);
DECLARE_SOA_COLUMN(TPCMomentumShift, tpcMomentumShift, float);
DECLARE_SOA_COLUMN(NormMultTPC, normMultTPC, float);
DECLARE_SOA_COLUMN(NClustersTPC, nClustersTPC, float);
DECLARE_SOA_COLUMN(ShiftX, shiftx, float);

DECLARE_SOA_COLUMN(CID, cid, float);
DECLARE_SOA_COLUMN(PID, pid, float);
DECLARE_SOA_COLUMN(TrackType, tracktype, float);
DECLARE_SOA_COLUMN(HasITS, hasITS, bool);                                                //! Track has the ITS
DECLARE_SOA_COLUMN(HasTPC, hasTPC, bool);                                                //! Track has the TPC
DECLARE_SOA_COLUMN(HasTRD, hasTRD, bool);                                                //! Track has the TRD
DECLARE_SOA_COLUMN(HasTOF, hasTOF, bool);  

DECLARE_SOA_COLUMN(X, x, float);
DECLARE_SOA_COLUMN(Y, y, float);
DECLARE_SOA_COLUMN(Z, z, float);
DECLARE_SOA_COLUMN(Alpha, alpha, float);
DECLARE_SOA_COLUMN(SNP, snp, float);
DECLARE_SOA_COLUMN(TGL, tgl, float);
DECLARE_SOA_COLUMN(Q2Pt, q2pt, float);
DECLARE_SOA_COLUMN(Pt, pt, float);
DECLARE_SOA_COLUMN(P, p, float);
DECLARE_SOA_COLUMN(Eta, eta, float);
DECLARE_SOA_COLUMN(Phi, phi, float);

DECLARE_SOA_COLUMN(X2, x2, float);
DECLARE_SOA_COLUMN(Y2, y2, float);
DECLARE_SOA_COLUMN(Z2, z2, float);
DECLARE_SOA_COLUMN(Alpha2, alpha2, float);
DECLARE_SOA_COLUMN(SNP2, snp2, float);
DECLARE_SOA_COLUMN(TGL2, tgl2, float);
DECLARE_SOA_COLUMN(Q2Pt2, q2pt2, float);
DECLARE_SOA_COLUMN(Pt2, pt2, float);
DECLARE_SOA_COLUMN(P2, p2, float);
DECLARE_SOA_COLUMN(Eta2, eta2, float);
DECLARE_SOA_COLUMN(Phi2, phi2, float);

DECLARE_SOA_COLUMN(BG, bg, float);



} // namespace tpc 
DECLARE_SOA_TABLE(TPCTree, "AOD", "TPCTREE",
                  o2::aod::track::TPCInnerParam,
                  o2::aod::track::TPCTgl,
                  o2::aod::track::TPCSigned1Pt,
                  tpc::TPCMomentumDiff,
                  tpc::TPCMomentumShift,
                  tpc::NormMultTPC,
                  tpc::NClustersTPC,
                  o2::aod::track::TPCNClsdEdx,
                  o2::aod::track::TPCSignal,
                  tpc::BG,
                  o2::aod::track::Phi,
                  track::DcaXY,
                  track::DcaZ,
                  tpc::ShiftX,
                  tpc::CID,
                  tpc::PID,
                  tpc::TrackType,
                  tpc::HasITS,
                  tpc::HasTPC,
                  tpc::HasTRD,
                  tpc::HasTOF,
                  tpc::X,
                  tpc::Y,
                  tpc::Z,
                  tpc::Alpha,
                  tpc::SNP,
                  tpc::TGL,
                  tpc::Q2Pt,
                  tpc::Pt,
                  tpc::P,
                  tpc::Eta,
                  tpc::Phi,
                  tpc::X2,
                  tpc::Y2,
                  tpc::Z2,
                  tpc::Alpha2,
                  tpc::SNP2,
                  tpc::TGL2,
                  tpc::Q2Pt2,
                  tpc::Pt2,
                  tpc::P2,
                  tpc::Eta2,
                  tpc::Phi2);

} // namespace aod

} // namespace 02


void customize(std::vector<o2::framework::ConfigParamSpec>& workflowOptions)
{
  std::vector<ConfigParamSpec> options{{"add-qa", VariantType::Int, 0, {"Produce TPC PID QA histograms"}}};
  std::swap(workflowOptions, options);
}

#include "Framework/runDataProcessing.h"

/// Task to produce the response table
struct tpcPidFullPropTracks {
  using Trks = soa::Join<aod::Tracks, aod::TracksExtra, aod::FullTracks,  aod::TrackSelection, aod::TracksExtended>;
  using Coll = soa::Join<aod::Collisions, aod::Mults>;

  // Tables to produce
  Produces<o2::aod::pidTPCFullEl> tablePIDEl;
  Produces<o2::aod::pidTPCFullMu> tablePIDMu;
  Produces<o2::aod::pidTPCFullPi> tablePIDPi;
  Produces<o2::aod::pidTPCFullKa> tablePIDKa;
  Produces<o2::aod::pidTPCFullPr> tablePIDPr;
  Produces<o2::aod::pidTPCFullDe> tablePIDDe;
  Produces<o2::aod::pidTPCFullTr> tablePIDTr;
  Produces<o2::aod::pidTPCFullHe> tablePIDHe;
  Produces<o2::aod::pidTPCFullAl> tablePIDAl;

  Produces<o2::aod::TPCTree> rowTPCTree;
  // TPC PID Response
  o2::pid::tpc::Response response;
  o2::pid::tpc::Response* responseptr = nullptr;
  // Network correction for TPC PID response
  Network network;
  int runNumber = -1;

  // Input parameters
    // Input parameters
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  o2::base::Propagator::MatCorrType matCorr;
  o2::parameters::GRPObject* grpo = nullptr;
  o2::base::MatLayerCylSet* lut = nullptr;
  static constexpr float MaxSnp = 0.9;
  Configurable<std::string> paramfile{"param-file", "", "Path to the parametrization object, if emtpy the parametrization is not taken from file"};
  Configurable<std::string> url{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> ccdbPath{"ccdbPath", "Analysis/PID/TPC/Response", "Path of the TPC parametrization on the CCDB"};
  Configurable<long> ccdbTimestamp{"ccdb-timestamp", 0, "timestamp of the object used to query in CCDB the detector response. Exceptions: -1 gets the latest object, 0 gets the run dependent timestamp"};
  // Parameters for loading network from a file / downloading the file
  Configurable<int> useNetworkCorrection{"useNetworkCorrection", 0, "Using the network correction for the TPC dE/dx signal"};
  Configurable<int> downloadNetworkFromAlien{"downloadNetworkFromAlien", 0, "Download network from AliEn (1) or use a local file (filepath must be provided by --networkPathLocally /path/to/file) (0)"};
  Configurable<std::string> networkPathAlien{"networkPathAlien", "alien:///alice/cern.ch/user/c/csonnabe/tpc_network_testing/net_onnx_0.onnx", "Path to .onnx file containing the network on AliEn"};
  Configurable<std::string> networkPathLocally{"networkPathLocally", "network.onnx", "Path to local .onnx file containing the network"};
  Configurable<std::string> lutPath{"lutPath", "GLO/Param/MatLUT", "Path of the Lut parametrization"};
  Configurable<std::string> geoPath{"geoPath", "GLO/Config/GeometryAligned", "Path of the geometry file"};
  Configurable<std::string> grpPath{"grpPath", "GLO/GRP/GRP", "Path of the grp file"};
  // Configuration flags to include and exclude particle hypotheses
  Configurable<int> pidEl{"pid-el", -1, {"Produce PID information for the Electron mass hypothesis, overrides the automatic setup: the corresponding table can be set off (0) or on (1)"}};
  Configurable<int> pidMu{"pid-mu", -1, {"Produce PID information for the Muon mass hypothesis, overrides the automatic setup: the corresponding table can be set off (0) or on (1)"}};
  Configurable<int> pidPi{"pid-pi", -1, {"Produce PID information for the Pion mass hypothesis, overrides the automatic setup: the corresponding table can be set off (0) or on (1)"}};
  Configurable<int> pidKa{"pid-ka", -1, {"Produce PID information for the Kaon mass hypothesis, overrides the automatic setup: the corresponding table can be set off (0) or on (1)"}};
  Configurable<int> pidPr{"pid-pr", -1, {"Produce PID information for the Proton mass hypothesis, overrides the automatic setup: the corresponding table can be set off (0) or on (1)"}};
  Configurable<int> pidDe{"pid-de", -1, {"Produce PID information for the Deuterons mass hypothesis, overrides the automatic setup: the corresponding table can be set off (0) or on (1)"}};
  Configurable<int> pidTr{"pid-tr", -1, {"Produce PID information for the Triton mass hypothesis, overrides the automatic setup: the corresponding table can be set off (0) or on (1)"}};
  Configurable<int> pidHe{"pid-he", -1, {"Produce PID information for the Helium3 mass hypothesis, overrides the automatic setup: the corresponding table can be set off (0) or on (1)"}};
  Configurable<int> pidAl{"pid-al", -1, {"Produce PID information for the Alpha mass hypothesis, overrides the automatic setup: the corresponding table can be set off (0) or on (1)"}};

  // Paramatrization configuration
  bool useCCDBParam = false;

  void init(o2::framework::InitContext& initContext)
  {
    // Checking the tables are requested in the workflow and enabling them
    auto enableFlag = [&](const std::string particle, Configurable<int>& flag) {
      const std::string table = "pidTPCFull" + particle;
      if (isTableRequiredInWorkflow(initContext, table)) {
        if (flag < 0) {
          flag.value = 1;
          LOG(info) << "Auto-enabling table: " + table;
        } else if (flag > 0) {
          flag.value = 1;
          LOG(info) << "Table enabled: " + table;
        } else {
          LOG(info) << "Table disabled: " + table;
        }
      }
    };
    enableFlag("El", pidEl);
    enableFlag("Mu", pidMu);
    enableFlag("Pi", pidPi);
    enableFlag("Ka", pidKa);
    enableFlag("Pr", pidPr);
    enableFlag("De", pidDe);
    enableFlag("Tr", pidTr);
    enableFlag("He", pidHe);
    enableFlag("Al", pidAl);

    ccdb->setURL(url.value);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    lut = o2::base::MatLayerCylSet::rectifyPtrFromFile(ccdb->get<o2::base::MatLayerCylSet>(lutPath));
    matCorr = o2::base::Propagator::MatCorrType::USEMatCorrLUT;
    if (!o2::base::GeometryManager::isGeometryLoaded()) {
      ccdb->get<TGeoManager>(geoPath);
    }

    const TString fname = paramfile.value;
    if (fname != "") { // Loading the parametrization from file
      LOGP(info, "Loading TPC response from file {}", fname);
      try {
        std::unique_ptr<TFile> f(TFile::Open(fname, "READ"));
        f->GetObject("Response", responseptr);
        response.SetParameters(responseptr);
      } catch (...) {
        LOGF(fatal, "Loading the TPC PID Response from file {} failed!", fname);
      };
    } else {
      useCCDBParam = true;
      const std::string path = ccdbPath.value;
      const auto time = ccdbTimestamp.value;
      ccdb->setURL(url.value);
      ccdb->setTimestamp(time);
      ccdb->setCaching(true);
      ccdb->setLocalObjectValidityChecking();
      ccdb->setCreatedNotAfter(std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count());
      response.SetParameters(ccdb->getForTimeStamp<o2::pid::tpc::Response>(path, time));
      LOGP(info, "Loading TPC response from CCDB, using path: {} for ccdbTimestamp {}", path, time);
      response.PrintAll();
    }


    if (!useNetworkCorrection) {
      return;
    } else {
      Network temp_net(networkPathLocally.value,
                       downloadNetworkFromAlien.value,
                       networkPathAlien.value,
                       true);
      network = temp_net;
    }
  }

  void initCCDB(aod::BCsWithTimestamps::iterator const& bc)
  {
    if (runNumber == bc.runNumber()) {
      return;
    }
    grpo = ccdb->getForTimeStamp<o2::parameters::GRPObject>(grpPath, bc.timestamp());
    LOGF(info, "Setting magnetic field to %d kG for run %d from its GRP CCDB object", grpo->getNominalL3Field(), bc.runNumber());
    o2::base::Propagator::initFieldFromGRP(grpo);
    o2::base::Propagator::Instance()->setMatLUT(lut);
    runNumber = bc.runNumber();
  }

  void process(Coll const& collisions, Trks const& tracks,
               aod::BCsWithTimestamps const& bcs)
  {
    if (bcs.size() == 0) {
      return;
    }
    initCCDB(bcs.begin());

    auto reserveTable = [&tracks](const Configurable<int>& flag, auto& table) {
      if (flag.value != 1) {
        return;
      }
      table.reserve(tracks.size());
    };
    // Prepare memory for enabled tables
    reserveTable(pidEl, tablePIDEl);
    reserveTable(pidMu, tablePIDMu);
    reserveTable(pidPi, tablePIDPi);
    reserveTable(pidKa, tablePIDKa);
    reserveTable(pidPr, tablePIDPr);
    reserveTable(pidDe, tablePIDDe);
    reserveTable(pidTr, tablePIDTr);
    reserveTable(pidHe, tablePIDHe);
    reserveTable(pidAl, tablePIDAl);

    auto reserveTableTPCTree = [&tracks](auto& tabletree) {
      tabletree.reserve(tracks.size());
    };
    reserveTableTPCTree(rowTPCTree);

    std::vector<float> network_prediction;

    if (useNetworkCorrection) {

      auto start_overhead = std::chrono::high_resolution_clock::now();
      std::vector<float> track_properties;
      for (int i = 0; i < 9; i++) { // Loop over particle number for which network correction is used
        for (auto const& trk : tracks) {
          std::vector<float> net_tensor = network.createInputFromTrack(trk, i);
          for (auto value : net_tensor) {
            track_properties.push_back(value);
          }
        }
      }
      const unsigned long track_prop_size = tracks.size() * 9;
      auto stop_overhead = std::chrono::high_resolution_clock::now();
      float duration_overhead = std::chrono::duration<float, std::ratio<1, 1000000000>>(stop_overhead - start_overhead).count();
      float time_per_track_overhead = duration_overhead / track_prop_size; // There are n (typically n=7) variables in each track which are being extracted in track_properties. Each network evaluation takes time_per_track_overhead/9 nano-seconds
      LOG(info) << "Time per track (overhead): " << time_per_track_overhead << "ns ; Overhead total: " << duration_overhead / 1000000000 << "s";

      auto start_network = std::chrono::high_resolution_clock::now();
      float* output_network = network.evalNetwork(track_properties);
      for (unsigned long i = 0; i < track_prop_size; i++) {
        network_prediction.push_back(output_network[i]);
      }
      track_properties.clear();
      auto stop_network = std::chrono::high_resolution_clock::now();
      float duration_network = std::chrono::duration<float, std::ratio<1, 1000000000>>(stop_network - start_network).count();
      float time_per_track_net = duration_network / track_prop_size;
      LOG(info) << "Time per track (net): " << time_per_track_net << "ns ; Network total: " << duration_network / 1000000000 << "s"; // The time per track but with 9 particle mass hypotheses: So actual time per track is (time_per_track_net / 9)

      // Uncomment if you want to check example-outputs of the netwwork:
      // for(int i=0; i<100; i++){
      //   LOG(info) << "Output " << i << ": " << network_prediction[i] << " ; Input: [" << track_properties[7*i + 0] << ", " << track_properties[7*i + 1] << ", " << track_properties[7*i + 2] << ", " << track_properties[7*i + 3] << ", " << track_properties[7*i + 4] << ", " << track_properties[7*i + 5] << ", " << track_properties[7*i + 6] << "]";
      // }
    }

    int lastCollisionId = -1; // Last collision ID analysed
    unsigned long count_tracks = 0;
    const int tracks_size = tracks.size();

    for (auto const& trk : tracks) {
      // Loop on Tracks
      if (useCCDBParam && ccdbTimestamp.value == 0 && trk.has_collision() && trk.collisionId() != lastCollisionId) { // Updating parametrization only if the initial timestamp is 0
        lastCollisionId = trk.collisionId();
        const auto& bc = collisions.iteratorAt(trk.collisionId()).bc_as<aod::BCsWithTimestamps>();
        response.SetParameters(ccdb->getForTimeStamp<o2::pid::tpc::Response>(ccdbPath.value, bc.timestamp()));
      }
      // Propagate tracks to the inner wall of the TPC
      float xtogo = 0;
      float Bz = o2::base::Propagator::Instance()->getNominalBz();
      auto trackParInit = getTrackPar(trk);
      
      auto trackPar = getTrackPar(trk);

      std::array<float, 9> Probability;
      Probability[0] = response.ComputeTPCProbability(collisions.iteratorAt(trk.collisionId()), trk, o2::track::PID::Electron);
      Probability[1] = response.ComputeTPCProbability(collisions.iteratorAt(trk.collisionId()), trk, o2::track::PID::Muon);
      Probability[2] = response.ComputeTPCProbability(collisions.iteratorAt(trk.collisionId()), trk, o2::track::PID::Pion);
      Probability[3] = response.ComputeTPCProbability(collisions.iteratorAt(trk.collisionId()), trk, o2::track::PID::Kaon);
      Probability[4] = response.ComputeTPCProbability(collisions.iteratorAt(trk.collisionId()), trk, o2::track::PID::Proton);
      Probability[5] = response.ComputeTPCProbability(collisions.iteratorAt(trk.collisionId()), trk, o2::track::PID::Deuteron);
      Probability[6] = response.ComputeTPCProbability(collisions.iteratorAt(trk.collisionId()), trk, o2::track::PID::Triton);
      Probability[7] = response.ComputeTPCProbability(collisions.iteratorAt(trk.collisionId()), trk, o2::track::PID::Helium3);
      Probability[8] = response.ComputeTPCProbability(collisions.iteratorAt(trk.collisionId()), trk, o2::track::PID::Alpha);

      if(!trk.hasTPC()){trackPar.setPID(o2::track::PID::Pion);}
      // find max probability
      Float_t max=0.,min=1.e9;
      Int_t pid=-1;
      for (Int_t i=0; i<9; i++) {
        if (Probability[i]>max) {pid=i; max=Probability[i];}
        if (Probability[i]<min) min=Probability[i];
      }
      if(pid==-1) {trackPar.setPID(o2::track::PID::Pion);}
      if (pid==0) { // dE/dx "crossing points" in the TPC
        Double_t p = trackPar.getP();
        if ((p>0.38)&&(p<0.48)) {
          if (Probability[0]<Probability[3]*10.) {trackPar.setPID(o2::track::PID::Kaon);}
          else {trackPar.setPID(o2::track::PID::Electron);}
        }
        else if ((p>0.75)&&(p<0.85)) {
          if (Probability[0]<Probability[4]*10.) {trackPar.setPID(o2::track::PID::Proton);}
          else {trackPar.setPID(o2::track::PID::Electron);}
        }
      }

      if(pid==1) {trackPar.setPID(o2::track::PID::Muon);}
      if(pid==2) {trackPar.setPID(o2::track::PID::Pion);}
      if(pid==3) {trackPar.setPID(o2::track::PID::Kaon);}
      if(pid==4) {trackPar.setPID(o2::track::PID::Proton);}
      if(pid==5) {trackPar.setPID(o2::track::PID::Deuteron);}
      if(pid==6) {trackPar.setPID(o2::track::PID::Triton);}
      if(pid==7) {trackPar.setPID(o2::track::PID::Helium3);}
      if(pid==8) {trackPar.setPID(o2::track::PID::Alpha);}

      if (min>=max) {trackPar.setPID(o2::track::PID::Pion);}

      if (!trackPar.getXatLabR(83., xtogo, Bz, DirOutward ) ||
      !o2::base::Propagator::Instance()->PropagateToXBxByBz(trackPar, xtogo, MaxSnp, 10., matCorr)) {
        LOG(debug) << "Propagation to inner TPC boundary X=" << xtogo << " failed.";
      }

      // Check and fill enabled tables
      auto makeTable = [&trk, &collisions, &network_prediction, &count_tracks, &tracks_size, this](const Configurable<int>& flag, auto& table, const o2::track::PID::ID pid) {
        if (flag.value != 1) {
          return;
        }

        if (useNetworkCorrection) {
          table(response.GetExpectedSigma(collisions.iteratorAt(trk.collisionId()), trk, pid),
                (trk.tpcSignal() - (network_prediction[count_tracks + tracks_size * pid]) * response.GetExpectedSignal(trk, pid)) / response.GetExpectedSigma(collisions.iteratorAt(trk.collisionId()), trk, pid));
        } else {
          table(response.GetExpectedSigma(collisions.iteratorAt(trk.collisionId()), trk, pid),
                response.GetNumberOfSigma(collisions.iteratorAt(trk.collisionId()), trk, pid));
        }
      };

      // const o2::pid::tpc::Response& response;
      makeTable(pidEl, tablePIDEl, o2::track::PID::Electron);
      makeTable(pidMu, tablePIDMu, o2::track::PID::Muon);
      makeTable(pidPi, tablePIDPi, o2::track::PID::Pion);
      makeTable(pidKa, tablePIDKa, o2::track::PID::Kaon);
      makeTable(pidPr, tablePIDPr, o2::track::PID::Proton);
      makeTable(pidDe, tablePIDDe, o2::track::PID::Deuteron);
      makeTable(pidTr, tablePIDTr, o2::track::PID::Triton);
      makeTable(pidHe, tablePIDHe, o2::track::PID::Helium3);
      makeTable(pidAl, tablePIDAl, o2::track::PID::Alpha);

      auto fillTPCTable = [&trk, trackPar, trackParInit, Probability, &collisions, this] (auto& tableTPCTree)
      {

        const double ncl = trk.tpcNClsFound();
        const int multTPC = collisions.iteratorAt(trk.collisionId()).multTPC();
        const float pDiff = trk.tpcInnerParam() - trackPar.getP();
        const float pShift = trackParInit.getP() - trackPar.getP();
        const float shiftx = trackParInit.getX() - trackPar.getX();
        const float bg = trk.tpcInnerParam()/o2::track::pid_constants::sMasses[trackPar.getPID()];


        // if (pDiff < -1. &&  trk.tpcInnerParam()!=0.)
        // {
        //   cout << "Track parameters before propagation" << endl;
        //   trackParInit.printParam();
        //   cout << "Track parameters after propagation" << endl;
        //   trackPar.printParam(); 
        //   cout << "tpcInnerParam: " << trk.tpcInnerParam() << endl;
        //   cout << "flags: " << trk.flags() << endl;
        //   cout << "hasITS: " << trk.hasITS() << endl;
        //   cout << "hasTPC: " << trk.hasTPC() << endl;
        //   cout << "Probability[0] " << Probability[0] << endl;
        //   cout << "Probability[1] " << Probability[1] << endl;
        //   cout << "Probability[2] " << Probability[2] << endl;
        //   cout << "Probability[3] " << Probability[3] << endl;
        //   cout << "Probability[4] " << Probability[4] << endl;
        //   cout << "Probability[5] " << Probability[5] << endl;
        //   cout << "Probability[6] " << Probability[6] << endl;
        //   cout << "Probability[7] " << Probability[7] << endl;
        //   cout << "Probability[8] " << Probability[8] << endl; 

        // }
        

        

        tableTPCTree(
                    trk.tpcInnerParam(),
                    trk.tpcTgl(),
                    trk.tpcSigned1Pt(),
                    pDiff,
                    pShift,
                    multTPC / 11000.,
                    ncl,
                    trk.tpcNClsdEdx(),
                    trk.tpcSignal(),
                    bg,
                    trk.phi(),
                    trk.dcaXY(),
                    trk.dcaZ(),
                    shiftx,
                    trk.collisionId(),
                    trackParInit.getPID(),
                    trk.trackType(),
                    trk.hasITS(),
                    trk.hasTPC(),
                    trk.hasTRD(),
                    trk.hasTOF(),
                    trackParInit.getX(),
                    trackParInit.getY(),
                    trackParInit.getZ(),
                    trackParInit.getAlpha(),
                    trackParInit.getSnp(),
                    trackParInit.getTgl(),
                    trackParInit.getQ2Pt(),
                    trackParInit.getPt(),
                    trackParInit.getP(),
                    trackParInit.getEta(),
                    trackParInit.getPhi(),
                    trackPar.getX(),
                    trackPar.getY(),
                    trackPar.getZ(),
                    trackPar.getAlpha(),
                    trackPar.getSnp(),
                    trackPar.getTgl(),
                    trackPar.getQ2Pt(),
                    trackPar.getPt(),
                    trackPar.getP(),
                    trackPar.getEta(),
                    trackPar.getPhi());
        
      };
      fillTPCTable(rowTPCTree);

      count_tracks++;
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  auto workflow = WorkflowSpec{adaptAnalysisTask<tpcPidFullPropTracks>(cfgc)};
  if (cfgc.options().get<int>("add-qa")) {
    workflow.push_back(adaptAnalysisTask<tpcPidQa>(cfgc));
  }
  return workflow;
}