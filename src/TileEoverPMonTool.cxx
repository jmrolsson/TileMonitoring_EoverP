// ********************************************************************
//
// NAME:     TileEoverPMonTool.cxx
// PACKAGE:  TileMonitoring 
//
// AUTHOR:   Joakim Olsson (joakim.olsson@cern.ch)
//
// ********************************************************************

#include "TileMonitoring/TileEoverPMonTool.h"

// Tracks
#include "TrkTrack/Track.h"
#include "TrkParameters/TrackParameters.h"
#include "TrkExInterfaces/IExtrapolator.h"
#include "xAODTracking/TrackParticleContainer.h"
#include "xAODTracking/VertexContainer.h"
#include "xAODTracking/VertexAuxContainer.h"
#include "InDetTrackSelectionTool/IInDetTrackSelectionTool.h"

// Extrapolation to the calo
#include "TrkCaloExtension/CaloExtension.h"
#include "TrkCaloExtension/CaloExtensionCollection.h"
#include "TrkParametersIdentificationHelpers/TrackParametersIdHelper.h"
#include "CaloDetDescr/CaloDepthTool.h"

// Calo and cell information
#include "TileEvent/TileContainer.h"
#include "TileIdentifier/TileTBID.h"
#include "CaloEvent/CaloCellContainer.h"
#include "CaloTrackingGeometry/ICaloSurfaceHelper.h"
#include "TrkSurfaces/DiscSurface.h"
#include "GeoPrimitives/GeoPrimitives.h"
#include "CaloEvent/CaloClusterContainer.h"
#include "CaloEvent/CaloCluster.h"
#include "CaloUtils/CaloClusterSignalState.h"
#include "CaloEvent/CaloClusterCellLinkContainer.h"
#include "xAODCaloEvent/CaloClusterChangeSignalState.h"

// Other xAOD incudes
#include "xAODEventInfo/EventInfo.h"

#include <math.h>
#include <string>
#include <vector>
#include <utility>
#include <limits>

/*---------------------------------------------------------*/
TileEoverPMonTool::TileEoverPMonTool(const std::string & type, const std::string & name, const IInterface* parent) : 
  TileFatherMonTool(type, name, parent),
  m_prefix(""),
  m_eventInfoContainerName("EventInfo"),
  m_vertexContainerName("PrimaryVertices"),
  m_trackContainerName("InDetTrackParticles"),
  m_caloClusterContainerName("CaloCalTopoClusters"),
  m_cellContainerName("AllCalo"),
  m_extrapolator("Trk::Extrapolator"),
  m_trackExtrapolatorTool("Trk::ParticleCaloExtensionTool"),
  m_trackParametersIdHelper(new Trk::TrackParametersIdHelper),
  m_surfaceHelper("CaloSurfaceHelper/CaloSurfaceHelper"),
  m_trackSelectionTool("InDet::InDetTrackSelectionTool/TrackSelectionTool", this),
  m_trkIsoDRmax(0.4),
  m_trkClusMatchDRmax(0.2),
  m_doTrkExtrapolation(false),
  m_tileTBID(0) {
    /*---------------------------------------------------------*/

    declareInterface<IMonitorToolBase>(this);
    declareProperty("Prefix", m_prefix);
    declareProperty("EventContainerName", m_eventInfoContainerName);
    declareProperty("VertexContainerName", m_vertexContainerName);
    declareProperty("TrackContainerName", m_trackContainerName);
    declareProperty("CaloClusterContainerName", m_caloClusterContainerName);
    declareProperty("CellContainerName", m_cellContainerName);
    declareProperty("Extrapolator", m_extrapolator);
    declareProperty("TrackExtrapolatorTool", m_trackExtrapolatorTool);
    declareProperty("TrackSelectionTool", m_trackSelectionTool);
    declareProperty("TrackIsolationDeltaR", m_trkIsoDRmax);
    declareProperty("TrackClusterMatchingDeltaR", m_trkClusMatchDRmax);
    declareProperty("DoTrkExtrapolation", m_doTrkExtrapolation);
    declareProperty("DoBasicTrkSelection", m_doBasicTrkSelection);

    m_path = "/Tile/EoverP";

    m_first_event = true;
  }

/*---------------------------------------------------------*/
TileEoverPMonTool::~TileEoverPMonTool() {
/*---------------------------------------------------------*/

}

/*---------------------------------------------------------*/
StatusCode TileEoverPMonTool::initialize() {
/*---------------------------------------------------------*/

  ATH_MSG_INFO( "in TileEoverPMonTool::initialize()" );

  // Track extrapolation tool
  if (m_doTrkExtrapolation) {
    ATH_CHECK( m_extrapolator.retrieve() );
    ATH_CHECK( m_trackExtrapolatorTool.retrieve() );
    ATH_CHECK( m_surfaceHelper.retrieve() );
    // Get the test beam identifier for the MBTS
    ATH_CHECK( detStore()->retrieve(m_tileTBID) );
  }
  // Track selection tool
  if (m_doBasicTrkSelection) {
    ATH_CHECK( m_trackSelectionTool.retrieve() );
  }

  return TileFatherMonTool::initialize();
}

// BookHistogram is called at every event block, lumi block and run
// For the E/p monitoring we only want to book every run, therefore 
// the separate method 'bookEoverPHistograms' is used below
/*---------------------------------------------------------*/
StatusCode TileEoverPMonTool::bookHistograms() {
/*---------------------------------------------------------*/

  // ATH_MSG_INFO( "in TileEoverPMonTool::bookHistograms()" );

  return StatusCode::SUCCESS;
}


/*---------------------------------------------------------*/
StatusCode TileEoverPMonTool::fillHistograms() {
/*---------------------------------------------------------*/

  ATH_MSG_DEBUG( "in TileEoverPMonTool::fillHistograms()" );

  fillEvtInfo();
  if (m_first_event) {
    m_first_event = false;
    if (bookEoverPHistograms() != StatusCode::SUCCESS) {
      ATH_MSG_ERROR("Something wrong when booking histograms.");
    }
  }	

  // Retrieve primary vertex
  const xAOD::Vertex *pv(nullptr);;
  const xAOD::VertexContainer *vtxs(nullptr);
  CHECK(evtStore()->retrieve(vtxs, m_vertexContainerName));
  for( auto vtx_itr : *vtxs) {
    if(vtx_itr->vertexType() != xAOD::VxType::VertexType::PriVtx) continue;
    pv = vtx_itr;
  }

  // Retrieve tracks
  const xAOD::TrackParticleContainer* trks = 0;
  CHECK(evtStore()->retrieve(trks, m_trackContainerName));

  // Retrieve calo clusters
  const xAOD::CaloClusterContainer* clusters = 0;
  CHECK(evtStore()->retrieve(clusters, m_caloClusterContainerName));

  // Retrieve calo cells
  const CaloCellContainer *cells= 0;
  CHECK(evtStore()->retrieve(cells, m_cellContainerName));

  m_trk_n_all->Fill(trks->size()); 
  int trk_n_pass_trksel = 0;
  int trk_n_pass_extrapol = 0;
  int trk_n_pass_p = 0;
  int trk_n_pass_eta = 0;
  int trk_n_pass_isolation = 0;
  int trk_n_pass_larEmax = 0;
  int trk_n_pass_tileEfrac = 0;

  // Iterate over all tracks
  xAOD::TrackParticleContainer::const_iterator trk_itr = trks->begin();
  xAOD::TrackParticleContainer::const_iterator trk_end = trks->end();
  xAOD::TrackParticleContainer::const_iterator trk2_itr = trks->begin();
  xAOD::TrackParticleContainer::const_iterator trk2_end = trks->end();
  for( ; trk_itr != trk_end; ++trk_itr ) {
    const xAOD::TrackParticle* trk = (*trk_itr);

    // Apply basic track selections (defined in TileMon_jobOptions.py)
    if (m_doBasicTrkSelection && !m_trackSelectionTool->accept(*trk, pv)) continue;
    trk_n_pass_trksel++;

    /*---------E/p Decorators---------*/
    /*---------Calo Sample layer Variables---------*/
    // PreSamplerB=0, EMB1, EMB2, EMB3,       // LAr barrel
    // PreSamplerE, EME1, EME2, EME3,         // LAr EM endcap
    // HEC0, HEC1, HEC2, HEC3,                // Hadronic end cap cal.
    // TileBar0, TileBar1, TileBar2,          // Tile barrel
    // TileGap1, TileGap2, TileGap3,          // Tile gap (ITC & scint)
    // TileExt0, TileExt1, TileExt2,          // Tile extended barrel
    // FCAL0, FCAL1, FCAL2,                   // Forward EM endcap (excluded)
    // Unknown

    // Extrapolate tracks to calorimeter
    if (m_doTrkExtrapolation) {
      // TODO Implement track extrapolation 
      // (essentially just need to copy some code from the E/p derivation package)
    }
    trk_n_pass_extrapol++;

    float trk_p = 0.;
    if (fabs((trk->qOverP())>0.)) trk_p = (1./fabs(trk->qOverP()))/1e3; 
    float trk_pt = trk->pt()/1e3;
    float trk_eta = trk->eta();
    float trk_phi = trk->phi();
    float d0 = trk->d0();
    float pvz = 0.;
    if (pv) pvz = pv->z(); 
    float z0 = trk->z0() + trk->vz() - pvz;

    // Track p cut to increase fraction of tracks reaching TileCal
    if (trk_p < 2.0) continue;
    trk_n_pass_p++;

    // Limit |eta| to within TileCal
    if (fabs(trk_eta) > 1.7) continue;
    trk_n_pass_eta++;

    // Track isolation: Require there to be no other tracks within a cone of DR=0.4 around a selected track
    bool trk_not_isolated = false;
    trk2_itr = trks->begin();
    for( ; trk2_itr != trk2_end; ++trk2_itr ) {
      if (trk_itr == trk2_itr) continue; // do not double count the selected track 
      const xAOD::TrackParticle* trk2 = (*trk2_itr);
      if (trk->p4().DeltaR(trk2->p4()) <= m_trkIsoDRmax) trk_not_isolated = true;
    }
    if (trk_not_isolated) continue;
    trk_n_pass_isolation++;

    // These selections require looking at different layers in the calorimeter
    if (m_doTrkExtrapolation) {
      // TODO Need to implement track extrapolation in order for these to work
      // (need to know the fraction of matched clusters in TileCal layers vs. LAr EM layers)
    }
    trk_n_pass_larEmax++;
    trk_n_pass_tileEfrac++;

    // Fill trk hists
    m_trk_pt->Fill(trk_pt);
    m_trk_p->Fill(trk_p);
    m_trk_eta->Fill(trk_eta);
    m_trk_phi->Fill(trk_phi);
    m_trk_d0->Fill(d0);
    m_trk_z0->Fill(z0);
    m_trk_z0sinT->Fill(z0*sin(trk->theta()));
    m_trk_chi2prob->Fill(TMath::Prob(trk->chiSquared(), trk->numberDoF()));
    m_trk_charge->Fill(trk->charge());

    // Get the energy sum of clusters matched to the track (within a cone of DR=0.2)
    float matched_clus_sum_e = getTrkMatchedClusterSumEnergy(trk, clusters, m_trkClusMatchDRmax);
    float matched_cell_sum_e = getTrkMatchedCellSumEnergy(trk, cells, m_trkClusMatchDRmax);
    float clus_eop = matched_clus_sum_e/trk_p;
    float cell_eop= matched_cell_sum_e/trk_p;

    // Fill E/p hists
    m_clus_eop->Fill(clus_eop);
    m_clus_avgeop_vs_trkP->Fill(trk_p, clus_eop);
    m_clus_avgeop_vs_trkEta->Fill(trk_eta, clus_eop);
    m_clus_avgeop_vs_trkPhi->Fill(trk_phi, clus_eop);
    m_cell_eop->Fill(cell_eop);
    m_cell_avgeop_vs_trkP->Fill(trk_p, cell_eop);
    m_cell_avgeop_vs_trkEta->Fill(trk_eta, cell_eop);
    m_cell_avgeop_vs_trkPhi->Fill(trk_phi, cell_eop);

  } // END looping trk

  m_trk_n_pass_trksel->Fill(trk_n_pass_trksel);
  // m_trk_n_pass_extrapol->Fill(trk_n_pass_extrapol);
  m_trk_n_pass_p->Fill(trk_n_pass_p);
  m_trk_n_pass_eta->Fill(trk_n_pass_eta);
  m_trk_n_pass_isolation->Fill(trk_n_pass_isolation);
  // m_trk_n_pass_larEmax->Fill(trk_n_pass_larEmax);
  // m_trk_n_pass_tileEfrac->Fill(trk_n_pass_tileEfrac);

  return StatusCode::SUCCESS;
}

/*---------------------------------------------------------*/
StatusCode TileEoverPMonTool::procHistograms() {
/*---------------------------------------------------------*/

  // ATH_MSG_INFO( "in TileEoverPMonTool::procHistograms()" );

  return StatusCode::SUCCESS;
}

/*---------------------------------------------------------*/
StatusCode TileEoverPMonTool::bookEoverPHistograms() {
/*---------------------------------------------------------*/

  ATH_MSG_INFO( "in TileEoverPMonTool::bookEoverPHistograms()" );

  // track hists
  unsigned int nBinsTrkN = 200; float minTrkN = -0.5;         float maxTrkN = 199.5;
  unsigned int nBinsP = 100;    float minP = 0;               float maxP = 50;
  unsigned int nBinsPhi = 64;   float minPhi = -TMath::Pi();  float maxPhi = TMath::Pi(); 
  unsigned int nBinsEta = 50;   float minEta = -2.5;          float maxEta = 2.5;

  m_trk_n_all = book1F(m_prefix, "trk_n_all", "trk_n_all", nBinsTrkN, minTrkN, maxTrkN);
  m_trk_n_all-> GetXaxis()->SetTitle("N_{trk}");
  m_trk_n_pass_trksel = book1F(m_prefix, "trk_n_pass_trksel", "trk_n_pass_trksel", nBinsTrkN, minTrkN, maxTrkN);
  m_trk_n_pass_trksel-> GetXaxis()->SetTitle("N_{trk}");
  // m_trk_n_pass_extrapol = book1F(m_prefix, "trk_n_pass_extrapol", "trk_n_pass_extrapol", nBinsTrkN, minTrkN, maxTrkN);
  // m_trk_n_pass_extrapol-> GetXaxis()->SetTitle("N_{trk}");
  m_trk_n_pass_p = book1F(m_prefix, "trk_n_pass_p", "trk_n_pass_p", nBinsTrkN, minTrkN, maxTrkN);
  m_trk_n_pass_p-> GetXaxis()->SetTitle("N_{trk}");
  m_trk_n_pass_eta = book1F(m_prefix, "trk_n_pass_eta", "trk_n_pass_eta", nBinsTrkN, minTrkN, maxTrkN);
  m_trk_n_pass_eta-> GetXaxis()->SetTitle("N_{trk}");
  m_trk_n_pass_isolation = book1F(m_prefix, "trk_n_pass_isolation", "trk_n_pass_isolation", nBinsTrkN, minTrkN, maxTrkN);
  m_trk_n_pass_isolation-> GetXaxis()->SetTitle("N_{trk}");
  // m_trk_n_pass_larEmax = book1F(m_prefix, "trk_n_pass_larEmax", "trk_n_pass_larEmax", nBinsTrkN, minTrkN, maxTrkN);
  // m_trk_n_pass_larEmax-> GetXaxis()->SetTitle("N_{trk}");
  // m_trk_n_pass_tileEfrac = book1F(m_prefix, "trk_n_pass_tileEfrac", "trk_n_pass_tileEfrac", nBinsTrkN, minTrkN, maxTrkN); 
  // m_trk_n_pass_tileEfrac-> GetXaxis()->SetTitle("N_{trk}");

  m_trk_pt = book1F(m_prefix, "trk_pt", "trk_pt", nBinsP, minP, maxP);
  m_trk_pt-> GetXaxis()->SetTitle("p_{T, trk} [GeV]");
  m_trk_p = book1F(m_prefix, "trk_p", "trk_p", nBinsP, minP, maxP);
  m_trk_p->GetXaxis()->SetTitle("p_{trk} [GeV]");
  m_trk_eta = book1F(m_prefix, "trk_eta", "trk_eta", nBinsEta, minEta, maxEta);
  m_trk_eta->GetXaxis()->SetTitle("#eta_{trk}");
  m_trk_phi = book1F(m_prefix, "trk_phi", "trk_phi", nBinsPhi, minPhi, maxPhi);
  m_trk_phi->GetXaxis()->SetTitle("#phi_{trk}");
  m_trk_d0 = book1F(m_prefix, "trk_d0", "trk_d0", 100, -5.0, 5.0);
  m_trk_d0->GetXaxis()->SetTitle("d_{0} [mm]");
  m_trk_z0 = book1F(m_prefix, "trk_z0", "trk_z0", 100, -5.0, 5.0);
  m_trk_z0->GetXaxis()->SetTitle("z_{0} [mm]");
  m_trk_z0sinT = book1F(m_prefix, "trk_z0sinT", "trk_z0sinT", 100, -5.0, 5.0);
  m_trk_z0sinT->GetXaxis()->SetTitle("z_{0}sin(#theta) [mm]");
  m_trk_chi2prob = book1F(m_prefix, "trk_chi2prob", "trk_chi2prob", 100, -0.01, 1.0);
  m_trk_chi2prob->GetXaxis()->SetTitle("chi2prob [mm]");
  m_trk_charge = book1F(m_prefix, "trk_charge", "trk_charge", 3, -1.5, 1.5);
  m_trk_charge->GetXaxis()->SetTitle("charge");

  // cluster hists
  unsigned int nBinsClusN = 125;      float minClusN = -0.5;      float maxClusN = 499.5;
  unsigned int nBinsClusMatchN = 20;  float minClusMatchN = -0.5; float maxClusMatchN = 19.5;
  unsigned int nBinsClusE = 100;      float minClusE = 0;         float maxClusE = 100;
  unsigned int nBinsClusE_l = 100;    float minClusE_l = 0;       float maxClusE_l = 500;

  m_clus_n_all = book1F(m_prefix, "clus_n_all", "clus_n_all", nBinsClusN, minClusN, maxClusN);
  m_clus_n_all -> GetXaxis()->SetTitle("N_{cluster}");
  m_clus_e_all = book1F(m_prefix, "clus_e_all", "clus_e_all", nBinsClusE, minClusE, maxClusE);
  m_clus_e_all->GetXaxis()->SetTitle("E_{cluster} [GeV]");
  m_clus_e_sum_all = book1F(m_prefix, "clus_e_sum_all", "clus_e_sum_all", nBinsClusE_l, minClusE_l, maxClusE_l);
  m_clus_e_sum_all->GetXaxis()->SetTitle("E_{cluster} [GeV]");
  m_clus_eta_all = book1F(m_prefix, "clus_eta_all", "clus_eta_all", nBinsEta, minEta, maxEta);
  m_clus_eta_all->GetXaxis()->SetTitle("#eta_{cluster}");
  m_clus_phi_all = book1F(m_prefix, "clus_phi_all", "clus_phi_all", nBinsPhi, minPhi, maxPhi);
  m_clus_phi_all->GetXaxis()->SetTitle("#phi_{cluster}");
  m_clus_n_matched = book1F(m_prefix, "clus_n_matched", "clus_n_matched", nBinsClusMatchN, minClusMatchN, maxClusMatchN);
  m_clus_n_matched -> GetXaxis()->SetTitle("N_{cluster}");
  m_clus_e_matched = book1F(m_prefix, "clus_e_matched", "clus_e_matched", nBinsClusE, minClusE, maxClusE);
  m_clus_e_matched->GetXaxis()->SetTitle("E_{cluster} [GeV]");
  m_clus_e_sum_matched = book1F(m_prefix, "clus_e_sum_matched", "clus_e_sum_matched", nBinsClusE, minClusE, maxClusE);
  m_clus_e_sum_matched->GetXaxis()->SetTitle("E_{cluster} [GeV]");
  m_clus_eta_matched = book1F(m_prefix, "clus_eta_matched", "clus_eta_matched", nBinsEta, minEta, maxEta);
  m_clus_eta_matched->GetXaxis()->SetTitle("#eta_{cluster}");
  m_clus_phi_matched = book1F(m_prefix, "clus_phi_matched", "clus_phi_matched", nBinsPhi, minPhi, maxPhi);
  m_clus_phi_matched->GetXaxis()->SetTitle("#phi_{cluster}");

  // cell hists
  unsigned int nBinsCellMatchN = 125; float minCellMatchN = -0.5; float maxCellMatchN = 999.5;
  m_cell_e_all = book1F(m_prefix, "cell_e_all", "cell_e_all", nBinsClusE, minClusE, maxClusE);
  m_cell_e_all->GetXaxis()->SetTitle("E_{cell} [GeV]");
  m_cell_e_sum_all = book1F(m_prefix, "cell_e_sum_all", "cell_e_sum_all", nBinsClusE_l, minClusE_l, maxClusE_l);
  m_cell_e_sum_all->GetXaxis()->SetTitle("E_{cell} [GeV]");
  m_cell_eta_all = book1F(m_prefix, "cell_eta_all", "cell_eta_all", nBinsEta, minEta, maxEta);
  m_cell_eta_all->GetXaxis()->SetTitle("#eta_{cell}");
  m_cell_phi_all = book1F(m_prefix, "cell_phi_all", "cell_phi_all", nBinsPhi, minPhi, maxPhi);
  m_cell_phi_all->GetXaxis()->SetTitle("#phi_{cell}");
  m_cell_n_matched = book1F(m_prefix, "cell_n_matched", "cell_n_matched", nBinsCellMatchN, minCellMatchN, maxCellMatchN);
  m_cell_n_matched -> GetXaxis()->SetTitle("N_{cell}");
  m_cell_e_matched = book1F(m_prefix, "cell_e_matched", "cell_e_matched", nBinsClusE, minClusE, maxClusE);
  m_cell_e_matched->GetXaxis()->SetTitle("E_{cell} [GeV]");
  m_cell_e_sum_matched = book1F(m_prefix, "cell_e_sum_matched", "cell_e_sum_matched", nBinsClusE, minClusE, maxClusE);
  m_cell_e_sum_matched->GetXaxis()->SetTitle("E_{cell} [GeV]");
  m_cell_eta_matched = book1F(m_prefix, "cell_eta_matched", "cell_eta_matched", nBinsEta, minEta, maxEta);
  m_cell_eta_matched->GetXaxis()->SetTitle("#eta_{cell}");
  m_cell_phi_matched = book1F(m_prefix, "cell_phi_matched", "cell_phi_matched", nBinsPhi, minPhi, maxPhi);
  m_cell_phi_matched->GetXaxis()->SetTitle("#phi_{cell}");

  // E/p hists
  unsigned int nBinsEop = 120;     float minEop = -2;            float maxEop = 10;
  m_clus_eop = book1F(m_prefix, "clus_eop", "clus_eop", nBinsEop, minEop, maxEop);
  m_clus_eop->GetXaxis()->SetTitle("E/p");
  m_clus_avgeop_vs_trkP = bookProfile(m_prefix, "clus_avgeop_vs_trkP", "clus_avgeop_vs_trkP", nBinsP, minP, maxP, minEop, maxEop);
  m_clus_avgeop_vs_trkP->GetXaxis()->SetTitle("p_{T, trk} [GeV]");
  m_clus_avgeop_vs_trkP->GetYaxis()->SetTitle("<E/p>");
  m_clus_avgeop_vs_trkEta = bookProfile(m_prefix, "clus_avgeop_vs_trkEta", "clus_avgeop_vs_trkEta", nBinsEta, minEta, maxEta, minEop, maxEop);
  m_clus_avgeop_vs_trkEta->GetXaxis()->SetTitle("#eta_{trk}");
  m_clus_avgeop_vs_trkEta->GetYaxis()->SetTitle("<E/p>");
  m_clus_avgeop_vs_trkPhi = bookProfile(m_prefix, "clus_avgeop_vs_trkPhi", "clus_avgeop_vs_trkPhi", nBinsPhi, minPhi, maxPhi, minEop, maxEop);
  m_clus_avgeop_vs_trkPhi->GetXaxis()->SetTitle("#phi_{trk}");
  m_clus_avgeop_vs_trkPhi->GetYaxis()->SetTitle("<E/p>");

  m_cell_eop = book1F(m_prefix, "cell_eop", "cell_eop", nBinsEop, minEop, maxEop);
  m_cell_eop->GetXaxis()->SetTitle("E/p");
  m_cell_avgeop_vs_trkP = bookProfile(m_prefix, "cell_avgeop_vs_trkP", "cell_avgeop_vs_trkP", nBinsP, minP, maxP, minEop, maxEop);
  m_cell_avgeop_vs_trkP->GetXaxis()->SetTitle("p_{T, trk} [GeV]");
  m_cell_avgeop_vs_trkP->GetYaxis()->SetTitle("<E/p>");
  m_cell_avgeop_vs_trkEta = bookProfile(m_prefix, "cell_avgeop_vs_trkEta", "cell_avgeop_vs_trkEta", nBinsEta, minEta, maxEta, minEop, maxEop);
  m_cell_avgeop_vs_trkEta->GetXaxis()->SetTitle("#eta_{trk}");
  m_cell_avgeop_vs_trkEta->GetYaxis()->SetTitle("<E/p>");
  m_cell_avgeop_vs_trkPhi = bookProfile(m_prefix, "cell_avgeop_vs_trkPhi", "cell_avgeop_vs_trkPhi", nBinsPhi, minPhi, maxPhi, minEop, maxEop);
  m_cell_avgeop_vs_trkPhi->GetXaxis()->SetTitle("#phi_{trk}");
  m_cell_avgeop_vs_trkPhi->GetYaxis()->SetTitle("<E/p>");

  return StatusCode::SUCCESS;
}

/*---------------------------------------------------------*/
float TileEoverPMonTool::getTrkMatchedClusterSumEnergy(const xAOD::TrackParticle* trk, const xAOD::CaloClusterContainer* clusters, float radius) {
/*---------------------------------------------------------*/

  m_clus_n_all->Fill(clusters->size());
  int clus_n_matched = 0;
  float clus_e_sum_all = 0.;
  float clus_e_sum_matched = 0.;

  for (const auto& clus : *clusters) {

    float clus_e = clus->e()/1e3;
    float clus_eta = clus->eta();
    float clus_phi = clus->phi();

    clus_e_sum_all += clus_e;
    m_clus_e_all->Fill(clus_e);
    m_clus_eta_all->Fill(clus_eta);
    m_clus_phi_all->Fill(clus_phi);

    if (trk->p4().DeltaR(clus->p4()) <= radius) {
      clus_n_matched++;
      clus_e_sum_matched += clus_e;
      m_clus_e_matched->Fill(clus_e);
      m_clus_eta_matched->Fill(clus_eta);
      m_clus_phi_matched->Fill(clus_phi);
    }
  }

  m_clus_n_matched->Fill(clus_n_matched);
  m_clus_e_sum_all->Fill(clus_e_sum_all);
  m_clus_e_sum_matched->Fill(clus_e_sum_matched);

  return clus_e_sum_matched;
}

/*---------------------------------------------------------*/
float TileEoverPMonTool::getTrkMatchedCellSumEnergy(const xAOD::TrackParticle* trk, const CaloCellContainer* cells, float radius) {
/*---------------------------------------------------------*/

  int cell_n_matched = 0;
  float cell_e_sum_all = 0.;
  float cell_e_sum_matched = 0.;

  for (const auto& cell : *cells) {

    float cell_e = cell->e()/1e3;
    float cell_eta = cell->eta();
    float cell_phi = cell->phi();

    cell_e_sum_all += cell_e;
    m_cell_e_all->Fill(cell_e);
    m_cell_eta_all->Fill(cell_eta);
    m_cell_phi_all->Fill(cell_phi);

    float deta = fabs(trk->eta() - cell_eta);
    float dphi = fabs(trk->phi() - cell_phi);
    if (dphi > TMath::Pi()) dphi = 2*TMath::Pi() - dphi;
    float dR = sqrt( pow(deta,2) + pow(dphi,2) );
    if (dR <= radius) {
      cell_n_matched++;
      cell_e_sum_matched += cell_e;
      m_cell_e_matched->Fill(cell_e);
      m_cell_eta_matched->Fill(cell_eta);
      m_cell_phi_matched->Fill(cell_phi);
    }
  }

  m_cell_n_matched->Fill(cell_n_matched);
  m_cell_e_sum_all->Fill(cell_e_sum_all);
  m_cell_e_sum_matched->Fill(cell_e_sum_matched);

  return cell_e_sum_matched;
}
