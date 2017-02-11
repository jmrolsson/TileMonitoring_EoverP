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

#include <string>
#include <vector>
#include <cmath>
#include <utility>
#include <limits>

/*---------------------------------------------------------*/
TileEoverPMonTool::TileEoverPMonTool(const std::string & type, const std::string & name, const IInterface* parent) : 
  TileFatherMonTool(type, name, parent),
  m_prefix(""),
  m_eventInfoContainerName("EventInfo"),
  m_vxContainerName("PrimaryVertices"),
  m_trackContainerName("InDetTrackParticles"),
  m_caloClusterContainerName("CaloCalTopoClusters"),
  m_cellContainerName("AllCalo"),
  m_extrapolator("Trk::Extrapolator"),
  m_theTrackExtrapolatorTool("Trk::ParticleCaloExtensionTool"),
  m_trackParametersIdHelper(new Trk::TrackParametersIdHelper),
  m_surfaceHelper("CaloSurfaceHelper/CaloSurfaceHelper"),
  m_tileTBID(0) {
/*---------------------------------------------------------*/
  declareInterface<IMonitorToolBase>(this);
  declareProperty("Prefix", m_prefix);
  declareProperty("EventContainer", m_eventInfoContainerName);
  declareProperty("VertexContainer", m_vxContainerName);
  declareProperty("TrackContainer", m_trackContainerName);
  declareProperty("CaloClusterContainer", m_caloClusterContainerName);
  declareProperty("CellContainer", m_cellContainerName);
  declareProperty("Extrapolator", m_extrapolator);
  declareProperty("TheTrackExtrapolatorTool", m_theTrackExtrapolatorTool);

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

  ATH_CHECK( m_extrapolator.retrieve() );
  ATH_CHECK( m_theTrackExtrapolatorTool.retrieve() );
  ATH_CHECK( m_surfaceHelper.retrieve() );
  // Get the test beam identifier for the MBTS
  ATH_CHECK( detStore()->retrieve(m_tileTBID) );

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

  ATH_MSG_INFO( "in TileEoverPMonTool::fillHistograms()" );

  fillEvtInfo();
  if (m_first_event) {
    m_first_event = false;
    if (bookEoverPHistograms() != StatusCode::SUCCESS) {
      ATH_MSG_ERROR("Something wrong when booking histograms.");
    }
  }	

  const xAOD::TrackParticleContainer* trks = 0;
  CHECK(evtStore()->retrieve(trks, m_trackContainerName));

  xAOD::TrackParticleContainer::const_iterator trk_itr = trks->begin();
  xAOD::TrackParticleContainer::const_iterator trk_end = trks->end();
  for( ; trk_itr != trk_end; ++trk_itr ) {

    const xAOD::TrackParticle* trk = (*trk_itr);
    float trk_p = 0;
    if (TMath::Abs(trk->qOverP())>0.) trk_p = (1./TMath::Abs(trk->qOverP()))/1e3; 
    float trk_pt = trk->pt();
    float trk_eta = trk->eta();
    float trk_phi = trk->phi();

    m_trk_pt->Fill(trk_pt);
    m_trk_p->Fill(trk_p);
    m_trk_eta->Fill(trk_eta);
    m_trk_phi->Fill(trk_phi);

  }
 
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

  unsigned int nBinsP = 100;   float minP = 0;               float maxP = 50;
  unsigned int nBinsPhi = 64;  float minPhi = -TMath::Pi();  float maxPhi = TMath::Pi(); 
  unsigned int nBinsEta = 50;  float minEta = -2.5;          float maxEta = 2.5;

  m_trk_pt = book1F(m_prefix, "trk_pt", "trk_pt", nBinsP, minP, maxP);
  m_trk_pt-> GetXaxis()->SetTitle("p_{T, trk} [GeV]");
  m_trk_p = book1F(m_prefix, "trk_p", "trk_p", nBinsP, minP, maxP);
  m_trk_p->GetXaxis()->SetTitle("p_{trk} [GeV]");
  m_trk_eta = book1F(m_prefix, "trk_eta", "trk_eta", nBinsEta, minEta, maxEta);
  m_trk_eta->GetXaxis()->SetTitle("#eta_{trk}");
  m_trk_phi = book1F(m_prefix, "trk_phi", "trk_phi", nBinsPhi, minPhi, maxPhi);
  m_trk_phi->GetXaxis()->SetTitle("#phi_{trk}");

  return StatusCode::SUCCESS;
}
