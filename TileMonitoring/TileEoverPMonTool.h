// ********************************************************************
//
// NAME:     TileEoverPMonTool.h
// PACKAGE:  TileMonitoring
//
// AUTHOR:   Joakim Olsson (joakim.olsson@cern.ch)
//
// ********************************************************************
#ifndef TILEMONITORING_TILEEOVERPMONTOOL_H
#define TILEMONITORING_TILEEOVERPMONTOOL_H

#include "TileMonitoring/TileFatherMonTool.h"

#include "AsgTools/ToolHandle.h"
#include "CaloIdentifier/CaloCell_ID.h"
#include "RecoToolInterfaces/IParticleCaloExtensionTool.h"
#include "xAODCaloEvent/CaloClusterContainer.h"
#include "xAODCaloEvent/CaloClusterChangeSignalState.h"
#include "CaloEvent/CaloCellContainer.h"

#include "TH1F.h"
#include "TH2F.h"
#include "TProfile2D.h"

#include <string>
#include <vector>

class TileTBID;
class ICaloSurfaceHelper;

namespace Trk {
  class IExtrapolator;
  class Surface;
  class TrackParametersIdHelper;
}

namespace InDet {
  class IInDetTrackSelectionTool;
}

/** @class TileEoverPMonTool
 *  @brief Class for E/p monitoring histograms and tools
 */

class TileEoverPMonTool: public TileFatherMonTool {

  public:
    TileEoverPMonTool(const std::string & type, const std::string & name, const IInterface* parent);

    ~TileEoverPMonTool();

    virtual StatusCode initialize();
    virtual StatusCode fillHistograms();
    virtual StatusCode bookHistograms();
    virtual StatusCode procHistograms();
    StatusCode bookEoverPHistograms();

    float getTrkMatchedClusterSumEnergy(const xAOD::TrackParticle* trk, const xAOD::CaloClusterContainer* clusters, float radius = 0.2);
    float getTrkMatchedCellSumEnergy(const xAOD::TrackParticle* trk, const CaloCellContainer* cells, float radius = 0.2);

  private:

    bool m_first_event;

    bool m_doTrkExtrapolation;
    bool m_doBasicTrkSelection;

    float m_trkIsoDRmax;
    float m_trkClusMatchDRmax;

    std::string m_prefix;
    std::string m_eventInfoContainerName;
    std::string m_vertexContainerName;
    std::string m_trackContainerName;
    std::string m_caloClusterContainerName;
    std::string m_cellContainerName;

    ToolHandle<Trk::IExtrapolator> m_extrapolator;
    ToolHandle<Trk::IParticleCaloExtensionTool> m_trackExtrapolatorTool;
    Trk::TrackParametersIdHelper* m_trackParametersIdHelper;
    ToolHandle<ICaloSurfaceHelper> m_surfaceHelper;
    const TileTBID* m_tileTBID; 

    // track selection tool
    ToolHandle< InDet::IInDetTrackSelectionTool > m_trackSelectionTool;

    //// declare histograms

    // track hists 
    TH1F* m_trk_n_all;
    TH1F* m_trk_n_pass_trksel;
    TH1F* m_trk_n_pass_extrapol;
    TH1F* m_trk_n_pass_p;
    TH1F* m_trk_n_pass_eta;
    TH1F* m_trk_n_pass_isolation;
    TH1F* m_trk_n_pass_larEmax;
    TH1F* m_trk_n_pass_tileEfrac;
    TH1F* m_trk_pt;
    TH1F* m_trk_p;
    TH1F* m_trk_eta;
    TH1F* m_trk_phi;
    TH1F* m_trk_d0;
    TH1F* m_trk_z0;
    TH1F* m_trk_z0sinT;
    TH1F* m_trk_chi2prob;
    TH1F* m_trk_charge;

    // cluster hists
    TH1F* m_clus_n_all;
    TH1F* m_clus_e_all;
    TH1F* m_clus_e_sum_all;
    TH1F* m_clus_eta_all;
    TH1F* m_clus_phi_all;
    TH1F* m_clus_n_matched;
    TH1F* m_clus_e_matched;
    TH1F* m_clus_e_sum_matched;
    TH1F* m_clus_eta_matched;
    TH1F* m_clus_phi_matched;

    // cell hists
    TH1F* m_cell_e_all;
    TH1F* m_cell_e_sum_all;
    TH1F* m_cell_eta_all;
    TH1F* m_cell_phi_all;
    TH1F* m_cell_n_matched;
    TH1F* m_cell_e_matched;
    TH1F* m_cell_e_sum_matched;
    TH1F* m_cell_eta_matched;
    TH1F* m_cell_phi_matched;

    // E/p hists
    TH1F* m_clus_eop;
    TProfile* m_clus_avgeop_vs_trkP;
    TProfile* m_clus_avgeop_vs_trkEta;
    TProfile* m_clus_avgeop_vs_trkPhi;

    TH1F* m_cell_eop;
    TProfile* m_cell_avgeop_vs_trkP;
    TProfile* m_cell_avgeop_vs_trkEta;
    TProfile* m_cell_avgeop_vs_trkPhi;

};

#endif
