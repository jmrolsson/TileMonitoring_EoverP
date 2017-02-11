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
  
  private:

    bool m_first_event;

    std::string m_prefix;
    std::string m_eventInfoContainerName;
    std::string m_vxContainerName;
    std::string m_trackContainerName;
    std::string m_caloClusterContainerName;
    std::string m_cellContainerName;

    ToolHandle<Trk::IExtrapolator> m_extrapolator;
    ToolHandle<Trk::IParticleCaloExtensionTool> m_theTrackExtrapolatorTool;
    Trk::TrackParametersIdHelper* m_trackParametersIdHelper;
    ToolHandle<ICaloSurfaceHelper> m_surfaceHelper;
    const TileTBID* m_tileTBID; 

    //// declare histograms

    // track kinematics
    TH1F* m_trk_pt;
    TH1F* m_trk_p;
    TH1F* m_trk_eta;
    TH1F* m_trk_phi;

};

#endif
