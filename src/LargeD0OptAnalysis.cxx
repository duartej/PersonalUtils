////////////////////////////////////////////////////////////////////////////////// 
/// Extract some information related the Space Points and reconstructed tracks
/// to compare with different (tracking) algorithm options
///
/// Author: Jordi Duarte-Campderrros
/// Tel Aviv University, Jan. 14th, 2016
///
/// DESCRIPTION:
///
///
//////////////////////////////////////////////////////////////////////////////////

#include "PersonalUtils/LargeD0OptAnalysis.h"
// the first two come for free when using AthAlgorithm
//#include "GaudiKernel/AlgFactory.h"
//#include "GaudiKernel/IToolSvc.h"

// Gaudi
#include "GaudiKernel/ServiceHandle.h"
#include "GaudiKernel/ITHistSvc.h"

// Athena Specific
//#include "xAODEventInfo/EventInfo.h"

#include "TrkSpacePoint/SpacePointContainer.h"
#include "SiSpacePoint/SCT_SpacePoint.h"

//#include "TrkTrack/TrackCollection.h"

// ROOT
#include "TTree.h"

// System
#include <math.h>
#include <functional>
#include <iostream>

//////////////////////////////////////////////////////////////////////////////////////
/// Constructor

LargeD0OptAnalysis::LargeD0OptAnalysis(const std::string& name, ISvcLocator* pSvcLocator) :
    AthAlgorithm(name, pSvcLocator),
    m_processedEvts(0),
    m_pixelSPName("PixelSpacePoints"),
    m_SCTSPName("SCT_SpacePoints"),
    m_OverlapSPName("OverlapSpacePoints"),
    m_SiSeededTrackName("SiSpSeededLargeD0Tracks"),
    m_ResolvedTrackName("ResolvedLargeD0Tracks"),
    m_ExtendedTrackName("ExtendedLargeD0Tracks"),
    m_streamHist("SP_OPT"),
    m_tree(nullptr)
{
    /** switches to control the analysis through job options */
    declareProperty( "PixelSpacePoints",   m_pixelSPName   );
    declareProperty( "SCTSpacePoints",     m_SCTSPName     );
    declareProperty( "OverlapSpacePoints", m_OverlapSPName );

    declareProperty( "SiSeededTracks",     m_SiSeededTrackName );
    declareProperty( "ResolvedTracks",     m_ResolvedTrackName );
    declareProperty( "Extendedracks",     m_ExtendedTrackName );
}

/////////////////////////////////////////////////////////////////////////////////////
/// Destructor - check up memory allocation
/// delete any memory allocation on the heap

LargeD0OptAnalysis::~LargeD0OptAnalysis() {}

////////////////////////////////////////////////////////////////////////////////////
StatusCode LargeD0OptAnalysis::beginRun() 
{ 
    return StatusCode::SUCCESS;
} 

StatusCode LargeD0OptAnalysis::initialize() 
{
    ATH_MSG_DEBUG("Initializing LargeD0OptAnalysis");

    // The histogram service
    ServiceHandle<ITHistSvc> tHistSvc("THistSvc",this->name());
    StatusCode sc = tHistSvc.retrieve();
    if( sc.isFailure() )
    {
        ATH_MSG_FATAL( "Unable to retrieve pointer to THistSvc" );
        return sc;
    }

    // Initialization and registration of the histograms/tree
    m_tree = new TTree("SpacePoint_STEP","SpacePoint finder step" );
    CHECK( tHistSvc->regTree(std::string("/"+m_streamHist+"/SPacePoints_STEP").c_str(),m_tree) );

    m_tree->Branch("sct_SpacePoints",&m_lightSP,"x/F:y/F:z/F:r/F:phi/F");

    return StatusCode::SUCCESS;
}		 

///////////////////////////////////////////////////////////////////////////////////
/// Finalize - delete any memory allocation from the heap

StatusCode LargeD0OptAnalysis::finalize() 
{
    return StatusCode::SUCCESS;
}

//////////////////////////////////////////////////////////////////////////////////
/// Execute - on event by event

StatusCode LargeD0OptAnalysis::execute() 
{
    ++m_processedEvts;
    ATH_MSG_DEBUG(" in execute()");

    StatusCode sc;
    
    const SpacePointContainer *  sct_sp_container = nullptr;
    sc = evtStore()->retrieve( sct_sp_container, m_SCTSPName);
    if( sc.isFailure()  ||  ! sct_sp_container ) 
    {
        ATH_MSG_FATAL("No SpacePointContainer found with name '" 
                << m_SCTSPName << " !");
        return StatusCode::FAILURE;
    }

    for(auto & sct_sp_collection : *sct_sp_container )
    {
        for(auto & sct_sp : *sct_sp_collection)
        {
            m_lightSP.x = sct_sp->globalPosition().x();
            m_lightSP.y = sct_sp->globalPosition().y();
            m_lightSP.z = sct_sp->globalPosition().z();
            m_lightSP.r = sct_sp->r();
            m_lightSP.phi = sct_sp->phi();

            m_tree->Fill();
        }
    }
  
    ATH_MSG_DEBUG("End processed event " << m_processedEvts);
  
    return StatusCode::SUCCESS;
}

