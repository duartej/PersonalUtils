////////////////////////////////////////////////////////////////////////////////// 
/// Just extract some info from PRD multitruth
///
/// Author: Jordi Duarte-Campderrros
/// Tel Aviv University, Setember 13th, 2015
///
/// DESCRIPTION: Extract some info from the simulated hits per event
///
///
//////////////////////////////////////////////////////////////////////////////////

#include "PersonalUtils/ExtractTruthPRDInfo.h"
// the first two come for free when using AthAlgorithm
//#include "GaudiKernel/AlgFactory.h"
//#include "GaudiKernel/IToolSvc.h"

// Gaudi
#include "GaudiKernel/ServiceHandle.h"
#include "GaudiKernel/ITHistSvc.h"

// Event
#include "EventInfo/EventInfo.h"
#include "EventInfo/EventID.h"

// HepMC
#include "HepMC/GenVertex.h"
#include "HepMC/GenParticle.h"
#include "GeneratorObjects/McEventCollection.h"
//#include "HepPDT/ParticleDataTable.hh"

// Athena
#include "InDetSimEvent/SiHitCollection.h"
#include "InDetSimEvent/TRTUncompressedHitCollection.h"

#include "InDetPrepRawData/PixelClusterContainer.h"
#include "InDetPrepRawData/SCT_ClusterContainer.h"
#include "InDetPrepRawData/TRT_DriftCircleContainer.h"

#include "xAODEventInfo/EventInfo.h"

//#include "TrkTruthData/PRD_MultiTruthCollection.h"
//

// ROOT 
#include "TTree.h"

// std
#include<map>
#include<math.h>
#include<algorithm>
#include<string>

//////////////////////////////////////////////////////////////////////////////////////
/// Constructor

ExtractTruthPRDInfo::ExtractTruthPRDInfo(const std::string& name, ISvcLocator* pSvcLocator) :
    AthAlgorithm(name, pSvcLocator),
    m_processedEvts(0),
    m_recoMode(true),
    m_prdTruthPixel("PRD_MultiTruthPixel"),
    m_prdTruthSCT("PRD_MultiTruthSCT"),
    m_prdTruthTRT("PRD_MultiTruthTRT"),
    m_prdPixel("PixelClusters"),
    m_prdSCT("SCT_Clusters"),
    m_prdTRT("TRT_DriftCircles"),
    m_simHitsPixel("PixelHits"),
    m_simHitsSCT("SCT_Hits"),
    m_simHitsTRT("TRTUncompressedHits"),
    m_streamHist("SIMHITS"),
    m_tree(nullptr),
    m_hits({0,0,0,0,0,0,0,0,0,0,0,0})
{
    /** switches to control the analysis through job options */
    // Note that if the reco mode is set to false, it is uses the SimHit collections
    // and therefore a HITS/RDO input file is needed
    declareProperty( "IsRecoMode",   m_recoMode );
    declareProperty( "PRDTruthNamePixel",   m_prdTruthPixel );
    declareProperty( "PRDTruthNameSCT",   m_prdTruthSCT );
    declareProperty( "PRDTruthNameTRT",   m_prdTruthTRT );
    declareProperty( "SimHitNamePixel", m_simHitsPixel );
    declareProperty( "SimHitNameSCT",   m_simHitsSCT );
    declareProperty( "SimHitNameTRT",   m_simHitsTRT );
}

/////////////////////////////////////////////////////////////////////////////////////
/// Destructor - check up memory allocation
/// delete any memory allocation on the heap

ExtractTruthPRDInfo::~ExtractTruthPRDInfo() {}

////////////////////////////////////////////////////////////////////////////////////
StatusCode ExtractTruthPRDInfo::beginRun() 
{ 
    return StatusCode::SUCCESS;
} 

StatusCode ExtractTruthPRDInfo::initialize() 
{
    ATH_MSG_DEBUG("Initializing ExtractTruthPRDInfo");

    // Histogram service
    ServiceHandle<ITHistSvc> tHistSvc("THistSvc",this->name());
    StatusCode sc = tHistSvc.retrieve();
    if( sc.isFailure() )
    {
        ATH_MSG_FATAL( "Unable to retrieve pointer to THistSvc" );
        return sc;
    }
    // Initialization and registration of the histograms/tree
    m_tree = new TTree("TruthHits","SimHit/PRD truth" );
    CHECK( tHistSvc->regTree(std::string("/"+m_streamHist+"/TruthHits").c_str(),m_tree) );
    
    std::string branchName("simHits");
    if(m_recoMode)
    {
        branchName="prdTruth";
    }
    m_tree->Branch(branchName.c_str(),&m_hits,"eventNumber/I:runNumber/I:"\
            "linked_pixel/I:total_pixel/I:linked_sct/I:total_sct/I:"\
            "linked_silicon/I:total_silicon/I:linked_trt/I:total_trt/I:"\
            "linked_all/I:total_all/I");

    return StatusCode::SUCCESS;
}		 

///////////////////////////////////////////////////////////////////////////////////
/// Finalize - delete any memory allocation from the heap
StatusCode ExtractTruthPRDInfo::finalize() 
{
    return StatusCode::SUCCESS;
}

//////////////////////////////////////////////////////////////////////////////////
/// Execute - on event by event
StatusCode ExtractTruthPRDInfo::execute() 
{
    ++m_processedEvts;
    ATH_MSG_DEBUG(" in execute()");
    
    
    // Event info related and hit processing
    if( m_recoMode )
    {
        const xAOD::EventInfo * evtInfo(0);
        ATH_CHECK( evtStore()->retrieve(evtInfo, "EventInfo") );
        const int eN = evtInfo->eventNumber();
        const int rN = evtInfo->runNumber();
        
        ATH_CHECK( processRecoInput(eN,rN) );
    }
    else
    {
        const EventInfo * evtInfo(0);
        ATH_CHECK( evtStore()->retrieve(evtInfo, "McEventInfo") );
        const int eN = evtInfo->event_ID()->event_number();
        const int rN = evtInfo->event_ID()->run_number();
        ATH_CHECK( processHITSInput(eN,rN) );
    }

    // Event related stuff
    const McEventCollection* mcColl(0);
    ATH_CHECK( evtStore()->retrieve(mcColl, "TruthEvent") ); //m_mcCollName) );
    for(auto & evtItr : *mcColl)
    {
        if( evtItr->event_number() == -1 )
        {
            ATH_MSG_INFO("============== SEPARATOR EVENT =============================");
            continue;
        }
        //m_current_EvtNumber = evtItr->event_number();

        ATH_MSG_INFO("+ MC Event Number: " << evtItr->event_number() 
                << " Signal process ID: " << evtItr->signal_process_id()
                << " Number of particles: " << evtItr->particles_size()
                << " Number of vertices: " << evtItr->vertices_size());

        for(HepMC::GenEvent::particle_const_iterator it_p = evtItr->particles_begin();
                              it_p != evtItr->particles_end(); ++it_p) 
        {
            // Find the LLP
            if( (*it_p) && abs((*it_p)->pdg_id()) == 999 )
            {
                ATH_MSG_INFO("GEANTINO FOUND!!");
            }
        }
    }
  
    ATH_MSG_DEBUG("End processed event " << m_processedEvts);

  return StatusCode::SUCCESS;
}

StatusCode ExtractTruthPRDInfo::processRecoInput(const int & evtNumber, const int & runNumber)
{
    // Clusters PRD nd Truth hits
    // ---  PIXEL 
    // - truth
    const PRD_MultiTruthCollection * prdTruthColPixel = nullptr;
    ATH_CHECK( evtStore()->retrieve(prdTruthColPixel, m_prdTruthPixel) );
    // - reco 
    const InDet::PixelClusterContainer * prdColPixel = nullptr;
    ATH_CHECK( evtStore()->retrieve(prdColPixel, m_prdPixel) );

    // --- SCT
    // - truth
    const PRD_MultiTruthCollection * prdTruthColSCT = nullptr;
    ATH_CHECK( evtStore()->retrieve(prdTruthColSCT, m_prdTruthSCT) );
    // - reco
    const InDet::SCT_ClusterContainer * prdColSCT = nullptr;
    ATH_CHECK( evtStore()->retrieve(prdColSCT, m_prdSCT) );
    
    // --- TRT
    // - truth
    const InDet::TRT_DriftCircleContainer * prdColTRT = nullptr;
    ATH_CHECK( evtStore()->retrieve(prdColTRT, m_prdTRT) );
    // - reco
    const PRD_MultiTruthCollection * prdTruthColTRT = nullptr;
    ATH_CHECK( evtStore()->retrieve(prdTruthColTRT, m_prdTruthTRT) );

    // The number of link-no linked
    auto pixHits = getLinkedAndTotalPRDs(prdColPixel,prdTruthColPixel);
    auto sctHits = getLinkedAndTotalPRDs(prdColSCT,prdTruthColSCT);
    auto trtHits = getLinkedAndTotalPRDs(prdColTRT,prdTruthColTRT);

    this->fillHitsTree(evtNumber,runNumber,pixHits,sctHits,trtHits);
    
    ATH_MSG_DEBUG( this->getSummary() );
    
    return StatusCode::SUCCESS;
}

const std::pair<int,int> ExtractTruthPRDInfo::getLinkedAndTotalPRDs(const PRD_MultiTruthCollection * const prdCol)
{
    int nPRDs=0;
    int nPRDsLinked=0;
    for(auto & cluster_genpart: *prdCol)
    {
        ++nPRDs;
        if( cluster_genpart.second.isValid() )
        {
            ATH_MSG_VERBOSE("PRD-ID:" << cluster_genpart.first << " PARTICLE: " 
                    << cluster_genpart.second.cptr()->pdg_id() << " [barcode:" 
                    << cluster_genpart.second.cptr()->barcode() << "]" );
            ++nPRDsLinked;
        }
        else
        {
            ATH_MSG_VERBOSE("PRD-ID:" << cluster_genpart.first << " NOT FOUND MC-GEN LINK, noise?");
        }
    }
    return std::pair<int,int>(nPRDsLinked,nPRDs);
}

template<class PRD_Container> 
const std::pair<int,int> ExtractTruthPRDInfo::getLinkedAndTotalPRDs(const PRD_Container * const prdContainer,
        const PRD_MultiTruthCollection * const prdTruthCol)
{
    std::map<Identifier,int> noiseMap;
    int nPRDs=0;
    int nPRDsLinked=0;
    for(auto & prdCollection: *prdContainer)
    {
        for(auto & prd: *prdCollection)
        {
            ++nPRDs;
            const auto truthGen_it =  prdTruthCol->find(prd->identify());
            // Is the PRD associated to any truth prd?
            if( truthGen_it == prdTruthCol->end() )
            {
                ATH_MSG_VERBOSE("PRD-ID:" << truthGen_it->first << " NOT FOUND MC-GEN PRD LINK, noise?");
                // -- Getting the number of hits associated to the same PRD-ID, probably
                //    it is meaning that the hit was produced with a pile-up event (filtering)
                if( noiseMap.find(truthGen_it->first) != noiseMap.end() )
                {
                    ++noiseMap[truthGen_it->first];
                }
                else
                {
                    noiseMap.emplace(truthGen_it->first,1);
                }
                continue;
            }
            ++nPRDsLinked;
            ATH_MSG_VERBOSE("PRD-ID:" << truthGen_it->first << " PARTICLE: " 
                    << truthGen_it->second.cptr()->pdg_id() << " [barcode:" 
                    << truthGen_it->second.cptr()->barcode() << "]" );
        }
    }
    if( msg().level() <= MSG::DEBUG )
    {
        std::string messageNoise(std::to_string(noiseMap.size())+" NOISE HITS, IDs ==> ");
        for(auto & idAcc_it: noiseMap)
        {
            messageNoise += " "+idAcc_it.first.getString()+" ["+std::to_string(idAcc_it.second)+"]";
        }
        ATH_MSG_DEBUG(messageNoise);
    }

    return std::pair<int,int>(nPRDsLinked,nPRDs);
}

StatusCode ExtractTruthPRDInfo::processHITSInput(const int & evtNumber,const int & runNumber)
{
    std::map<int,bool> pixMap;
    int linkPixel = 0;
    std::map<int,bool> sctMap;
    int linkSCT = 0;
    std::map<int,bool> trtMap;
    int linkTRT = 0;
    // --- Sim hits: Pixel
    const SiHitCollection * simHitColPixel = nullptr;
    ATH_CHECK( evtStore()->retrieve(simHitColPixel, m_simHitsPixel) );
    ATH_MSG_VERBOSE("PIXEL -------------------------");
    for(auto & simhit: * simHitColPixel)
    {
        if( pixMap.emplace(simhit.identify(),simhit.particleLink().isValid()).second )
        {
            if( simhit.particleLink().isValid() )
            {
                ++linkPixel;
            }
            ATH_MSG_VERBOSE("SIM-ID:" << simhit.identify() 
                << " Has_Valid_ParticleLink: " << simhit.particleLink().isValid());
        }
    }
    // --- Sim hits: SCT
    const SiHitCollection * simHitColSCT = nullptr;
    ATH_CHECK( evtStore()->retrieve(simHitColSCT, m_simHitsSCT) );
    ATH_MSG_VERBOSE("SCT -------------------------");
    for(auto & simhit: * simHitColSCT)
    {
        if( sctMap.emplace(simhit.identify(),simhit.particleLink().isValid()).second )
        {
            if( simhit.particleLink().isValid() )
            {
                ++linkSCT;
            }
            ATH_MSG_VERBOSE("SIM-ID:" << simhit.identify() 
                    << " Has_Valid_ParticleLink: " << simhit.particleLink().isValid());
        }
    }
    // --- Sim hits: TRT
    const TRTUncompressedHitCollection * simHitColTRT = nullptr;
    ATH_CHECK( evtStore()->retrieve(simHitColTRT, m_simHitsTRT) );
    ATH_MSG_VERBOSE("TRT -------------------------");
    for(auto & simhit: * simHitColTRT)
    {
        if( trtMap.emplace(simhit.GetHitID(),simhit.particleLink().isValid()).second )
        {
            if( simhit.particleLink().isValid() )
            {
                ++linkTRT;
            }
            ATH_MSG_VERBOSE("SIM-ID:" << simhit.GetHitID() 
                    << " Has_Valid_ParticleLink: " << simhit.particleLink().isValid());
        } 
    }
    // --- FIXME:: TO BE DEPRECATED
    //ATH_MSG_INFO(" ------------------: " <<  linkPixel << " " << linkSCT << " (" 
    //        << (linkPixel+linkSCT) << ") " << linkTRT << " : " <<(linkPixel+linkSCT+linkTRT) );
    //ATH_MSG_INFO(" ------------------: " << pixMap.size() << " " << sctMap.size() << " (" 
    //        << (pixMap.size()+sctMap.size()) << ") " << trtMap.size() << " : " <<(pixMap.size()+sctMap.size()+trtMap.size()) );

    this->fillHitsTree(evtNumber,runNumber,
            linkPixel,pixMap.size(),linkSCT,sctMap.size(),linkTRT,trtMap.size());

    ATH_MSG_DEBUG( this->getSummary() );

    return StatusCode::SUCCESS;
}

const std::string ExtractTruthPRDInfo::getSummary() const
{
    std::ostringstream oss;
    // Assuming is called when the event is already filled
    std::string message( "Event: "+std::to_string(m_hits.eventNumber)+" RunNumber: "+
            std::to_string(m_hits.runNumber)+"\n" );
    message += "---------------------------------------------------------------\n";
    if(m_hits.total_pixel != 0 )
    {
        message += "Pixel: "+std::to_string(m_hits.linked_pixel)+"/"+std::to_string(m_hits.total_pixel)+" ["+
            std::to_string((float(m_hits.linked_pixel)/float(m_hits.total_pixel))*100.0)+"%]\n";
    }
    if(m_hits.total_sct != 0 )
    {
        message += "SCT: "+std::to_string(m_hits.linked_sct)+"/"+std::to_string(m_hits.total_sct)+" ["+
            std::to_string((float(m_hits.linked_sct)/float(m_hits.total_sct))*100.)+"%]\n";
    }
    if(m_hits.total_silicon != 0 )
    {
        message += "++ TOTAL SILICON: "+std::to_string(m_hits.linked_silicon)+"/"+std::to_string(m_hits.total_silicon)+" ["+
            std::to_string((float(m_hits.linked_silicon)/float(m_hits.total_silicon))*100.0)+"%] ++\n";
    }
    if(m_hits.total_trt != 0 )
    {
        message += "TRT: "+std::to_string(m_hits.linked_trt)+"/"+std::to_string(m_hits.total_trt)+" ["+
            std::to_string((float(m_hits.linked_trt)/float(m_hits.total_trt))*100.0)+"%]\n";
    }
    if(m_hits.total_all != 0 )
    {
        message += "+++++ ALL DETECTORS: "+std::to_string(m_hits.linked_all)+"/"+std::to_string(m_hits.total_all)+" ["+
            std::to_string((float(m_hits.linked_all)/float(m_hits.total_all))*100.0)+"%]+++++\n";
    }

    return message;    
}

void ExtractTruthPRDInfo::fillHitsTree(const int & evtNumber, const int & runNumber,
        const int & linked_pixel,const int & total_pixel,
        const int & linked_sct,const int & total_sct,
        const int & linked_trt,const int & total_trt)
{
    m_hits.eventNumber = evtNumber;
    m_hits.runNumber   = runNumber;
    m_hits.linked_pixel= linked_pixel;
    m_hits.total_pixel = total_pixel;
    m_hits.linked_sct  = linked_sct;
    m_hits.total_sct   = total_sct;
    m_hits.linked_trt  = linked_trt;
    m_hits.total_trt   = total_trt;
    // Add-up
    m_hits.linked_silicon = linked_pixel+linked_sct;
    m_hits.total_silicon  = total_pixel+total_sct;
    m_hits.linked_all     = m_hits.linked_silicon+linked_trt;
    m_hits.total_all      = m_hits.total_silicon+total_trt;
    // TO BE DPRECATED
    //ATH_MSG_INFO("INSIDE fillHITSTREE: "  
    //        << linked_pixel << " " << total_pixel << " "
    //        << linked_sct << " " << total_sct << " "
    //        << (linked_pixel+linked_sct) << " " << (total_pixel+total_sct) << " "
    //        << linked_trt << " " << total_trt << " "
    //        << (linked_pixel+linked_sct+linked_trt) << " " << (total_pixel+total_sct+total_trt));

    m_tree->Fill();
}

void ExtractTruthPRDInfo::fillHitsTree(const int & evtNumber, const int & runNumber,
        const std::pair<int,int> & pixel,
        const std::pair<int,int> & sct,
        const std::pair<int,int> & trt)
{
    fillHitsTree(evtNumber,runNumber,
            pixel.first,pixel.second,sct.first,sct.second,trt.first,trt.second);
}
