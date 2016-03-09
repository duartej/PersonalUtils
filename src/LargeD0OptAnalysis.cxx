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
#include "GaudiKernel/IPartPropSvc.h"


// HepMC
#include "HepMC/GenVertex.h"
#include "HepMC/GenParticle.h"
#include "GeneratorObjects/McEventCollection.h"
#include "HepPDT/ParticleDataTable.hh"

// Athena Specific
#include "xAODEventInfo/EventInfo.h"
//

#include "TrkSpacePoint/SpacePointContainer.h"
#include "SiSpacePoint/SCT_SpacePoint.h"
#include "SiSpacePoint/PixelSpacePoint.h"
#include "TrkTruthData/PRD_MultiTruthCollection.h"
#include "TrkPrepRawData/PrepRawData.h"

//#include "TrkTrack/TrackCollection.h"
//#include "TrkTrack/TrackStateOnSurface.h"
#include "TrkTrackSummary/TrackSummary.h"
#include "TrkTruthData/TrackTruthCollection.h"

// ROOT
#include "TTree.h"

// System
#include <math.h>
#include <functional>
#include <iostream>
#include <algorithm>

// Auxiliary function to fill the SpacePoints related variable
namespace spacepoints_functions
{
    void fillSpacePointTree(const SpacePointContainer * sp_container,
            const PRD_MultiTruthCollection * prdTruthCollection, 
            std::unordered_set<const HepMC::GenParticle *> daughters_llp,
            const int & evtNumber,SPPosition & sp_struct,TTree * tree)
    {
        for(auto & sp_collection : *sp_container )
        {
            for(auto & sp : *sp_collection)
            {
                sp_struct.evtNumber = evtNumber;
                sp_struct.x = sp->globalPosition().x();
                sp_struct.y = sp->globalPosition().y();
                sp_struct.z = sp->globalPosition().z();
                sp_struct.r = sp->r();
                sp_struct.phi = sp->phi();
                sp_struct.clusterFromSignal = 0;
                // covariance matrix
                sp_struct.sigma_xx = sp->globCovariance()(0,0);
                sp_struct.sigma_xy = sp->globCovariance()(0,1);            
                sp_struct.sigma_xz = sp->globCovariance()(0,2);
                sp_struct.sigma_yy = sp->globCovariance()(1,1);
                sp_struct.sigma_yz = sp->globCovariance()(1,2);
                sp_struct.sigma_zz = sp->globCovariance()(2,2);
                
                // Get used clusters to build the SP
                const std::pair<const Trk::PrepRawData*,const Trk::PrepRawData*> clusters = sp->clusterList();
                // were these clusters generated by our signal particles?
                // ------ note: the returned iterators are HepMcParticleLink, you can retrieve the MC particle 
                //              pointer by the cptr method, for instance
                for(PRD_MultiTruthCollection::const_iterator mc = prdTruthCollection->find(clusters.first->identify()); 
                        mc != prdTruthCollection->end(); ++mc)
                {
                    // the cluster was generated at least by one signal particle, 
                    // and the particleso don't care if it was 
                    // more than one signal particle... break the for
                    if( daughters_llp.find( (*mc).second.cptr()  ) != daughters_llp.end() )
                    {
                        sp_struct.clusterFromSignal = 1; 
                    }
                }
                tree->Fill();
            }  // for SPs
        }  // for SP container
    }  // end fillSpacePointTree
}  // namespace

//////////////////////////////////////////////////////////////////////////////////////
/// Constructor

LargeD0OptAnalysis::LargeD0OptAnalysis(const std::string& name, ISvcLocator* pSvcLocator) :
    AthAlgorithm(name, pSvcLocator),
    m_processedEvts(0),
    m_mcCollName("TruthEvent"),
    m_pixelSPName("PixelSpacePoints"),
    m_SCTSPName("SCT_SpacePoints"),
    m_OverlapSPName("OverlapSpacePoints"),
    m_PixelPRDTruth("PRD_MultiTruthPixel"),
    m_SCTPRDTruth("PRD_MultiTruthSCT"),
    m_SiSeededTrackName("SiSpSeededLargeD0Tracks"),
    m_SiSeededTrackTruthName("SiSpSeededLargeD0TracksTruthCollection"),
    m_ResolvedTrackName("ResolvedLargeD0Tracks"),
    m_ResolvedTrackTruthName("ResolvedLargeD0TracksTruthCollection"),
    m_ExtendedTrackName("ExtendedLargeD0Tracks"),
    m_ExtendedTrackTruthName("ExtendedLargeD0TracksTruthCollection"),
    m_trackCollectionName(""),
    m_trackTruthCollectionName(""),
    m_streamHist("SP_OPT"),
    m_tree(nullptr),
    m_tree_pixel(nullptr),
    m_track_evtNumber(nullptr),
    m_track_genMatched(nullptr),
    m_track_charge(nullptr),
    m_track_pdgId(nullptr),
    m_track_d0(nullptr),
    m_track_z0(nullptr),
    m_track_pt(nullptr),
    m_track_phi0(nullptr),
    m_track_eta(nullptr),
    m_track_prob(nullptr),
    m_track_radiusFirstHit(nullptr),
    m_track_hits_nBLayers(nullptr),
    m_track_hits_nPixelHoles(nullptr),
    m_track_hits_nPixelHits(nullptr),
    m_track_hits_nGangedPixel(nullptr),
    m_track_hits_nSCTHits(nullptr),
    m_track_hits_nSCTHoles(nullptr),
    m_track_hits_nOutliers(nullptr),
    m_track_hits_nBLayerSharedHits(nullptr),
    m_track_hits_nPixelSharedHits(nullptr),
    m_track_hits_nSCTSharedHits(nullptr),
    m_track_hits_nSCTDoubleHoles(nullptr),
    m_track_hits_nContribPixelLayers(nullptr),
    m_track_hits_nBLayerOutliers(nullptr),
    m_vector_reg_tracks_int{&m_track_genMatched, &m_track_charge, 
                     &m_track_pdgId,&m_track_evtNumber,
                     &m_track_hits_nBLayers, &m_track_hits_nPixelHoles,
                     &m_track_hits_nPixelHits,&m_track_hits_nGangedPixel, 
                     &m_track_hits_nSCTHits, &m_track_hits_nSCTHoles,
                     &m_track_hits_nOutliers,&m_track_hits_nBLayerSharedHits, 
                     &m_track_hits_nPixelSharedHits,&m_track_hits_nSCTSharedHits,
                     &m_track_hits_nSCTDoubleHoles,&m_track_hits_nContribPixelLayers,
                     &m_track_hits_nBLayerOutliers},
    m_vector_reg_tracks_float{&m_track_d0, &m_track_z0, &m_track_pt, &m_track_phi0, 
                     &m_track_prob, &m_track_eta,&m_track_radiusFirstHit},
    m_tree_recoTracks(nullptr),
    m_dv_x(nullptr),
    m_dv_y(nullptr),
    m_dv_z(nullptr),
    m_dv_betagamma(nullptr),
    m_genpart_evtNumber(nullptr),
    m_genpart_trackMatched(nullptr),
    m_genpart_pdgId(nullptr),
    m_genpart_charge(nullptr),
    m_genpart_dvID(nullptr),
    m_genpart_vx(nullptr),
    m_genpart_vy(nullptr),
    m_genpart_vz(nullptr),
    m_genpart_pt(nullptr),
    m_genpart_eta(nullptr),
    m_genpart_phi(nullptr), 
    m_vector_reg_genpart_int{&m_genpart_evtNumber,
                                   &m_genpart_trackMatched, &m_genpart_pdgId, 
                                   &m_genpart_dvID},
    m_vector_reg_genpart_float{&m_dv_x, &m_dv_y, &m_dv_z, &m_dv_betagamma,
                                   &m_genpart_vx, &m_genpart_vy, 
                                   &m_genpart_vz, &m_genpart_pt, &m_genpart_eta,
                                   &m_genpart_phi},
    m_tree_mc(nullptr),
    m_populatedCache(0),
    m_pdg_LLP(1000022),
    m_pdt(nullptr)
{
    /** switches to control the analysis through job options */
    declareProperty( "MCTruthCollection",  m_mcCollName   );
    declareProperty( "PdgLLP",             m_pdg_LLP );
    
    declareProperty( "PixelSpacePoints",   m_pixelSPName   );
    declareProperty( "SCTSpacePoints",     m_SCTSPName     );
    declareProperty( "OverlapSpacePoints", m_OverlapSPName );

    // truth PRD
    declareProperty( "PixelClustersTruth", m_PixelPRDTruth   );
    declareProperty( "SCTClustersTruth",   m_SCTPRDTruth   );

    declareProperty( "SiSeededTracks",     m_SiSeededTrackName );
    declareProperty( "SiSeededTruthTracks",m_SiSeededTrackTruthName );
    declareProperty( "ResolvedTracks",     m_ResolvedTrackName );
    declareProperty( "ResolvedTruthTracks",m_ResolvedTrackTruthName );
    declareProperty( "ExtendedTracks",     m_ExtendedTrackName );
    declareProperty( "ExtendedTruthTracks",m_ExtendedTrackTruthName );
    
    //declareProperty( "SelectedTracks",     m_trackCollectionName = m_SiSeededTrackName );
    //declareProperty( "SelectedTruthTracks",m_trackTruthCollectionName = m_SiSeededTrackTruthName );
    declareProperty( "SelectedTracks",     m_trackCollectionName = "CombinedInDetTracks" );
    declareProperty( "SelectedTruthTracks",m_trackTruthCollectionName = "CombinedInDetTracksTruthCollection" );
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
    
    // Get the Particle Properties Service
    ServiceHandle<IPartPropSvc> partPropSvc("PartPropSvc", this->name());
    if ( !partPropSvc.retrieve().isSuccess() ) 
    {
        ATH_MSG_ERROR(" Could not initialize Particle Properties Service");
        return StatusCode::FAILURE;
    }      
 
    m_pdt = partPropSvc->PDT();
    if ( m_pdt == nullptr ) 
    {
        ATH_MSG_ERROR("Could not retrieve HepPDT::ParticleDataTable from "\
                "ParticleProperties Service !!");
        return StatusCode::FAILURE;
    }

    // Initialization and registration of the histograms/tree
    m_tree = new TTree("SpacePoint_STEP","SpacePoint finder step" );
    CHECK( tHistSvc->regTree(std::string("/"+m_streamHist+"/SPacePoints_STEP").c_str(),m_tree) );

    m_tree->Branch("sct",&m_lightSP,"evtNumber/I:clusterFromSignal/I:"\
            "x/F:y/F:z/F:r/F:phi/F:sigma_xx/F:sigma_xy/F:sigma_xz/F:sigma_yy/F:sigma_yz/F:sigma_zz/F");
    
    m_tree_pixel = new TTree("SpacePoint_PIXEL","SpacePoint finder step for pixel" );
    CHECK( tHistSvc->regTree(std::string("/"+m_streamHist+"/SPacePoints_PIXEL").c_str(),m_tree_pixel) );

    m_tree_pixel->Branch("pixel",&m_lightSP_Pixel,"evtNumber/I:clusterFromSignal/I:"\
            "x/F:y/F:z/F:r/F:phi/F:sigma_xx/F:sigma_xy/F:sigma_xz/F:sigma_yy/F:sigma_yz/F:sigma_zz/F");
    
    // Initialization and registration of the histograms/tree
    m_tree_mc = new TTree("GenTree","MC Truth" );
    CHECK( tHistSvc->regTree(std::string("/"+m_streamHist+"/mctruth").c_str(),m_tree_mc) );

    m_tree_mc->Branch("evtNumber",&m_genpart_evtNumber);
    m_tree_mc->Branch("dv_x",&m_dv_x);
    m_tree_mc->Branch("dv_y",&m_dv_y);
    m_tree_mc->Branch("dv_z",&m_dv_z);
    m_tree_mc->Branch("dv_betagamma",&m_dv_betagamma);
    
    m_tree_mc->Branch("trackMatched",&m_genpart_trackMatched);
    m_tree_mc->Branch("pdgId",&m_genpart_pdgId);
    m_tree_mc->Branch("charge",&m_genpart_charge);
    m_tree_mc->Branch("dvID",&m_genpart_dvID);
    m_tree_mc->Branch("vx",&m_genpart_vx);
    m_tree_mc->Branch("vy",&m_genpart_vy);
    m_tree_mc->Branch("vz",&m_genpart_vz);
    m_tree_mc->Branch("pt",&m_genpart_pt);
    m_tree_mc->Branch("eta",&m_genpart_eta);
    m_tree_mc->Branch("phi",&m_genpart_phi);

    // Initialization and registration: recoTracks
    m_tree_recoTracks = new TTree("RecoTracks",m_trackCollectionName.c_str() );
    CHECK( tHistSvc->regTree(std::string("/"+m_streamHist+"/recoTracks").c_str(),m_tree_recoTracks) );

    m_tree_recoTracks->Branch("evtNumber",&m_track_evtNumber);
    m_tree_recoTracks->Branch("genMatched",&m_track_genMatched);
    m_tree_recoTracks->Branch("charge",&m_track_charge);
    m_tree_recoTracks->Branch("pdgId",&m_track_pdgId);
    m_tree_recoTracks->Branch("d0",&m_track_d0);
    m_tree_recoTracks->Branch("z0",&m_track_z0);
    m_tree_recoTracks->Branch("pt",&m_track_pt);
    m_tree_recoTracks->Branch("phi0",&m_track_phi0);
    m_tree_recoTracks->Branch("eta",&m_track_eta);
    m_tree_recoTracks->Branch("prob",&m_track_prob);
    // hits
    m_tree_recoTracks->Branch("radiusFirstHit",&m_track_radiusFirstHit);
    m_tree_recoTracks->Branch("nBLayers",&m_track_hits_nBLayers);
    m_tree_recoTracks->Branch("nPixelHoles",&m_track_hits_nPixelHoles);
    m_tree_recoTracks->Branch("nPixelHits",&m_track_hits_nPixelHits);
    m_tree_recoTracks->Branch("nGangedPixel",&m_track_hits_nGangedPixel);
    m_tree_recoTracks->Branch("nSCTHits",&m_track_hits_nSCTHits); 
    m_tree_recoTracks->Branch("nSCTHoles",&m_track_hits_nSCTHoles);
    m_tree_recoTracks->Branch("nOutliers",&m_track_hits_nOutliers);
    m_tree_recoTracks->Branch("nBLayerOutliers",&m_track_hits_nBLayerOutliers);
    m_tree_recoTracks->Branch("nBLayerSharedHits",&m_track_hits_nBLayerSharedHits); 
    m_tree_recoTracks->Branch("nPixelSharedHits",&m_track_hits_nPixelSharedHits);
    m_tree_recoTracks->Branch("nSCTSharedHits",&m_track_hits_nSCTSharedHits);
    m_tree_recoTracks->Branch("nSCTDoubleHoles",&m_track_hits_nSCTDoubleHoles);
    m_tree_recoTracks->Branch("nContribPixelLayers",&m_track_hits_nContribPixelLayers);


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
    // Empty and cleaning variables
    m_daughters_llp.clear();
    m_populatedCache = 0;
    m_current_EvtNumber = -1;
    m_particleTrack_map.clear();

    ATH_MSG_DEBUG(" in execute()");

    StatusCode sc;

    // Event info related 
    const xAOD::EventInfo * evtInfo(0);
    sc = evtStore()->retrieve(evtInfo, "EventInfo");
    if(sc.isFailure()) 
    {
        ATH_MSG_FATAL("unable to retrieve the EventInfo");
        return StatusCode::FAILURE;
    }
    m_current_EvtNumber = evtInfo->eventNumber();

    
    // FIXME CHECK IF IS MC
    const McEventCollection* mcColl(0);
    sc = evtStore()->retrieve(mcColl, m_mcCollName);
    if(sc.isFailure()) 
    {
        ATH_MSG_FATAL("unable to retrieve MC coll");
        return StatusCode::FAILURE;
    }
    
    // Track-related --------------------------------------------------------

    //// Obtain the reco tracks
    const TrackCollection * trkCol(0);
    sc = evtStore()->retrieve(trkCol, m_trackCollectionName);
    if(sc.isFailure()) 
    {
        ATH_MSG_FATAL("unable to retrieve the Track Collection " << m_trackCollectionName);
        return StatusCode::FAILURE;
    }
    
    //// and the track-truth collection
    const TrackTruthCollection* truthMap(0);
    sc = evtStore()->retrieve(truthMap, m_trackTruthCollectionName);
    if(sc.isFailure()) 
    {
        ATH_MSG_FATAL("unable to retrieve the TrackTruth Collection " << m_trackTruthCollectionName);
        return StatusCode::FAILURE;
    }

    // --> XXX Be careful the order calling the functions is important!!!
    //     The algorithm should be re-designed, there is inconsistencies...
    //     FIXME: Anyway, the algorithm should call tools (to be implemented)
    
    // Create the map: GenParticle (HepMCParticleLink) to reco Track
    this->createGenParticleTrackMap(trkCol,truthMap);

    // Monte Carlo truth related ------------------------------------------------
    // --- Get the number of MC-truth final charged particles 
    //     from the displaced vertex (with some kinetic cuts)
    const std::unordered_set<const HepMC::GenParticle *> dv_gen_charged = getSignalChargedParticles(mcColl);

    
    // Store reco tracks info
    this->storeRecoTracksInfo(trkCol,truthMap);
    
    
    // --- Truth hits: SCT
    const PRD_MultiTruthCollection * prdTruthCollection = nullptr;
    sc = evtStore()->retrieve(prdTruthCollection, m_SCTPRDTruth);
    if(sc.isFailure())
    {
        ATH_MSG_FATAL("SCT PRD_MultiTruthCollection "<<m_SCTPRDTruth<<" NOT found");
        return StatusCode::FAILURE;
    } 
    else 
    {
        ATH_MSG_DEBUG ("Got SCT PRD_MultiTruthCollection "<<m_SCTPRDTruth);
    }

    // --------------------------------------------------------------------------
    // SpacePoints: SCT
    const SpacePointContainer *  sct_sp_container = nullptr;
    sc = evtStore()->retrieve( sct_sp_container, m_SCTSPName);
    if( sc.isFailure()  ||  ! sct_sp_container ) 
    {
        ATH_MSG_FATAL("No SpacePointContainer found with name '" 
                << m_SCTSPName << " !");
        return StatusCode::FAILURE;
    }

    // Retrieve space point info SCT
    spacepoints_functions::fillSpacePointTree(sct_sp_container,
            prdTruthCollection,m_daughters_llp,
            m_current_EvtNumber,m_lightSP,m_tree);
    
    // --- Truth hits: Pixel
    const PRD_MultiTruthCollection * pixel_prdTruthCollection = nullptr;
    sc = evtStore()->retrieve(pixel_prdTruthCollection, m_PixelPRDTruth);
    if(sc.isFailure())
    {
        ATH_MSG_FATAL("Pixel PRD_MultiTruthCollection "<<m_PixelPRDTruth<<" NOT found");
        return StatusCode::FAILURE;
    } 
    else 
    {
        ATH_MSG_DEBUG ("Got Pixel PRD_MultiTruthCollection "<<m_PixelPRDTruth);
    }
    // SpacePoints: Pixel
    const SpacePointContainer *  pixel_sp_container = nullptr;
    sc = evtStore()->retrieve( pixel_sp_container, m_pixelSPName);
    if( sc.isFailure()  ||  ! pixel_sp_container ) 
    {
        ATH_MSG_FATAL("No SpacePointContainer found with name '" 
                << m_pixelSPName << " !");
        return StatusCode::FAILURE;
    }
    // Retrieve space point info Pixel
    spacepoints_functions::fillSpacePointTree(pixel_sp_container,
            pixel_prdTruthCollection,
            m_daughters_llp,m_current_EvtNumber,
            m_lightSP_Pixel,m_tree_pixel);
    
    
    ATH_MSG_DEBUG("End processed event " << m_processedEvts);
  
    return StatusCode::SUCCESS;
}


const std::unordered_set<const HepMC::GenParticle *> LargeD0OptAnalysis::getSignalChargedParticles(const McEventCollection * mcColl)
{
    if( m_populatedCache > 1 )
    {
        return m_daughters_llp;
    }
    
    // -- allocate and set initial values to MC truth tree variables
    this->prepareMCTruthTreeRelated();

    // Event loop
    for(auto & evtItr : *mcColl)
    {
        if( evtItr->event_number() == -1 )
        {
            continue;
        }
        //m_current_EvtNumber = evtItr->event_number();

        ATH_MSG_DEBUG(" + MC Event Number: " << evtItr->event_number());
        ATH_MSG_DEBUG(" + Signal process ID: " << evtItr->signal_process_id());
        ATH_MSG_DEBUG(" + Number of particles: " << evtItr->particles_size());
        ATH_MSG_DEBUG(" + Number of vertices: " << evtItr->vertices_size());

        for(HepMC::GenEvent::particle_const_iterator it_p = evtItr->particles_begin();
                              it_p != evtItr->particles_end(); ++it_p) 
        {
            // Find the LLP
            if( (*it_p) && m_pdg_LLP == abs((*it_p)->pdg_id()) && 
                    (*it_p)->end_vertex()->particles_out_size () > 1 && 
                    (*it_p)->end_vertex()->particles_in_size () == 1)
            {
                // Store displaced-vertex,
                storeDV( (*it_p)->end_vertex() );
                // Store the daughters
                this->cacheDaughters( (*it_p) );
                ++m_populatedCache;
            }
        }
    }
    m_tree_mc->Fill();
    
    ATH_MSG_DEBUG("Number of final state, charged, pt > 1 GeV, gen-particles decayed from a LLP (" 
            << m_pdg_LLP << "): " << m_daughters_llp.size());
    //ATH_MSG_DEBUG("Number of final state, charged, pt > 1 GeV, gen-particles decayed from a LLP (" 
    //        << m_pdg_LLP << ") with a reconstructed track: " << m_particleTrack_map.size());
    ATH_MSG_DEBUG("Number of final state, charged, pt > 1 GeV, gen-particles decayed from a LLP (" 
            << m_pdg_LLP << ") with a reconstructed track: " 
            << std::count_if(m_genpart_trackMatched->begin(), m_genpart_trackMatched->end(), [](int i) {return i == 1;}));

    // -- deallocate variables related to MC truth tree
    this->deallocateMCTruthTreeRelated();
    
    return m_daughters_llp;
} 

void LargeD0OptAnalysis::cacheDaughters(const HepMC::GenParticle * p)
{
    //just basic check...
    if( !p ) 
    {
        return;
    }
    // store this particle (if it is not the main mother), 
    // pass minimum cuts (pt, charged, ...) and it wasn't stored before
    if( abs(p->pdg_id()) != m_pdg_LLP && p->is_undecayed() && passKinCuts(p) 
            && m_daughters_llp.find(p) == m_daughters_llp.end() )
    {
        m_daughters_llp.emplace(p);
        // store in the tree, 
        this->storeGenParticleInfo(p);
        ATH_MSG_DEBUG(*p);
    }
    // get the decay vertex to obtain the parent
    const HepMC::GenVertex * decayVtx = p->end_vertex();
    if( decayVtx == nullptr )
    {
        return;
    }
    // daughters recursive calls
    for(HepMC::GenVertex::particles_out_const_iterator it_daughter =
            decayVtx->particles_out_const_begin();      
            it_daughter != decayVtx->particles_out_const_end(); ++it_daughter) 
    {
        cacheDaughters(*it_daughter);
    }
}

bool LargeD0OptAnalysis::passKinCuts(const HepMC::GenParticle * p)
{
    // pt in MeV
    if( p->momentum().perp() < 1000.0 || abs(p->momentum().eta()) > 5.0 )
    {
        return false;
    }

    // Not so far away from the vertex ?
    /*const float vx_rel = p->production_vertex()->point3d().x()-(*m_dv_x)[m_current_dvID];
    const float vy_rel = p->production_vertex()->point3d().y()-(*m_dv_y)[m_current_dvID];
    const float radius = sqrt(vx_rel*vx_rel+vy_rel*vy_rel);
    //FIXME:: declare as properties
    const float m_rmin = 0.0;
    const float m_rmax = 30.0; // mm
    if( radius < m_rmin || radius > m_rmax )
    {
        return false;
    }*/
    // Only charged particles
    const HepPDT::ParticleData * pd = m_pdt->particle(abs(p->pdg_id()));
    if( pd == nullptr)
    {
        return false;
    }

    if( abs(pd->charge()) < 0.5 )
    {
        return false;
    }

    return true;
}

void LargeD0OptAnalysis::createGenParticleTrackMap(const TrackCollection * recoTracks, const TrackTruthCollection *truthMap)
{
    // TrackTruthCollection is inheriting from std::map<Trk::TrackTruthKey, TrackTruth> 
    for(auto & kv_pair: *truthMap)
    {
        size_t trackIndex = kv_pair.first.index();
        m_particleTrack_map.emplace( kv_pair.second.particleLink().cptr(), recoTracks->at(trackIndex) ); 
    }

    return;
}

void LargeD0OptAnalysis::storeRecoTracksInfo(const TrackCollection * recoTracks, const TrackTruthCollection * truthMap)
{
    this->prepareTrackTreeRelated(recoTracks->size());

    // How to avoid duplicated tracks??
    // XXX: Be carefull!! Seems duplicated tracks because the front-TSOS
    // contains exactly the same values, but the track is defined by 
    // more things, right? It should be checked
    for(size_t i = 0; i < recoTracks->size(); ++i)
    {
        const Trk::Track * track = (*recoTracks)[i];
        if( track == nullptr )
        {
            continue;
        }

        // event number
        m_track_evtNumber->push_back( m_current_EvtNumber );

        const DataVector< const Trk::TrackStateOnSurface > * tsos = track->trackStateOnSurfaces();
        // first hit
        bool fhfilled = false;
        for( auto & _ts : *tsos)
        {
            if( _ts->type(Trk::TrackStateOnSurface::Measurement) )
            {
                // Asumming from inside to outside
                const Amg::Vector3D pfh = _ts->trackParameters()->position();
                m_track_radiusFirstHit->push_back(sqrt(pfh[0]*pfh[0]+pfh[1]*pfh[1]));
                fhfilled = true;
                break;
            }
        }
        if( ! fhfilled )
        {
            m_track_radiusFirstHit->push_back(-1.);
        }

        //const Trk::Perigee * perigee = track->perigeeParameters();
        // Check if the perigee is available, if not we need to extrapolate it
        //if( perigee == nullptr )
        bool pfilled = false;
        if( tsos != nullptr )
        {
            // Fill the radius first hit for the first measurement type
            if(tsos->front() != nullptr)
            {
                const Trk::TrackParameters * perigee = tsos->front()->trackParameters();
                if( perigee != nullptr )
                {
                    m_track_d0->push_back( perigee->parameters()[Trk::d0] );
                    m_track_z0->push_back( perigee->parameters()[Trk::z0] );
                    m_track_pt->push_back( perigee->pT()*1e-3 );
                    m_track_phi0->push_back( perigee->parameters()[Trk::phi0] );
                    m_track_eta->push_back( perigee->eta() );
                    m_track_charge->push_back( perigee->charge() );
                    pfilled = true;
                }
            }
        }
        if( ! pfilled )
        //if( perigee == nullptr )
        {
            ATH_MSG_DEBUG("Perigee track parameters not found in " << i << "-track");
            m_track_d0->push_back( 1e10 );
            m_track_z0->push_back( 1e10 );
            m_track_pt->push_back( -1. );
            m_track_phi0->push_back( 1e10 );
            m_track_eta->push_back( 1e10 );
            m_track_charge->push_back( 0 );
        }
        /*else
        {
            m_track_d0->push_back( perigee->parameters()[Trk::d0] );
            m_track_z0->push_back( perigee->parameters()[Trk::z0] );
            m_track_pt->push_back( perigee->pT()*1e-3 );
            m_track_phi0->push_back( perigee->parameters()[Trk::phi0] );
            m_track_eta->push_back( perigee->eta() );
            m_track_charge->push_back( perigee->charge() );
        }*/
        
        // Has this track a gen-particle associated?
        // using the TrackTruth, where the containing particleLink is the 'bestMatch'
        // (when several GenParticle contribute to the hits of the track, the particle
        // with largest energy deposit is chosen )
        // Note that if there is no genparticle link (or equivalently the particle barcode
        // is 0), the track is considered coming from Noise (or pile-up),
        // The probability is related with a TruthTrack (no necessarely related to a GenParticle)?
        TrackTruthCollection::const_iterator it_truthMap = truthMap->find( Trk::TrackTruthKey(i) );
        if( it_truthMap != truthMap->end() )
        {
            if( (*it_truthMap).second.particleLink().isValid() )
            {
                m_track_genMatched->push_back(1);
                m_track_pdgId->push_back( (*it_truthMap).second.particleLink()->pdg_id() );
                m_track_prob->push_back( (*it_truthMap).second.probability() );
            }
            else
            {
                ATH_MSG_DEBUG("No GenParticle link found! Barcode: " 
                        << (*it_truthMap).second.particleLink().barcode());
                m_track_genMatched->push_back(0);
                m_track_pdgId->push_back( 0 );
                m_track_prob->push_back( (*it_truthMap).second.probability() );

            }
        }
        else
        {
            ATH_MSG_DEBUG("No GenParticle link found, noLink to " << i << "-track"); 
            m_track_genMatched->push_back(0);
            m_track_pdgId->push_back( 0 );
            m_track_prob->push_back( -1. );
        }
        // Hits and other stuff, ...
        const Trk::TrackSummary * track_summary = track->trackSummary();
        if( track_summary != nullptr )
        {
            m_track_hits_nBLayers->push_back(track_summary->get(Trk::SummaryType::numberOfBLayerHits));
            m_track_hits_nPixelHoles->push_back(track_summary->get(Trk::SummaryType::numberOfPixelHoles));
            m_track_hits_nPixelHits->push_back(track_summary->get(Trk::SummaryType::numberOfPixelHits));
            m_track_hits_nGangedPixel->push_back(track_summary->get(Trk::SummaryType::numberOfGangedPixels));
            m_track_hits_nSCTHits->push_back(track_summary->get(Trk::SummaryType::numberOfSCTHits));
            m_track_hits_nSCTHoles->push_back(track_summary->get(Trk::SummaryType::numberOfSCTHoles));
            m_track_hits_nOutliers->push_back(track_summary->get(Trk::SummaryType::numberOfOutliersOnTrack));
            m_track_hits_nBLayerSharedHits->push_back(track_summary->get(Trk::SummaryType::numberOfBLayerSharedHits));
            m_track_hits_nPixelSharedHits->push_back(track_summary->get(Trk::SummaryType::numberOfPixelSharedHits));
            m_track_hits_nSCTSharedHits->push_back(track_summary->get(Trk::SummaryType::numberOfSCTSharedHits));
            m_track_hits_nSCTDoubleHoles->push_back(track_summary->get(Trk::SummaryType::numberOfSCTDoubleHoles));
            m_track_hits_nContribPixelLayers->push_back(track_summary->get(Trk::SummaryType::numberOfContribPixelLayers));
            m_track_hits_nBLayerOutliers->push_back(track_summary->get(Trk::SummaryType::numberOfBLayerOutliers));
            //  numberOfBLayerOutliers = 
        }
    }
    m_tree_recoTracks->Fill();

    this->deallocateTrackTreeRelated();
}


void LargeD0OptAnalysis::storeDV( const HepMC::GenVertex * dv)
{
    // Store only status 2 vertex (in particles: 1 , out particles > 1)
    if( dv->particles_out_size() < 2 || dv->particles_in_size() != 1 )
    {
        return;
    }
    ATH_MSG_DEBUG("Displaced vertex at d0= " << dv->point3d().perp()
            << " mm, z0= " << dv->point3d().z());
    m_dv_x->push_back(dv->point3d().x());
    m_dv_y->push_back(dv->point3d().y());
    m_dv_z->push_back(dv->point3d().z());
    const HepMC::GenParticle * llp =  *(dv->particles_in_const_begin());
    m_dv_betagamma->push_back(llp->momentum().rho()/llp->momentum().m());
    ++m_current_dvID;
}

void LargeD0OptAnalysis::storeGenParticleInfo(const HepMC::GenParticle * p)
{
    /// evtNumber
    m_genpart_evtNumber->push_back(m_current_EvtNumber);
    /// -> found a track which origin is this gen particle
    if( m_particleTrack_map.find(p) != m_particleTrack_map.end() )
    {
        m_genpart_trackMatched->push_back(1);
    }
    else
    {
        m_genpart_trackMatched->push_back(0);
    }
    // m_genpart_trackMatched  = new std::vector<int>;
    m_genpart_pdgId->push_back(p->pdg_id());

    // Note that if we are here means that it is a charged particle
    if( p->pdg_id() < 0 )
    {
        m_genpart_charge->push_back(-1);
    }
    else
    {
        m_genpart_charge->push_back(1);
    }
    m_genpart_dvID->push_back(m_current_dvID);
    m_genpart_vx->push_back(p->production_vertex()->point3d().x() );
    m_genpart_vy->push_back(p->production_vertex()->point3d().y() );
    m_genpart_vz->push_back(p->production_vertex()->point3d().z()); 
    m_genpart_pt->push_back(p->momentum().perp()*1e-3);
    m_genpart_eta->push_back(p->momentum().eta());
    m_genpart_phi->push_back(p->momentum().phi()); 
}

void LargeD0OptAnalysis::prepareTrackTreeRelated(const unsigned int & recoSize)
{
    for(auto & trackPointerVectorI: m_vector_reg_tracks_int)
    {
        (*trackPointerVectorI) = new std::vector<int>;
        (*trackPointerVectorI)->reserve(recoSize);
    }
    for(auto & trackPointerVectorF: m_vector_reg_tracks_float)
    {
        (*trackPointerVectorF) = new std::vector<float>;
        (*trackPointerVectorF)->reserve(recoSize);
    }
}

void LargeD0OptAnalysis::deallocateTrackTreeRelated()
{
    for(auto & trackPointerVectorI: m_vector_reg_tracks_int)
    {
        if( *trackPointerVectorI )
        {
            delete *trackPointerVectorI;
            (*trackPointerVectorI) = nullptr;
        }
    }
    for(auto & trackPointerVectorF: m_vector_reg_tracks_float)
    {
        if( *trackPointerVectorF )
        {
            delete (*trackPointerVectorF);
            (*trackPointerVectorF) = nullptr;
        }
    }
}

void LargeD0OptAnalysis::prepareMCTruthTreeRelated()
{
    m_current_dvID = -1;
    for(auto & genPointerVectorI: m_vector_reg_genpart_int)
    {
        (*genPointerVectorI) = new std::vector<int>;
    }
    for(auto & genPointerVectorF: m_vector_reg_genpart_float)
    {
        (*genPointerVectorF) = new std::vector<float>;
    }
}

void LargeD0OptAnalysis::deallocateMCTruthTreeRelated()
{
    for(auto & genPointerVectorI: m_vector_reg_genpart_int)
    {
        if( *genPointerVectorI )
        {
            delete *genPointerVectorI;
            (*genPointerVectorI) = nullptr;
        }
    }
    for(auto & genPointerVectorF: m_vector_reg_genpart_float)
    {
        if( *genPointerVectorF )
        {
            delete (*genPointerVectorF);
            (*genPointerVectorF) = nullptr;
        }
    }
}
