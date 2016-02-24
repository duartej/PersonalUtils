#ifndef LARGED0OPTANALYSIS_H
#define LARGED0OPTANALYSIS_H
/////////////////////////////////////////////////////////////////////////////////////////////////////
/// Name    : LargeD0OptAnalysis.h
/// Package : 
/// Author  : Jordi Duarte-Campderros
/// Created : Nov. 2015
///
/// DESCRIPTION: Simple class to extract some information of the tracks presents in a file
///
/// 
///
///////////////////////////////////////////////////////////////////////////////////////////////////////

//#include "StoreGate/StoreGateSvc.h"
#include "AthenaBaseComps/AthAlgorithm.h"

#include "TrkTrack/TrackCollection.h"
//#include "GeneratorObjects/HepMcParticleLink.h"

#include <string>
#include <unordered_set>
#include <map>

class TTree;

namespace HepPDT
{
    class ParticleDataTable;
}

namespace HepMC
{
    class GenVertex;
    class GenParticle;
}

class McEventCollection;

namespace Trk 
{
    class Track;
}

class TrackTruthCollection;

// Auxiliary structure for the position variables of the
// SpacePoints
struct SPPosition 
{
    int   evtNumber;
    int   nClustersFromSignal;
    float x;
    float y;
    float z;
    float r;
    float phi;
    float sigma_xx;
    float sigma_xy;
    float sigma_xz;
    float sigma_yy;
    float sigma_yz;
    float sigma_zz;    
};


class LargeD0OptAnalysis : public AthAlgorithm  
{
    public:
        LargeD0OptAnalysis(const std::string& name, ISvcLocator* pSvcLocator);
        ~LargeD0OptAnalysis();
        
        virtual StatusCode beginRun();
        virtual StatusCode initialize();
        virtual StatusCode finalize();
        virtual StatusCode execute();
        
    private:
        // -------- Some internal useful methods
        /** get the undecayed MC truth particles from a LLP */
        const std::unordered_set<const HepMC::GenParticle *> 
            getSignalChargedParticles(const McEventCollection * mcColl);
        
        /** populate the cache with the daughter particle of a llp */
        void cacheDaughters(const HepMC::GenParticle * p);
        /** particle cuts */
        bool passKinCuts(const HepMC::GenParticle *p);
        
        /** mctruth tree related */
        void createGenParticleTrackMap(const TrackCollection * recoTracks, const TrackTruthCollection * truthMap);
        void storeDV(const HepMC::GenVertex * dv);
        void storeGenParticleInfo(const HepMC::GenParticle *p);
        void prepareMCTruthTreeRelated();
        void deallocateMCTruthTreeRelated();
        
        /** reco track tree related */
        void storeRecoTracksInfo(const TrackCollection * recoTracks, const TrackTruthCollection * truthMap);
        void prepareTrackTreeRelated(const unsigned int & trackSize);
        void deallocateTrackTreeRelated();


        /** Number of processed events */
        int m_processedEvts;
        int m_current_EvtNumber;

        /** the key of the MC truth collection */
        std::string m_mcCollName;
        
        /** the key of the SpacePoint containers */
        std::string m_pixelSPName;
        std::string m_SCTSPName;
        std::string m_OverlapSPName;

        /** the key of the PRD truth container */
        std::string m_PixelPRDTruth;
        std::string m_SCTPRDTruth;

        /** the key of the Track Containers */
        std::string m_SiSeededTrackName;
        std::string m_SiSeededTrackTruthName;
        std::string m_ResolvedTrackName;
        std::string m_ResolvedTrackTruthName;
        std::string m_ExtendedTrackName;
        std::string m_ExtendedTrackTruthName;
        std::string m_trackCollectionName;
        std::string m_trackTruthCollectionName;

        /** SpacePoint Tree: SCT */
        std::string m_streamHist;
        SPPosition m_lightSP;
        TTree * m_tree;
        /** Pixel */
        SPPosition m_lightSP_Pixel;
        TTree * m_tree_pixel;
        
        /** recotrack tree */
        std::vector<int>   * m_track_evtNumber;
        std::vector<int>   * m_track_genMatched;
        std::vector<int>   * m_track_charge;
        std::vector<int>   * m_track_pdgId;
        std::vector<float> * m_track_d0;
        std::vector<float> * m_track_z0;
        std::vector<float> * m_track_pt;
        std::vector<float> * m_track_phi0;
        std::vector<float> * m_track_eta;
        std::vector<float> * m_track_prob;
        // hit related
        std::vector<int>   * m_track_hits_nBLayers;
        std::vector<int>   * m_track_hits_nPixelHoles;
        std::vector<int>   * m_track_hits_nPixelHits;
        std::vector<int>   * m_track_hits_nGangedPixel;
        std::vector<int>   * m_track_hits_nSCTHits;
        std::vector<int>   * m_track_hits_nSCTHoles;
        std::vector<int>   * m_track_hits_nOutliers;
        std::vector<int>   * m_track_hits_nBLayerSharedHits;
        std::vector<int>   * m_track_hits_nPixelSharedHits;
        std::vector<int>   * m_track_hits_nSCTSharedHits;
        std::vector<int>   * m_track_hits_nSCTDoubleHoles;
        std::vector<int>   * m_track_hits_nContribPixelLayers;
        std::vector<int>   * m_track_hits_nBLayerOutliers;
        std::vector<std::vector<int> **>   m_vector_reg_tracks_int; 
        std::vector<std::vector<float> **> m_vector_reg_tracks_float; 
        TTree * m_tree_recoTracks;

        /** MC truth tree */
        std::vector<float> * m_dv_x;
        std::vector<float> * m_dv_y;
        std::vector<float> * m_dv_z;
        std::vector<int>   * m_genpart_evtNumber;
        std::vector<int>   * m_genpart_trackMatched;
        std::vector<int>   * m_genpart_pdgId;
        std::vector<int>   * m_genpart_charge;
        std::vector<int>   * m_genpart_dvID;
        std::vector<float> * m_genpart_vx;
        std::vector<float> * m_genpart_vy;
        std::vector<float> * m_genpart_vz;
        std::vector<float> * m_genpart_pt;
        std::vector<float> * m_genpart_eta;
        std::vector<float> * m_genpart_phi; 
        std::vector<std::vector<int> **>   m_vector_reg_genpart_int; 
        std::vector<std::vector<float> **> m_vector_reg_genpart_float; 
        TTree * m_tree_mc;

        /** MC Truth particles info */
        std::unordered_set<const HepMC::GenParticle *> m_daughters_llp;
        int m_populatedCache;
        int m_pdg_LLP;
        int m_current_dvID;
        /** Gen particle link to reco track map */
        std::map<const HepMC::GenParticle *,const Trk::Track *> m_particleTrack_map;

        /** Particle data table */
        const HepPDT::ParticleDataTable * m_pdt;
};

#endif // LARGED0OPTANALYSIS_H

