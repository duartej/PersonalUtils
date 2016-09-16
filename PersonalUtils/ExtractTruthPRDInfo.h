#ifndef EXTRACTTRUTHPRDINFO_H
#define EXTRACTTRUTHPRDINFO_H
/////////////////////////////////////////////////////////////////////////////////////////////////////
/// Name    : ExtractTruthPRDInfo.h
/// Package : PersonalUtils
/// Author  : Jordi Duarte-Campderros
/// Created : Sep. 2016
///
/// DESCRIPTION: Simple class to extract some information of the simulated hits (in particular
//               the truth particles that generated them, ... )
/// 
///
///////////////////////////////////////////////////////////////////////////////////////////////////////

//#include "StoreGate/StoreGateSvc.h"
#include "AthenaBaseComps/AthAlgorithm.h"

#include "TrkTruthData/PRD_MultiTruthCollection.h"
#include "InDetPrepRawData/PixelClusterContainer.h"
#include "InDetPrepRawData/SCT_ClusterContainer.h"
#include "InDetPrepRawData/TRT_DriftCircleContainer.h"


#include <string>

class TTree;

// Auxiliary structure to store truth-hit related info,
struct TruthHit
{
    int eventNumber;
    int runNumber;
    int linked_pixel;
    int total_pixel;
    int linked_sct;
    int total_sct;
    int linked_silicon;
    int total_silicon;
    int linked_trt;
    int total_trt;
    int linked_all;
    int total_all;
};
 
class ExtractTruthPRDInfo : public AthAlgorithm  
{
    public:
        ExtractTruthPRDInfo(const std::string& name, ISvcLocator* pSvcLocator);
        ~ExtractTruthPRDInfo();
        
        virtual StatusCode beginRun();
        virtual StatusCode initialize();
        virtual StatusCode finalize();
        virtual StatusCode execute();
        
    private:
        /** process reco data **/ 
        StatusCode processRecoInput(const int & evtNumber, const int & runNumber);
        StatusCode processHITSInput(const int & evtNumber, const int & runNumber);
        /** extract the number of sim hits with a gen-particle generated, and the total
            number of hits */
        const std::pair<int,int> getLinkedAndTotalPRDs(const PRD_MultiTruthCollection * const prdCol);
        /** extract the number of reco hits with a gen-particle generated, and the total
            number of hits */
        template<class PRD_Container> 
            const std::pair<int,int> getLinkedAndTotalPRDs(const PRD_Container * const prdContainer,
                    const PRD_MultiTruthCollection * const prdTruthCol);
        /** Print out a summary of hits in the event */
        const std::string getSummary() const;

        /** Fill tree related functions */
        void fillHitsTree(const int & evtNumber, const int & runNumber,
                const int & linked_pixel,const int & total_pixel,
                const int & linked_sct,const int & total_sct,
                const int & linked_trt,const int & total_trt);
        void fillHitsTree(const int & evtNumber, const int & runNumber,
                const std::pair<int,int> & pixel,
                const std::pair<int,int> & sct,
                const std::pair<int,int> & trt);

        /** Number of processed events */
        int m_processedEvts;

        // whether or not using a RECO file (or a HITS)
        bool m_recoMode;
        
        /** the key of the PRD Truth of the different subdetectors*/
        std::string m_prdTruthPixel;
        std::string m_prdTruthSCT;
        std::string m_prdTruthTRT;
        /** the key of the PRD of the different subdetectors*/
        std::string m_prdPixel;
        std::string m_prdSCT;
        std::string m_prdTRT;
        /** the key of the Sim hits of the different subdetectors*/
        std::string m_simHitsPixel;
        std::string m_simHitsSCT;
        std::string m_simHitsTRT;

        /* Tree related */
        std::string m_streamHist;
        TTree * m_tree;
        TruthHit m_hits;
};

#endif // EXTRACTTRUTHPRDINFO_H

