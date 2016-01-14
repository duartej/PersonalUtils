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

#include <string>

class TTree;

// Auxiliary structure for the position variables of the
// SpacePoints
struct SPPosition 
{
    float x;
    float y;
    float z;
    float r;
    float phi;
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
        /** Number of processed events */
        int m_processedEvts;
        
        /** the key of the SpacePoint containers */
        std::string m_pixelSPName;
        std::string m_SCTSPName;
        std::string m_OverlapSPName;

        /** the key of the Track Containers */
        std::string m_SiSeededTrackName;
        std::string m_ResolvedTrackName;
        std::string m_ExtendedTrackName;

        /** Tree */
        std::string m_streamHist;
        SPPosition m_lightSP;
        TTree * m_tree;
};

#endif // LARGED0OPTANALYSIS_H

