#ifndef SIMPLETRACKCHECKS_H
#define SIMPLETRACKCHECKS_H
/////////////////////////////////////////////////////////////////////////////////////////////////////
/// Name    : SimpleTrackChecks.h
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

// upgrade to inherit from AthAlgorithm
 
class SimpleTrackChecks : public AthAlgorithm  
{
    public:
        SimpleTrackChecks(const std::string& name, ISvcLocator* pSvcLocator);
        ~SimpleTrackChecks();
        
        virtual StatusCode beginRun();
        virtual StatusCode initialize();
        virtual StatusCode finalize();
        virtual StatusCode execute();
        
    private:
        /** Number of processed events */
        int m_processedEvts;
        
        /** the key of the Track Container to retrieve from the AOD */
        std::string m_tracksContainerName;
};

#endif // SIMPLETRACKCHECKS_H

