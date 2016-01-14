#include "PersonalUtils/SimpleTrackChecks.h"
#include "PersonalUtils/LargeD0OptAnalysis.h"

#include "GaudiKernel/DeclareFactoryEntries.h"

DECLARE_ALGORITHM_FACTORY( SimpleTrackChecks )
DECLARE_ALGORITHM_FACTORY( LargeD0OptAnalysis )

DECLARE_FACTORY_ENTRIES( PersonalUtils ) {
  DECLARE_ALGORITHM( SimpleTrackChecks )
  DECLARE_ALGORITHM( LargeD0OptAnalysis )
}

