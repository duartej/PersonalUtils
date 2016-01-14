////////////////////////////////////////////////////////////////////////////////// 
/// Just minor checks in teh Trk::Tracks
///
/// Author: Jordi Duarte-Campderrros
/// Tel Aviv University, November 25, 2015
///
/// DESCRIPTION: Perform some checks in the Trk::Tracks presents 
///
///
//////////////////////////////////////////////////////////////////////////////////

#include "PersonalUtils/SimpleTrackChecks.h"
// the first two come for free when using AthAlgorithm
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/IToolSvc.h"

//#include "xAODEventInfo/EventInfo.h"

#include "TrkTrack/TrackCollection.h"
#include "TrkSurfaces/Surface.h"
#include "TrkSurfaces/ConeSurface.h"
#include "TrkSurfaces/CylinderSurface.h"
#include "TrkSurfaces/DiscSurface.h"
#include "TrkSurfaces/PerigeeSurface.h"
#include "TrkSurfaces/PlaneSurface.h"
//#include "TrkSurfaces/LineSurface.h" 
//#include "TrkSurfaces/CurvilinearSurface.h"   

#include <math.h>
#include <functional>
#include <iostream>

//////////////////////////////////////////////////////////////////////////////////////
/// Constructor

SimpleTrackChecks::SimpleTrackChecks(const std::string& name, ISvcLocator* pSvcLocator) :
  AthAlgorithm(name, pSvcLocator),
  m_processedEvts(0),
  m_tracksContainerName("Tracks")
{
  /** switches to control the analysis through job options */
      declareProperty( "InputTracks",   m_tracksContainerName );
}

/////////////////////////////////////////////////////////////////////////////////////
/// Destructor - check up memory allocation
/// delete any memory allocation on the heap

SimpleTrackChecks::~SimpleTrackChecks() {}

////////////////////////////////////////////////////////////////////////////////////
StatusCode SimpleTrackChecks::beginRun() 
{ 
  return StatusCode::SUCCESS;
} 

StatusCode SimpleTrackChecks::initialize() 
{
  ATH_MSG_DEBUG("Initializing SimpleTrackChecks");
  return StatusCode::SUCCESS;
}		 

///////////////////////////////////////////////////////////////////////////////////
/// Finalize - delete any memory allocation from the heap

StatusCode SimpleTrackChecks::finalize() 
{  
  return StatusCode::SUCCESS;
}

//////////////////////////////////////////////////////////////////////////////////
/// Execute - on event by event

StatusCode SimpleTrackChecks::execute() 
{
  ++m_processedEvts;
  ATH_MSG_DEBUG(" in execute()");

  StatusCode sc;
  
  const TrackCollection *  trackcol = 0;
  sc = evtStore()->retrieve( trackcol, m_tracksContainerName);
  if( sc.isFailure()  ||  !trackcol ) 
  {
     ATH_MSG_ERROR("No Trk::TrackCollection found with name '" 
             << m_tracksContainerName << " !");
     return StatusCode::FAILURE;
  }

  int i_track = 0;
  for(auto & track: (*trackcol) )
  {
      ATH_MSG_INFO(i_track << "-track");
      //std::cout << *track << std::endl;
      int i_tsos = 0;
      for(auto & tsos: *(track->trackStateOnSurfaces()) )
      { 
          ATH_MSG_INFO("@track " << i_track << " " << i_tsos << "-TSOS");
          if( tsos->materialEffectsOnTrack() != 0 )
          {
              // indications that the problem is located in the associated surface of
              // the material effects class
              const Trk::Surface * surface = &(tsos->materialEffectsOnTrack()->associatedSurface());

              // which type? --> No problem with disc...
              ATH_MSG_INFO( "Surface Type: " << surface->type() );
              /*if( surface->type() == Trk::Surface::SurfaceType::Disc )
              {
                  std::cout << (*surface) << std::endl;
              }*/
              // The problem seems to be focalized at the Bounds of the cylinder... Right, NULL pointer.
              if( surface->type() == Trk::Surface::SurfaceType::Cylinder )
              {
                  std::cout << "Address for CylinderBounds: " << &(surface->bounds()) << std::endl;
              }
          }
          ++i_tsos;
      }
      ++i_track;
  }
  ATH_MSG_INFO("End processed event " << m_processedEvts);

  return StatusCode::SUCCESS;
}

