# Author: Jordi Duarte-Campderros
# Tel Aviv University, November 25, 2015

# Event selector, service messenger, common config flags
import AthenaPoolCnvSvc.ReadAthenaPool
from AthenaCommon.AppMgr import ServiceMgr
from AthenaCommon.AthenaCommonFlags import athenaCommonFlags

athenaCommonFlags.FilesInput = [ "root://eosatlas.cern.ch//eos/atlas/user/k/kmotohas/DisplacedVertex/valid1.159064.ParticleGenerator_K0S_Pt10.recon.ESD.e3099_s2578_r7297/ESD.07031345._000040.pool.root.1" ]

# Auto config
import AthenaPython.ConfigLib as apcl
cfg = apcl.AutoCfg(name = 'MyAutoConfig', input_files=athenaCommonFlags.FilesInput())
cfg.configure_job()

# Dealing with the detector elements associated to the hits (TSoS in tracks)
include("RecExCond/AllDet_detDescr.py")

from AthenaCommon.AlgSequence import AlgSequence
topSequence = AlgSequence()

from PersonalUtils.PersonalUtilsConf import SimpleTrackChecks
simpleTrackChecks = SimpleTrackChecks('SimpleTrackChecks')
############# The properties of the SimpleTrackChecks Algorithm
simpleTrackChecks.InputTracks       = "Tracks"
#simpleTrackChecks.OutputLevel       = VERBOSE
topSequence += simpleTrackChecks

# Set output level threshold (2=DEBUG, 3=INFO, 4=WARNING, 5=ERROR, 6=FATAL )
ServiceMgr.MessageSvc.OutputLevel = INFO
ServiceMgr.MessageSvc.defaultLimit = 9999999

## Configure input
ServiceMgr.EventSelector.InputCollections = athenaCommonFlags.FilesInput()
theApp.EvtMax = 10
ServiceMgr.EventSelector.SkipEvents=0
