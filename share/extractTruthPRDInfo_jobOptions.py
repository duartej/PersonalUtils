# Author: Jordi Duarte-Campderros
# Tel Aviv University, November 25, 2015

# Event selector, service messenger, common config flags
import AthenaPoolCnvSvc.ReadAthenaPool
from AthenaCommon.AppMgr import ServiceMgr
from AthenaCommon.AthenaCommonFlags import athenaCommonFlags

athenaCommonFlags.FilesInput = ['root://eosatlas.cern.ch//eos/atlas/user/d/duarte/datasets/mc15_13TeV.403069.PythiaRhad_AUET2BCTEQ6L1_gen_gluino_p1_1400_qq_1300_1ns.recon.ESD.e5220_s2601_s2952_r8318/ESD.09247061._000012.pool.root.1' ]

# Auto config
import AthenaPython.ConfigLib as apcl
cfg = apcl.AutoCfg(name = 'MyAutoConfig', input_files=athenaCommonFlags.FilesInput())
cfg.configure_job()

# Dealing with the detector elements associated to the hits (TSoS in tracks)
include("RecExCond/AllDet_detDescr.py")

from AthenaCommon.AlgSequence import AlgSequence
topSequence = AlgSequence()

from PersonalUtils.PersonalUtilsConf import ExtractTruthPRDInfo
extractPRDTruthInfo = ExtractTruthPRDInfo("ExtractTruthPRDInfo")
############# The properties of the ExtractTruthPRDInfo Algorithm
## extractPRDTruthInfo.PRDTruthNamePixel = "PRD_MultiTruthPixel"
# ...
topSequence += extractPRDTruthInfo

# Set output level threshold (2=DEBUG, 3=INFO, 4=WARNING, 5=ERROR, 6=FATAL )
ServiceMgr.MessageSvc.OutputLevel = INFO
ServiceMgr.MessageSvc.defaultLimit = 9999999

## Configure input
ServiceMgr.EventSelector.InputCollections = athenaCommonFlags.FilesInput()
theApp.EvtMax = 10
ServiceMgr.EventSelector.SkipEvents=0
