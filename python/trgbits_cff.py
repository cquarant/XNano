import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.common_cff import *

trgTable = cms.EDProducer( "TrgBitTableProducer",
                           hltresults = cms.InputTag("TriggerResults::HLT"),

                           # interesting paths for 2016 
                           # paths      = cms.vstring(
                           #                   "HLT_Dimuon16_Jpsi",                        # inclusive dimuon jpsi
                           #                   "HLT_Dimuon10_Jpsi_Barrel",                 # inclusive dimuon jpsi
                           #                   "HLT_DoubleMu4_JpsiTrk_Displaced",          # displaced jpistrk 
                           #                   "HLT_Dimuon20_JpsiHLT_Dimuon20_Jpsi",       # 
                           #                   "HLT_DoubleMu4_3_Jpsi_Displaced",           # prescaled
                           #                   "HLT_DoubleMu4_3_Jpsi",                     # prescaled
                           #                   "HLT_DoubleMu4_Jpsi_Displaced",             # prescaled
                           #                   "HLT_Dimuon13_PsiPrime",                    # inclusive dimuon psi2s
                           #                   "HLT_Dimuon8_PsiPrime_Barrel",              # inclusive dimuon psi2s
                           #                   "HLT_DoubleMu4_PsiPrimeTrk_Displaced",      # displaced psi2s trk 
                           #                   "HLT_Dimuon0_Jpsi_Muon",                    # triple-mu (jpsi + muon)
                           #                   "HLT_Dimuon6_Jpsi_NoVertexing",             # jpsi no vertex very prescaled
                           #                   "HLT_Dimuon0er16_Jpsi_NoOS_NoVertexing",    # jpsi no vertexing paths
                           #                   "HLT_Dimuon0er16_Jpsi_NoVertexing"          # jpsi no vertexing paths
                           # ),

                           # interesting paths for 2017 and 2018
                           # paths      = cms.vstring(
                           #     "HLT_Dimuon25_Jpsi",                        # inclusive dimuon jpsi
                           #     "HLT_DoubleMu4_JpsiTrk_Displaced",          # displaced jpistrk or jpsi trktrk
                           #     "HLT_DoubleMu4_JpsiTrkTrk_Displaced",       # displaced jpistrk or jpsi trktrk
                           #     "HLT_Dimuon18_PsiPrime",                    # inclusive dimuon psi2s
                           #     "HLT_DoubleMu4_PsiPrimeTrk_Displaced"       # displaced psi2s trk 
                           # ),
                           
                           # All of Charmonium 2018 paths
                           # paths      = cms.vstring(
                           #     "HLT_Dimuon0_Jpsi3p5_Muon2",
                           #     "HLT_Dimuon0_Jpsi_L1_4R_0er1p5R",
                           #     "HLT_Dimuon0_Jpsi_L1_NoOS",
                           #     "HLT_Dimuon0_Jpsi_NoVertexing_L1_4R_0er1p5R",
                           #     "HLT_Dimuon0_Jpsi_NoVertexing_NoOS",
                           #     "HLT_Dimuon0_Jpsi_NoVertexing",
                           #     "HLT_Dimuon0_Jpsi",
                           #     "HLT_Dimuon0_LowMass_L1_0er1p5R",
                           #     "HLT_Dimuon0_LowMass_L1_0er1p5",
                           #     "HLT_Dimuon0_LowMass_L1_4R",
                           #     "HLT_Dimuon0_LowMass_L1_4",
                           #     "HLT_Dimuon0_LowMass",
                           #     "HLT_Dimuon10_PsiPrime_Barrel_Seagulls",
                           #     "HLT_Dimuon18_PsiPrime_noCorrL1",
                           #     "HLT_Dimuon18_PsiPrime",
                           #     "HLT_Dimuon20_Jpsi_Barrel_Seagulls",
                           #     "HLT_Dimuon25_Jpsi_noCorrL1",
                           #     "HLT_Dimuon25_Jpsi",
                           #     "HLT_DoubleMu2_Jpsi_DoubleTkMu0_Phi",
                           #     "HLT_DoubleMu2_Jpsi_DoubleTrk1_Phi1p05",
                           #     "HLT_DoubleMu4_3_Bs",
                           #     "HLT_DoubleMu4_3_Jpsi",
                           #     "HLT_DoubleMu4_JpsiTrkTrk_Displaced",
                           #     "HLT_DoubleMu4_JpsiTrk_Displaced",
                           #     "HLT_DoubleMu4_Jpsi_Displaced",
                           #     "HLT_DoubleMu4_Jpsi_NoVertexing",
                           #     "HLT_DoubleMu4_PsiPrimeTrk_Displaced",
                           #     "HLT_Mu30_TkMu0_Psi",
                           #     "HLT_Mu7p5_L2Mu2_Jpsi",
                           #     "HLT_Mu7p5_Track2_Jpsi",
                           #     "HLT_Mu7p5_Track3p5_Jpsi",
                           #     "HLT_Mu7p5_Track7_Jpsi",
                           # ),
                           
                           # BParkingNano paths (single mu)
                           paths      = cms.vstring(
                               "HLT_Mu7_IP4",
                               "HLT_Mu8_IP6",
                               "HLT_Mu8_IP5",
                               "HLT_Mu8_IP3",
                               "HLT_Mu8p5_IP3p5",
                               "HLT_Mu9_IP6",
                               "HLT_Mu9_IP5",
                               "HLT_Mu9_IP4",    
                               "HLT_Mu10p5_IP3p5",
                               "HLT_Mu12_IP6"
                           ),

)

trgTables = cms.Sequence(trgTable)



