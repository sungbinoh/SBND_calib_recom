pfpbranch = "rec.slc.reco.pfp."
trkbranch = pfpbranch + "trk."

trkhitbranches_perplane = lambda IPLANE : [
    trkbranch + "calo.%i.points.dedx"% IPLANE,
    trkbranch + "calo.%i.points.dqdx"% IPLANE,
    trkbranch + "calo.%i.points.pitch"% IPLANE,
    trkbranch + "calo.%i.points.integral"% IPLANE,
    trkbranch + "calo.%i.points.rr"% IPLANE,
    trkbranch + "calo.%i.points.wire"% IPLANE,
    trkbranch + "calo.%i.points.tpc"% IPLANE,
    trkbranch + "calo.%i.points.sumadc"% IPLANE,
    trkbranch + "calo.%i.points.t"% IPLANE,
    trkbranch + "calo.%i.points.x"% IPLANE,
    trkbranch + "calo.%i.points.y"% IPLANE,
    trkbranch + "calo.%i.points.z"% IPLANE,

    #trkbranch + "calo.%i.points.width"% IPLANE,
    #trkbranch + "calo.%i.points.mult"% IPLANE,
    #trkbranch + "calo.%i.points.tdc0"% IPLANE,

    trkbranch + "calo.%i.points.truth.h_e"% IPLANE,
    trkbranch + "calo.%i.points.truth.h_nelec"% IPLANE,
    trkbranch + "calo.%i.points.truth.pitch"% IPLANE,
    trkbranch + "calo.%i.points.truth.rr"% IPLANE,
]

trkhitbranches = trkhitbranches_perplane(2)
trkhitbranches_P1 = trkhitbranches_perplane(1)
trkhitbranches_P0 = trkhitbranches_perplane(0)

hdr_branches = ["rec.hdr.run", "rec.hdr.subrun", "rec.hdr.pot"]

mc_branches = [
    "rec.mc.nu.E",
    "rec.mc.nu.position.x",
    "rec.mc.nu.position.y",
    "rec.mc.nu.position.z",
    "rec.mc.nu.pdg",
    "rec.mc.nu.iscc",
    "rec.mc.nu.genie_inttype",
    "rec.mc.nu.genie_mode",
    "rec.mc.nu.momentum.x",
    "rec.mc.nu.momentum.y",
    "rec.mc.nu.momentum.z",
    'rec.mc.nu.nneutron', 'rec.mc.nu.npiminus', 'rec.mc.nu.npiplus', 'rec.mc.nu.npizero', 'rec.mc.nu.nproton', 
    'rec.mc.nu.nprim',
    'rec.mc.nu.genweight',
    'rec.mc.nu.time',
    'rec.mc.nu.parent_dcy_mode',
    'rec.mc.nu.parent_pdg',
    'rec.mc.nu.Q2', 'rec.mc.nu.q0', 'rec.mc.nu.q0_lab', 'rec.mc.nu.t', 'rec.mc.nu.w', 'rec.mc.nu.bjorkenX', 'rec.mc.nu.inelasticityY']


mc_prim_branches = ['rec.mc.nu.prim.pdg', 'rec.mc.nu.prim.genE',
    'rec.mc.nu.prim.start.x', 'rec.mc.nu.prim.start.y', 'rec.mc.nu.prim.start.z',
    'rec.mc.nu.prim.end.x','rec.mc.nu.prim.end.y','rec.mc.nu.prim.end.z',
    'rec.mc.nu.prim.startp.x','rec.mc.nu.prim.endp.x',
    'rec.mc.nu.prim.startp.y','rec.mc.nu.prim.endp.y',
    'rec.mc.nu.prim.startp.z', 'rec.mc.nu.prim.endp.z',
    'rec.mc.nu.prim.startE',
    'rec.mc.nu.prim.end_process', 'rec.mc.nu.prim.start_process'                 
                   ]

slc_branches = [
    "rec.slc.self",
    "rec.slc.charge",
    "rec.slc.truth.pdg",
    "rec.slc.producer",
    "rec.slc.nu_score",
    "rec.slc.vertex.x", "rec.slc.vertex.y", "rec.slc.vertex.z",
    "rec.slc.tmatch.eff", "rec.slc.tmatch.pur", "rec.slc.tmatch.index",
    "rec.slc.fmatch.score", "rec.slc.fmatch.time",
		]

slc_mc_branches = [
    "rec.slc.truth.E",
    "rec.slc.truth.position.x",
    "rec.slc.truth.position.y",
    "rec.slc.truth.position.z",
    "rec.slc.truth.pdg",
    "rec.slc.truth.iscc",
    "rec.slc.truth.genie_mode"
]

slc_mc_prim_branches = [
    "rec.slc.truth.prim.genE",
    "rec.slc.truth.prim.pdg",
]

pfp_trk_branches = [
    "rec.slc.reco.pfp.trk.start.x", "rec.slc.reco.pfp.trk.start.y", "rec.slc.reco.pfp.trk.start.z",
    "rec.slc.reco.pfp.trk.end.x", "rec.slc.reco.pfp.trk.end.y", "rec.slc.reco.pfp.trk.end.z",
    "rec.slc.reco.pfp.trk.dir.x", "rec.slc.reco.pfp.trk.dir.y", "rec.slc.reco.pfp.trk.dir.z",
    "rec.slc.reco.pfp.trk.phi", "rec.slc.reco.pfp.trk.costh",
    "rec.slc.reco.pfp.trk.len",
    "rec.slc.reco.pfp.trk.rangeP.p_muon",
    "rec.slc.reco.pfp.trk.mcsP.fwdP_muon",
    "rec.slc.reco.pfp.trk.mcsP.bwdP_muon",
    "rec.slc.reco.pfp.trk.mcsP.is_bwd_muon",
    "rec.slc.reco.pfp.trk.rangeP.p_pion",
    "rec.slc.reco.pfp.trk.mcsP.fwdP_pion",
    "rec.slc.reco.pfp.trk.rangeP.p_proton",
    "rec.slc.reco.pfp.trk.mcsP.fwdP_proton",
    "rec.slc.reco.pfp.trk.bestplane",
#    "rec.slc.reco.pfp.trk.dazzle.bestScore",
#    "rec.slc.reco.pfp.trk.dazzle.muonScore",
#    "rec.slc.reco.pfp.trk.dazzle.otherScore",
#    "rec.slc.reco.pfp.trk.dazzle.pdg",
#    "rec.slc.reco.pfp.trk.dazzle.pionScore",
#    "rec.slc.reco.pfp.trk.dazzle.protonScore",
    ]

pfp_branches = ["rec.slc.reco.pfp.trackScore",
                "rec.slc.reco.pfp.ndaughters"]

pfp_trk_chi2_branches = [
    "rec.slc.reco.pfp.trk.chi2pid.2.chi2_kaon", "rec.slc.reco.pfp.trk.chi2pid.2.chi2_muon", "rec.slc.reco.pfp.trk.chi2pid.2.chi2_pion", "rec.slc.reco.pfp.trk.chi2pid.2.chi2_proton",
    "rec.slc.reco.pfp.trk.chi2pid.1.chi2_kaon", "rec.slc.reco.pfp.trk.chi2pid.1.chi2_muon", "rec.slc.reco.pfp.trk.chi2pid.1.chi2_pion", "rec.slc.reco.pfp.trk.chi2pid.1.chi2_proton",
    "rec.slc.reco.pfp.trk.chi2pid.0.chi2_kaon", "rec.slc.reco.pfp.trk.chi2pid.0.chi2_muon", "rec.slc.reco.pfp.trk.chi2pid.0.chi2_pion", "rec.slc.reco.pfp.trk.chi2pid.0.chi2_proton",
]

pfp_trk_mc_branches_names = [
    "interaction_id",
    "parent",
    "pdg",
    "G4ID",
    "end_process",
    "start_process",
    "startE",
    "start.x", "start.y", "start.z",
    "startp.x", "startp.y", "startp.z",
    "end.x", "end.y", "end.z",
    "endp.x", "endp.y", "endp.z",
    "genp.x", "genp.y", "genp.z",
    "genE",
    "length",
    "cont_tpc",
]
pfp_trk_mc_branches = ["rec.slc.reco.pfp.trk.truth.p." + n for n in pfp_trk_mc_branches_names]

stub_branches = [
    "rec.slc.reco.stub.vtx.x",
    "rec.slc.reco.stub.vtx.y",
    "rec.slc.reco.stub.vtx.z",
    "rec.slc.reco.stub.end.x",
    "rec.slc.reco.stub.end.y",
    "rec.slc.reco.stub.end.z",
    "rec.slc.reco.stub.efield_vtx",
    "rec.slc.reco.stub.efield_end",
    "rec.slc.reco.stub.pfpid",
    "rec.slc.reco.stub.truth.p.pdg",
    "rec.slc.reco.stub.truth.p.genE",
    "rec.slc.reco.stub.truth.p.interaction_id",
]

stub_plane_branches = [
    "rec.slc.reco.stub.planes.p",
    "rec.slc.reco.stub.planes.hit_w",
    "rec.slc.reco.stub.planes.vtx_w",
    "rec.slc.reco.stub.planes.pitch",
    "rec.slc.reco.stub.planes.trkpitch",
]

stub_hit_branches = [
    "rec.slc.reco.stub.planes.hits.charge",
    "rec.slc.reco.stub.planes.hits.ontrack",
    "rec.slc.reco.stub.planes.hits.wire",
]


# # shower branches
# trueshwbranches = ['rec.slc.reco.shw.truth.p.pdg', 'rec.slc.reco.shw.truth.p.startE']