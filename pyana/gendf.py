import numpy as np
import pandas as pd
import uproot
import pickle
import os
import json
import argparse
from optparse import OptionParser
from tqdm import tqdm
from multiprocessing import Pool
from multiprocessing import get_context
from pandas_helpers import *
from branches import *

MASS_MUON = 0.105
MASS_PROTON = 0.938272
MASS_NEUTRON = 0.9395654 
MASS_PIZERO = 0.1349768 
MASS_PIPLUS = 0.13957039 
MASS_PIMINUS = 0.13957039 

def make_hdrdf(events):
    # hdr_df = events.arrays(hdr_branches, library="pd")
    hdrs = loadbranches(events, hdr_branches)
    _, ind = np.unique(hdrs.rec.hdr.subrun + hdrs.rec.hdr.run*100, return_index=True)
    pot = np.sum(hdrs.rec.hdr.pot[ind])
    return hdrs, pot


def make_nudf(events):
    nus = loadbranches(events, mc_branches)
    nus = nus.rec.mc.nu
    return nus


def make_nuprimdf(events):
    nuprims = loadbranches(events, mc_prim_branches)
    nuprims = nuprims.rec.mc.nu.prim
    return nuprims


def make_slcdf(events):
    slcs = loadbranches(events, slc_branches)
    slcs = slcs.rec.slc

    slcs_mc = loadbranches(events, slc_mc_branches)
    slcs_mc = slcs_mc.rec.slc.truth
    slcs_mc.rename(mapper = lambda x: f'mc_{x}', axis="columns", level=0, inplace=True)
    slcs = slcs.join(slcs_mc)

    slcs_mcprims = loadbranches(events, slc_mc_prim_branches)
    slcs_mcprims = slcs_mcprims.rec.slc.truth.prim
    slcs["mc_max_proton_KE"] = ((slcs_mcprims.genE - MASS_PROTON) * (slcs_mcprims.pdg == 2212)).groupby(level=[0,1]).max()
    slcs["mc_sum_proton_KE"] = ((slcs_mcprims.genE - MASS_PROTON) * (slcs_mcprims.pdg == 2212)).groupby(level=[0,1]).sum()
    return slcs


def make_pfpdf(events):
    pfptrks = loadbranches(events, pfp_trk_branches)
    pfptrks = pfptrks.rec.slc.reco.pfp.trk

    pfps = loadbranches(events, pfp_branches)
    pfps = pfps.rec.slc.reco #.pfp
    pfptrks = pfptrks.join(pfps)

    chi2s = loadbranches(events, pfp_trk_chi2_branches)
    chi2s = chi2s.rec.slc.reco.pfp.trk.chi2pid
    chi2s.rename(mapper = lambda x: f'chi2_{x}', axis="columns", level=0, inplace=True)
    pfptrks = pfptrks.join(chi2s)

    pfp_trk_mc_df = loadbranches(events, pfp_trk_mc_branches)
    pfp_trk_mc_df = pfp_trk_mc_df.rec.slc.reco.pfp.trk.truth.p
    pfp_trk_mc_df.rename(mapper = lambda x: f'mc_{x}', axis="columns", level=0, inplace=True)
    pfptrks = pfptrks.join(pfp_trk_mc_df)
    return pfptrks


def make_hitdf(events):
    hits = loadbranches(events, trkhitbranches)
    hits = hits.rec.slc.reco.pfp.trk.calo.I2.points
    return hits


def make_dfs(make_dfs_args):
    which_dfs, files, fidx = make_dfs_args
    if len(which_dfs) == 0:
        print("No dataframes requested")
        exit(1)

    events = uproot.open(files[fidx]+":recTree")

    keys = []
    dfs = []
    lastidx = []

    hdrs, pot = make_hdrdf(events)    
    hdrs["fidx"] = fidx
    hdrs.set_index("fidx", append=True, inplace=True)
    keys.append("pot")
    dfs.append(pot)
    keys.append("hdr")
    dfs.append(hdrs)

    if "nu" in which_dfs:
        nus = make_nudf(events)
        nus["fidx"] = fidx
        nus.set_index("fidx", append=True, inplace=True)
        keys.append("nu")
        dfs.append(nus)
        lastidx.append(nus.index.levels[0][-1])

    if "nuprim" in which_dfs:
        nuprims = make_nuprimdf(events)
        nuprims["fidx"] = fidx
        nuprims.set_index("fidx", append=True, inplace=True)
        keys.append("nuprim")
        dfs.append(nuprims)
        lastidx.append(nuprims.index.levels[0][-1])

    if "slc" in which_dfs:
        slcs = make_slcdf(events)
        slcs["fidx"] = fidx
        slcs.set_index("fidx", append=True, inplace=True)
        keys.append("slc")
        dfs.append(slcs)
        lastidx.append(slcs.index.levels[0][-1])
    
    if "pfp" in which_dfs:
        pfps = make_pfpdf(events)
        pfps["fidx"] = fidx
        pfps.set_index("fidx", append=True, inplace=True)
        keys.append("pfp")
        dfs.append(pfps)
        lastidx.append(pfps.index.levels[0][-1])

    if "hit" in which_dfs:
        hits = make_hitdf(events)
        hits["fidx"] = fidx
        hits.set_index("fidx", append=True, inplace=True)
        keys.append("hit")
        dfs.append(hits)
        lastidx.append(hits.index.levels[0][-1])

    return keys, dfs, (fidx, lastidx[0])


def concat_dfs(dfs_list, shift_list):
    dfs_list_shifted = []
    for df in dfs_list:
        this_fidx = df.reset_index().fidx.unique()[0]
        if this_fidx == 0: # don't shift the idxs of the first df
            entry_shift = 0
        else:
            entry_shift = shift_list[this_fidx-1]+this_fidx
        idx_names = df.index.names
        df = df.reset_index()
        df["entry"] = df["entry"] + entry_shift
        dfs_list_shifted.append(df)
    df = pd.concat(dfs_list_shifted)
    df.set_index(idx_names[:-1], inplace=True)
    df.sort_index(inplace=True)
    # df = df.sort_index(level=[0])
    return df


def main(nprocs, files, which_dfs, savepath, savename):
    dfs_list = []
    lastidx_dict = {}
    with Pool(processes=nprocs) as pool:
        make_dfs_args = [(which_dfs, files, fidx) for fidx in range(len(files))]
        for ret in tqdm(pool.imap_unordered(make_dfs, make_dfs_args)):
            if ret is not None:
                keys = ret[0] 
                dfs_list.append(ret[1])
                fidx, lastidx = ret[2]
                lastidx_dict[fidx] = lastidx
    dfs_list = {k: [dfs_list[f][idx] for f in range(len(files))] for idx, k in enumerate(keys)}
    lastidxs = sorted(lastidx_dict.items())
    shift_list = np.cumsum([lastidxs[i][1] for i in range(len(lastidxs))])
    dfs = {k: concat_dfs(dfs_list[k], shift_list) for k in keys[1:]}
    dfs["pot"] = np.sum(dfs_list["pot"])

    with open(os.path.join(savepath, savename), "wb") as f:
        pickle.dump(dfs, f)
    print("dfs saved to ", os.path.join(savepath, savename))


def get_args():
    parser = argparse.ArgumentParser(description='make dfs from flat caf files')
    parser.add_argument('-c', '--config', help='JSON with script configuration')
    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = get_args()
    with open(args.config, 'r') as fconf:
        config = json.loads(fconf.read());

    CPU_COUNT = os.cpu_count()
    nprocs = max(CPU_COUNT//2, config["ncpu"])
    print("using ", nprocs, " processes")

    with open(config["filelist"]) as f:
        files = f.readlines()
        files = [x.strip() for x in files]
    print("processing ", len(files), " files")

    which_dfs = config["which_dfs"]
    print("producing ", which_dfs, " dfs")
    savepath = config["savepath"]
    savename = config["savename"]

    main(nprocs, files, which_dfs, savepath, savename)