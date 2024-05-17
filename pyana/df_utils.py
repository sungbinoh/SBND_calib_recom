import numpy as np
import pandas as pd
from utils import InFV, mag

def merge_dfs(dfs):
    # load dfs
    pot = dfs["pot"]
    nu_df = dfs["nu"]
    nuprim_df = dfs["nuprim"]
    slc_df = dfs["slc"]
    pfp_df = dfs["pfp"]

    nu_df.rename(mapper = lambda x: f'nu_{x}', axis='columns', level=0, inplace=True)
    nuprim_df["totp"] = np.sqrt((nuprim_df.startp.x)**2+(nuprim_df.startp.y)**2+(nuprim_df.startp.z)**2) # |momentum| branch

    # slice truth matching
    slc_df.loc[np.invert(slc_df[("tmatch","eff")] > 0.5) & (slc_df[("tmatch","idx")] >= 0), ("tmatch","idx")] = np.nan
    slc_df.rename(mapper=lambda x: f'slc_{x}', axis='columns', level=0, inplace=True)
    match_df = pd.merge(slc_df.reset_index(), 
                        nu_df.reset_index(),
                        left_on=[("entry", ""), ("slc_tmatch", "idx")], 
                        right_on=[("entry", ""), ("rec.mc.nu..index", "")], 
                        how="left", # Keep all slices
                        ) 
    match_df = match_df.set_index(["entry", "rec.slc..index"], verify_integrity=True)
    assert(match_df.reset_index().groupby(["entry", "rec.mc.nu..index"])["slc_charge"].count().max() == 1)

    # calculate some track level stuff
    pfp_df["mc_totp"] = mag(pfp_df.mc_startp.x, pfp_df.mc_startp.y, pfp_df.mc_startp.z)
    pfp_df["is_contained"] = InFV(pfp_df.start) & InFV(pfp_df.end)
    pfp_df.loc[(pfp_df["is_contained"]), ("totp_muon", "")]  = pfp_df.loc[(pfp_df["is_contained"]), ("rangeP", "p_muon")]  #.rangeP.p_muon
    pfp_df.loc[np.invert(pfp_df["is_contained"]), ("totp_muon", "")] = pfp_df.loc[np.invert(pfp_df["is_contained"]), ("mcsP", "fwdP_muon")] #.mcsP.fwdP_muon
    pfp_df.loc[(pfp_df["is_contained"]), ("totp_proton", "")]  = pfp_df.loc[(pfp_df["is_contained"]), ("rangeP", "p_proton")] 
    pfp_df.loc[np.invert(pfp_df["is_contained"]), ("totp_proton", "")] = pfp_df.loc[np.invert(pfp_df["is_contained"]), ("mcsP", "fwdP_proton")]

    all_df = pd.merge(match_df.reset_index(), 
                    pfp_df.reset_index(),
                    left_on=[("entry", ""), ("rec.slc..index", "")], 
                    right_on=[("entry", ""), ("rec.slc..index", "")], 
                    how="right", # Keep every track
                    )
    all_df = all_df.set_index(["entry", "rec.slc..index", "rec.slc.reco.pfp..index"], verify_integrity=True)

    return all_df