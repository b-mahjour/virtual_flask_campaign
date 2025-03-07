import copy
import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors
import psycopg2
from virtual_flask_campaign.shared.reaction_class import returnReactionTemplates, VirtualFlask
from collections import Counter
from itertools import product
from pprint import pprint
import numpy as np
import matplotlib.pyplot as plt
from virtual_flask_campaign.shared.filters3 import apply_filters_local
from virtual_flask_campaign.shared.util import wells
import concurrent.futures
import pickle
from virtual_flask_campaign.shared.hpc.commands import precalculate_novelty_askcos


def connect_to_rds(
    host="10.63.5.43",
    port=5432,
    dbname="postgres",
    user="postgres",
    password="pass",
):
    """Create a connection to the PostgreSQL database"""
    try:
        conn = psycopg2.connect(
            host=host, port=port, dbname=dbname, user=user, password=password
        )
        # print("Connection established")
        return conn
    except Exception as e:
        print(f"Unable to connect to the database: {e}")
        print("No connection to RDS")
        return None


def get_mass_hits(cur, username, campaign_name, experiment_name, sample_name):
    cur.execute(
        "SELECT raw FROM campaign_analysis_packets WHERE username=%s AND campaign_name = %s AND experiment_name = %s AND method = %s AND sample_name = %s;",
        (
            username,
            campaign_name,
            experiment_name,
            "waters_lcms_quick_run_v1",
            sample_name,
        ),
    )
    hits = cur.fetchall()
    if len(hits) == 0:
        return [], [], []
    if len(hits) > 1:
        print("ERROR")
        return [], [], []
    hits = hits[0][0]
    uv = get_uv_trace(hits)

    msp = []
    msn = []
    for ii in hits:
        if ii["scan"] == "positive":
            msp.append(ii)
        if ii["scan"] == "negative":
            msn.append(ii)

    return msp, msn, uv


def get_forward_predictor_hits(
    cur, username, campaign_name, experiment_name, sample_name
):
    cur.execute(
        "SELECT raw FROM campaign_analysis_packets WHERE username=%s AND campaign_name = %s AND experiment_name = %s AND method = %s AND sample_name = %s;",
        (
            username,
            campaign_name,
            experiment_name,
            "forward_predictor_products",
            sample_name,
        ),
    )
    hits = cur.fetchall()
    if len(hits) == 0:
        return []
    if len(hits) > 1:
        print("ERROR")
        return []
    hits = hits[0][0]
    return hits


def get_reactants(cur, username, campaign_name, experiment_name, sample_name):
    cur.execute(
        "SELECT variable_value FROM campaign_experiments WHERE username=%s AND campaign_name=%s AND experiment_name = %s AND sample_name = %s AND variable_name = %s;",
        (username, campaign_name, experiment_name, sample_name, "substrate"),
    )
    hits = cur.fetchall()
    return hits


def plot_hit(hit_package, uvt):
    fig, ax = plt.subplots(2, 1, figsize=(2, 2), dpi=300, sharex=True)

    # set default font size and family to 6 and arial
    plt.rcParams.update({"font.size": 6, "font.family": "arial"})
    ax[0].yaxis.get_offset_text().set_fontsize(6)
    ax[1].yaxis.get_offset_text().set_fontsize(6)
    # set scientific notation for y axis
    ax[0].ticklabel_format(style="sci", axis="y", scilimits=(0, 0))

    mins = 4
    num_scans = len(hit_package["mz_trace"])
    xtl = [round(ii * mins / num_scans, 3) for ii in range(num_scans)]

    ax[0].plot(xtl, hit_package["mz_trace"], linewidth=0.5)
    ax[0].set_title(
        "EIC mz: "
        + str(round(hit_package["product_exact_mass"], 2))
        + f", found: {hit_package['target_matched_mz']}",
        fontsize=6,
        fontfamily="arial",
    )
    # ax[1].set_title("TWC", fontsize=6, fontfamily="arial")

    num_scans = len(uvt)
    xtl = [round(ii * mins / num_scans, 3) for ii in range(num_scans)]

    ax[1].plot(xtl, uvt, linewidth=0.5)
    ax[1].vlines(
        hit_package["min_start_uv"] * mins / num_scans,
        0,
        2e6,
        color="red",
        alpha=0.5,
        linewidth=0.5,
    )
    ax[1].vlines(
        hit_package["min_end_uv"] * mins / num_scans,
        0,
        2e6,
        color="green",
        alpha=0.5,
        linewidth=0.5,
    )
    ax[1].fill_between(
        [xx * mins / num_scans for xx in hit_package["uv_times"]],
        hit_package["uv_intensities"],
        color="darkblue",
        alpha=0.9,
        zorder=2,
        # linewidth=0.5
    )
    # set all to fontsize=6, fontfamily='arial'
    # ax[0].set_xticklabels(xtl, fontfamily="arial")
    ax[0].set_xlim(0, 3.8)
    # num_scans = len(uvt)
    ax[1].set_xlim(0, 3.8)
    # xtl = [round(ii * mins/num_scans,1) for ii in ax[1].get_xticks()]
    # ax[1].set_xticklabels(xtl, fontfamily="arial")
    ax[0].set_ylabel("EIC Intensity", fontsize=6, fontfamily="arial")
    ax[1].set_ylabel("TWC Intensity", fontsize=6, fontfamily="arial")
    ax[1].set_xlabel("r.t. (min.)", fontsize=6, fontfamily="arial")
    ax[0].set_yticks([])
    ax[1].set_yticks([])

    for a in ax:
        a.tick_params(labelsize=6)
        # set scale font size
        a.xaxis.set_tick_params(labelsize=6)
        a.yaxis.set_tick_params(labelsize=6)


def plot_well_results(well_res, well_name, look_for_mass=None):
    well_hit_data = well_res[0]
    uv_data = well_res[1]
    for idx, hit_frame in enumerate(well_hit_data):
        if look_for_mass != None:
            # print(round(hit_frame["product_exact_mass"],0), look_for_mass)
            if look_for_mass != round(hit_frame["product_exact_mass"], 0):
                continue
        plot_hit(hit_frame, uv_data)
        plt.subplots_adjust(hspace=0)
        # plt.tight_layout()
        # plt.show()
        plt.suptitle(well_name, fontsize=6, fontfamily="arial", fontweight="bold")
        plt.savefig(
            f"compiled_data/{well_name}-{idx}.png",
            dpi=300,
            bbox_inches="tight",
            pad_inches=0.01,
        )
        plt.close()
        # break


def return_mask(trace_in):
    current_mask = 1
    in_group = False
    mask = []
    for ii in trace_in:
        if ii == 0 and in_group == True:
            current_mask = current_mask + 1
            in_group = False
            mask.append(0)
        elif ii == 0 and in_group == False:
            mask.append(0)
        elif ii > 0 and in_group == True:
            mask.append(current_mask)
        elif ii > 0 and in_group == False:
            in_group = True
            mask.append(current_mask)
    return mask


def get_grouped_ms_frames(eic_trace_in, ms_scans):
    m = return_mask(eic_trace_in)
    prev = 0
    grouped_frames = {}
    for idx, ii in enumerate(m):
        if ii > 0:
            if ii not in grouped_frames:
                grouped_frames[ii] = []
            if ii != prev:
                prev = ii
                grouped_frames[prev] = []
            grouped_frames[ii].append(ms_scans[idx])
    return grouped_frames


def get_uv_scan_range_and_integration(uv_in, start_time, end_time, shift):
    min_start = 1e6
    min_end = 1e6
    scans_out = []
    for jj in uv_in:
        time_diff = np.abs(jj["start_time"] - start_time)
        if time_diff < min_start:
            min_start = time_diff
            min_start_uv = jj["scan"] + shift
        time_diff = np.abs(jj["start_time"] - end_time)
        if time_diff < min_end:
            min_end = time_diff
            min_end_uv = jj["scan"] + shift

    times = []
    intensities = []
    uv_frames = []
    for jj in uv_in:
        if jj["scan"] >= min_start_uv and jj["scan"] <= min_end_uv:
            times.append(jj["start_time"])
            intensities.append(jj["adj_bpi"])
            scans_out.append(jj["scan"])
            uv_frames.append(jj)
    integration = np.trapz(intensities, times)
    return min_start_uv, min_end_uv, integration, scans_out, intensities, uv_frames


def get_uv_integration(uv_array_in, ms_frames_in):

    uv_shift_scans = 21
    min_start_uv, min_end_uv, integration, times, intensities, uv_frames = (
        get_uv_scan_range_and_integration(
            uv_array_in,
            ms_frames_in[0]["start_time"],
            ms_frames_in[-1]["start_time"],
            uv_shift_scans,
        )
    )
    return integration, min_start_uv, min_end_uv, times, intensities, uv_frames


def get_uv_trace(mzml):
    uvs = []

    mecn_slope = (5.04e5 - 3.05e4) / (99.9 - 5)
    mecn_percent = (99.9 - 5) / (4535 - 934)
    nn = 0
    for k in mzml:
        # print(k["scan"], k)
        if k["scan"] == "pda":
            if k["start_time"] < 0.778:
                k["adj_bpi"] = k["bpi"]
            if k["start_time"] > 0.778 and k["start_time"] < 3.779:
                p_mecn = mecn_percent * nn - 19.6407
                k["adj_bpi"] = k["bpi"] - mecn_slope * p_mecn
            if k["start_time"] > 3.779:
                p_mecn = 99.9 - mecn_percent * (nn - 4535)
                k["adj_bpi"] = k["bpi"] - mecn_slope * p_mecn
            nn = nn + 1
            uvs.append(k)
    return uvs


def get_exact_masses(fp_smiles_in):
    masses_out = []
    for k in fp_smiles_in:
        mol = Chem.MolFromSmiles(k, sanitize=False)
        Chem.SanitizeMol(
            mol,
            sanitizeOps=Chem.SanitizeFlags.SANITIZE_NONE
            ^ Chem.SanitizeFlags.SANITIZE_SYMMRINGS
            ^ Chem.SanitizeFlags.SANITIZE_SETAROMATICITY,
        )
        exact_mass = Descriptors.ExactMolWt(mol)
        masses_out.append(exact_mass)
    return masses_out


def get_all_sets_of_triplets(substrates_in):
    all_triplets_s = list(product(substrates_in, repeat=3))
    seen = []
    outs = []
    for j in all_triplets_s:
        cc = Counter(j)
        if cc not in seen:
            seen.append(cc)
            outs.append(j)
    return outs


def get_vf_hits(substrates, mechs, returnNetwork=False):
    state_network = VirtualFlask(mechs)
    state_network.charge(substrates, "CS(C)=O")
    state_network.run_until_done(
        iters=5, thresh=50000, ring_filter=False, precalc_prods=[], time_limit=30
    )
    apply_filters_local(state_network)
    out_hits = []
    for k in state_network.nodes:
        if state_network.nodes[k].structural_failure == False:
            m = Chem.MolFromSmiles(state_network.nodes[k].other_data["target_molecule"])
            for at in m.GetAtoms():
                at.SetAtomMapNum(0)
            mz = Descriptors.ExactMolWt(m)
            sm = Chem.MolToSmiles(m)
            out_hits.append((sm, mz))
    if returnNetwork:
        return out_hits, state_network
    return out_hits


def get_intensities_at_mass(mz_target, mzml_raw, tol=0.5):
    intensities = []
    for k in mzml_raw:
        found = False
        for idx, mz in enumerate(k["mz_ar"]):
            if abs(mz - mz_target) < tol:
                intensities.append(k["in_ar"][idx])
                found = True
                break
        if not found:
            intensities.append(0)
    return intensities


def check_for_mass_hit(mass, mass_hits, mass_tolerance=0.5):
    sf = []
    for k in mass_hits:
        if abs(k["mz"] - mass) <= mass_tolerance:
            sf.append(k)
    return sf


def get_vf_hits_pipeline(
    substrates_in, mspa_in, msna_in, mechs_in, uva, triplets=True, verbose=False
):
    if triplets:
        all_triplets = get_all_sets_of_triplets(substrates_in)
    else:
        all_triplets = [tuple(substrates_in)]
    all_vf_hits = []
    adduct_to_mass_difference = {"groups_p1": 1, "groups_p23": 23, "groups_m1": -1}
    for trp in all_triplets:
        h0 = get_vf_hits(trp, mechs_in)
        for h in h0:
            eic_p1 = get_intensities_at_mass(h[1] + 1, mspa_in)
            eic_p23 = get_intensities_at_mass(h[1] + 23, mspa_in)
            eic_m1 = get_intensities_at_mass(h[1] - 1, mspa_in)
            groups_p1 = get_grouped_ms_frames(eic_p1, mspa_in)
            groups_p23 = get_grouped_ms_frames(eic_p23, mspa_in)
            groups_m1 = get_grouped_ms_frames(eic_m1, msna_in)
            all_vf_hits.append(
                {
                    "reaction_smiles": ".".join(trp) + ">>" + h[0],
                    "product_smiles": h[0],
                    "product_exact_mass": h[1],
                    "reaction_substrates": trp,
                    "trace_p1": eic_p1,
                    "trace_p23": eic_p23,
                    "trace_m1": eic_m1,
                    "groups_p1": groups_p1,
                    "groups_p23": groups_p23,
                    "groups_m1": groups_m1,
                }
            )

    final_hits = []
    for vf_hits in all_vf_hits:
        for adduct_group in ["groups_p1", "groups_p23", "groups_m1"]:
            for hit_frame in vf_hits[adduct_group]:
                if len(vf_hits[adduct_group][hit_frame]) <= 5:
                    continue

                (
                    uv_integration,
                    min_start_uv,
                    min_end_uv,
                    times,
                    intensities,
                    uv_frames,
                ) = get_uv_integration(uva, vf_hits[adduct_group][hit_frame])

                ms_to_intensity = {}
                for iii in vf_hits[adduct_group][hit_frame]:
                    for idxmz, imx in enumerate(iii["in_ar"]):
                        rounded_mz = round(iii["mz_ar"][idxmz], 1)
                        matched_mz = next(
                            (
                                mzz
                                for mzz in ms_to_intensity
                                if abs(mzz - rounded_mz) <= 0.5
                            ),
                            None,
                        )
                        if not matched_mz:
                            matched_mz = rounded_mz
                            ms_to_intensity[matched_mz] = 0
                        ms_to_intensity[matched_mz] += imx

                # sorted_ms_to_intensity = {k: v for k, v in sorted(ms_to_intensity.items(), key=lambda item: item[0], reverse=True)}
                most_abundant_mz = max(ms_to_intensity, key=ms_to_intensity.get)
                target_matched_mz = None
                closest_match = 100
                for mzz in ms_to_intensity:
                    this_match = abs(
                        mzz
                        - (
                            vf_hits["product_exact_mass"]
                            + adduct_to_mass_difference[adduct_group]
                        )
                    )
                    if this_match < closest_match:
                        target_matched_mz = mzz
                        closest_match = this_match

                final_hits.append(
                    {
                        "product_smiles": vf_hits["product_smiles"],
                        "product_exact_mass": vf_hits["product_exact_mass"],
                        "target_matched_mz": target_matched_mz,
                        "intensity_ratio": ms_to_intensity[target_matched_mz]
                        / ms_to_intensity[most_abundant_mz],
                        "uv_integration": uv_integration,
                        "ms_integration": ms_to_intensity[target_matched_mz],
                        "most_abundant_mz": most_abundant_mz,
                        "most_abundant_integration": ms_to_intensity[most_abundant_mz],
                        "uv_frames": uv_frames,
                        "mz_frames": vf_hits[adduct_group][hit_frame],
                        "adduct_group": adduct_group,
                        "mz_time_start": vf_hits[adduct_group][hit_frame][0][
                            "start_time"
                        ],
                        "mz_time_end": vf_hits[adduct_group][hit_frame][-1][
                            "start_time"
                        ],
                        "mz_trace": vf_hits["trace_" + adduct_group.split("_")[1]],
                        "min_start_uv": min_start_uv,
                        "min_end_uv": min_end_uv,
                        "uv_times": times,
                        "uv_intensities": intensities,
                    }
                )
                if verbose:
                    to_print = {}
                    for atr in final_hits[-1]:
                        if atr in [
                            "uv_frames",
                            "mz_frames",
                            "mz_trace",
                            "uv_intensities",
                            "uv_times",
                        ]:
                            to_print[f"{atr} length"] = len(final_hits[-1][atr])
                            continue
                        to_print[atr] = final_hits[-1][atr]

                    pprint(to_print)
                    print()
    return final_hits


def run_extract_frames(exp_name):
    mechs = returnReactionTemplates()

    conn = connect_to_rds()
    cur = conn.cursor()

    all_hits = []
    # for exp_name in wells_in:
    print(exp_name)
    hits, uvt = _run_extract_frames(exp_name, mechs, cur)
    print(exp_name, len(hits), "hits")
    print()
    all_hits.append((hits, uvt))
    # print(len(all_hits))
    cur.close()
    conn.close()
    return all_hits


def _run_extract_frames(exp_name, mechs, cur):
    mspa, msna, uva = get_mass_hits(exp_name, cur)
    substrates = [ii[0]["smiles"] for ii in get_reactants(exp_name, cur)]

    hits = get_vf_hits_pipeline(substrates, mspa, msna, mechs, uva, True, verbose=False)
    with open(f"compiled_data/results_{exp_name}.pkl", "wb") as f:
        pickle.dump((hits), f)

    return hits, [uu["adj_bpi"] for uu in uva]




def run_fp(reactants):
    try:
        ps, scs = precalculate_novelty_askcos(reactants, model="at")
        ps2, scs2 = precalculate_novelty_askcos(reactants, model="g2s")
        ps3, scs3 = precalculate_novelty_askcos(reactants, model="wln")
    except:
        return [], []

    ps = ps + ps2 + ps3
    scs = scs + scs2 + scs3
    return ps, scs


def upload_packet(
    username,
    campaign_name,
    experiment_name,
    sample_name,
    variable_name,
    variable_value,
):

    conn = connect_to_rds()
    cur = conn.cursor()
    cur.execute(
        "INSERT INTO campaign_experiments (username, campaign_name, experiment_name, sample_name, variable_value, variable_name) VALUES (%s, %s, %s, %s, %s,%s)",
        (
            username,
            campaign_name,
            experiment_name,
            sample_name,
            variable_value,
            variable_name,
        ),
    )

    conn.commit()
    cur.close()
    conn.close()
    return True


def run_parallel(rwells):
    with concurrent.futures.ProcessPoolExecutor() as executor:
        results = list(executor.map(run_extract_frames, rwells))  # Run in parallel


if __name__ == "__main__":
    rwells = wells[96]
    # rwells = ["A1"]
    run_parallel(rwells)
