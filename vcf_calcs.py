import numpy as np
import os
import argparse
import pysam


# helper func
def calc_hm_ht(vcf_files, n):
    homozygous_counts = np.zeros(n)
    heterozygous_counts = np.zeros(n)

    sample_names = []
    for i in range(n):
        sample_names.append("tsk_" + str(i))

    for x in vcf_files:
        vcf = pysam.VariantFile(x)
        reads = vcf.fetch()
        for v in reads:
            for i in range(n):
                gt = v.samples[sample_names[i]]["GT"]
                if gt[0] == gt[1] and gt[0] != 0 and gt[0] != '.':
                    homozygous_counts[i] += 1
                elif (gt[0] == 1 or gt[1] == 1) and gt[0] != gt[1]:
                    heterozygous_counts[i] += 1
                else:
                    continue

    num_f = len(vcf_files)
    homozygous_counts /= num_f
    heterozygous_counts /= num_f

    return homozygous_counts, heterozygous_counts


# calcs number of heterozygous or homozygous sites
def calc_individual_pgs(direc, n):
    vcf_files = []
    for x in os.listdir(direc):
        if x.endswith(".vcf"):
            vcf_files.append(x)


    onf_nm = [f for f in vcf_files if f.startswith("ONF__NM")]
    onf_nm_hm, onf_nm_ht = calc_hm_ht(onf_nm, n)
    np.save("homo_cts_onf_nm", onf_nm_hm)
    np.save("hetero_cts_onf_nm", onf_nm_ht)

    onf_ud = [f for f in vcf_files if f.startswith("ONF__UD")]
    onf_ud_hm, onf_ud_ht = calc_hm_ht(onf_ud, n)
    np.save("homo_cts_onf_ud", onf_ud_hm)
    np.save("hetero_cts_onf_ud", onf_ud_ht)

    onf_rc = [f for f in vcf_files if f.startswith("ONF__RC")]
    onf_rc_hm, onf_rc_ht = calc_hm_ht(onf_rc, n)
    np.save("homo_cts_onf_rc", onf_rc_hm)
    np.save("hetero_cts_onf_rc", onf_rc_ht)

    onf_dm = [f for f in vcf_files if f.startswith("ONF__DM")]
    onf_dm_hm, onf_dm_ht = calc_hm_ht(onf_dm, n)
    np.save("homo_cts_onf_dm", onf_dm_hm)
    np.save("hetero_cts_onf_dm", onf_dm_ht)

    abs_nm = [f for f in vcf_files if f.startswith("ABS__NM")]
    abs_nm_hm, abs_nm_ht = calc_hm_ht(abs_nm, n)
    np.save("homo_cts_abs_nm", abs_nm_hm)
    np.save("hetero_cts_abs_nm", abs_nm_ht)

    abs_ud = [f for f in vcf_files if f.startswith("ABS__UD")]
    abs_ud_hm, abs_ud_ht = calc_hm_ht(abs_ud, n)
    np.save("homo_cts_abs_ud", abs_ud_hm)
    np.save("hetero_cts_abs_ud", abs_ud_ht)

    abs_rc = [f for f in vcf_files if f.startswith("ABS__RC")]
    abs_rc_hm, abs_rc_ht = calc_hm_ht(abs_rc, n)
    np.save("homo_cts_abs_rc", abs_rc_hm)
    np.save("hetero_cts_abs_rc", abs_rc_ht)

    abs_dm = [f for f in vcf_files if f.startswith("ABS__DM")]
    abs_dm_hm, abs_dm_ht = calc_hm_ht(abs_dm, n)
    np.save("homo_cts_abs_dm", abs_dm_hm)
    np.save("hetero_cts_abs_dm", abs_dm_ht)

    cst_nm = [f for f in vcf_files if f.startswith("CST__NM")]
    cst_nm_hm, cst_nm_ht = calc_hm_ht(cst_nm, n)
    np.save("homo_cts_cst_nm", cst_nm_hm)
    np.save("hetero_cts_cst_nm", cst_nm_ht)

    cst_ud = [f for f in vcf_files if f.startswith("CST__UD")]
    cst_ud_hm, cst_ud_ht = calc_hm_ht(cst_ud, n)
    np.save("homo_cts_cst_ud", cst_ud_hm)
    np.save("hetero_cts_cst_ud", cst_ud_ht)

    cst_rc = [f for f in vcf_files if f.startswith("CST__RC")]
    cst_rc_hm, cst_rc_ht = calc_hm_ht(cst_rc, n)
    np.save("homo_cts_cst_rc", cst_rc_hm)
    np.save("hetero_cts_cst_rc", cst_rc_ht)

    cst_dm = [f for f in vcf_files if f.startswith("CST__DM")]
    cst_dm_hm, cst_dm_ht = calc_hm_ht(cst_dm, n)
    np.save("homo_cts_cst_dm", cst_dm_hm)
    np.save("hetero_cts_cst_dm", cst_dm_ht)


def main():
    parser = argparse.ArgumentParser(
        description="Compute number of Homozygous and Heterozygous sites in Simulations.")

    parser.add_argument("-dir", action="store", dest="dir", type=str,
                        default="/Users/mohamedabdelrahman/Desktop")
    parser.add_argument("-n", action="store", dest="n", type=int,
                        default=50)

    args = parser.parse_args()
    direc = args.dir
    n = args.n

    calc_individual_pgs(direc, n)

if __name__ == '__main__':
    main()