#!/usr/bin/env python
"""Identify causes of false negatives based on Sanger validation data.
"""
import os

import vcf
import numpy
import pandas

base_dir = "/usr/local/projects/xprize/quaternary"
prep_dir = os.path.join(base_dir, "NA12878_fosmid/work/prep")
fn_vcf = os.path.join(base_dir, "NA12878_validate/work/NA12878-fosfinal-ceval2-ceval2-discordant.vcf")
fosmid_vcf = os.path.join(prep_dir, "NA12878-allfos-nomnp-fullcombine-wrefs-cleaned-annotated-ffilter-nosv-nofilter-cfilter.vcf")
nohap_vcf = os.path.join(prep_dir, "NA12878-allfos-nomnp-fullcombine-nocall-wrefs-samplefix-nonhaploid.vcf")
hap_vcf = os.path.join(prep_dir, "NA12878-allfos-nomnp-fullcombine-nocall-wrefs-samplefix.vcf")

def main():
    nohap_out = "NA12878-validate-nohap.vcf"
    filtered_out = "NA12878-validate-filtered.vcf"
    with open(fn_vcf) as in_handle:
        fn_coords = read_fn_coords(in_handle)
    with open(nohap_vcf) as in_handle:
        with open(nohap_out, "w") as out_handle:
            fn_coords = report_matches(fn_coords, in_handle, out_handle)
    with open(fosmid_vcf) as in_handle:
        with open(filtered_out, "w") as out_handle:
            fn_coords = report_matches(fn_coords, in_handle, out_handle)
    print "No explanation", fn_coords
    print "Incorrect haps"
    with open(nohap_out) as in_handle:
        print_hap_info(in_handle)
    print "All non-hap"
    with open(nohap_vcf) as in_handle:
        print_hap_info(in_handle)
    print "All hets"
    with open(hap_vcf) as in_handle:
        print_hap_info(in_handle)

def print_hap_info(in_handle):
    vals = {"isindel": [],
            "homref": [],
            "homvar": []}
    for rec in vcf.Reader(in_handle):
        hets = rec.get_hets()
        if hets:
            vals["isindel"].append(int(rec.is_indel))
            vals["homref"].append(hets[0].data.PL[0] / 10.0)
            vals["homvar"].append(hets[0].data.PL[-1] / 10.0)
    df = pandas.DataFrame(vals)
    print "Indels"
    print df[df["isindel"] == 1].describe()
    print "SNPs"
    print df[df["isindel"] == 0].describe()

def _get_coords(line):
    parts = line.split("\t")
    return (parts[0], parts[1])

def report_matches(coords, in_handle, out_handle):
    for line in in_handle:
        if line.startswith("#"):
            out_handle.write(line)
        else:
            cur_coords = _get_coords(line)
            if cur_coords in coords:
                out_handle.write(line)
                coords.remove(cur_coords)
    return coords

def read_fn_coords(in_handle):
    coords = []
    for line in (x for x in in_handle if not x.startswith("#")):
        coords.append(_get_coords(line))
    return set(coords)

if __name__ == "__main__":
    main()
