#!/usr/bin/env python
"""Extract true positives and negatives with metrics based on validation results.

Also identifies incorrect true and false positives in training data sets to pull
out areas where the approach to creating training sets could improve.
"""
import os

base_dir = "/usr/local/projects/xprize/quaternary"
val_vcf = os.path.join(base_dir, "NA12878_validate/inputs/NA12878-sanger-validate-allowone.vcf")
prep_dir = os.path.join(base_dir, "NA12878_fosmid/work/prep")
full_vcf = os.path.join(prep_dir, "NA12878-allfos-nomnp-fullcombine-wrefs-cleaned-annotated.vcf")
mp_tps_vcf = os.path.join(prep_dir, "NA12878-allfos-nomnp-fullcombine-wrefs-cleaned-annotated-ffilter-nosv-nofilter-tps.vcf")
mp_fps_vcf = os.path.join(prep_dir, "NA12878-allfos-nomnp-fullcombine-wrefs-cleaned-annotated-ffilter-nosv-nofilter-fps.vcf")

def main():
    base = "NA12878-validate-%s.vcf"
    tp_out = base % "tps"
    fp_out = base % "fps"
    with open(val_vcf) as in_handle:
        tp_coords, fp_coords = read_validation(in_handle)
        coords = {"fp": fp_coords, "tp": tp_coords}
    with open(tp_out, "w") as tp_handle:
        with open(fp_out, "w") as fp_handle:
            with open(full_vcf) as in_handle:
                outs = {"fp": fp_handle, "tp": tp_handle}
                write_vcf_subsets(in_handle, outs, coords)
    wrongtp_out = base % "train-wrongtps"
    wrongfp_out = base % "train-wrongfps"
    with open(wrongtp_out, "w") as out_handle:
        with open(mp_tps_vcf) as in_handle:
            write_incorrect_train(in_handle, out_handle, coords["fp"])
    with open(wrongfp_out, "w") as out_handle:
        with open(mp_fps_vcf) as in_handle:
            write_incorrect_train(in_handle, out_handle, coords["tp"])

def _get_coords(line):
    parts = line.split("\t")
    return (parts[0], parts[1])

def write_incorrect_train(in_handle, out_handle, coords):
    for line in in_handle:
        if line.startswith("#"):
            out_handle.write(line)
        else:
            if _get_coords(line) in coords:
                out_handle.write(line)

def write_vcf_subsets(in_handle, outs, coords):
    for line in in_handle:
        if line.startswith("#"):
            for h in outs.itervalues():
                h.write(line)
        else:
            cur_coords = _get_coords(line)
            for k, v in coords.iteritems():
                if cur_coords in v:
                    outs[k].write(line)

def read_validation(in_handle):
    tp_coords = []
    fp_coords = []
    for line in (x for x in in_handle if not x.startswith("#")):
        if line.rstrip().endswith("0"):
            fp_coords.append(_get_coords(line))
        else:
            tp_coords.append(_get_coords(line))
    return set(tp_coords), set(fp_coords)

if __name__ == "__main__":
    main()
