#!/usr/bin/env python
"""Explore classification metrics for validated true positives and negatives.

This identifies useful classification metrics to distinguish
true and false variants.
"""
import os
import shutil

import vcf
import pandas
import ramp
import sklearn

def main():
    tp_vcf = "NA12878-validate-tps.vcf"
    fp_vcf = "NA12878-validate-fps.vcf"
    out_decision = "decision-tree-%s.graphviz"
    #metrics = ['DP', 'Entropy', 'FS', 'GC', 'HRun', 'HaplotypeScore', 'MFE',
    #           'MQ', 'NBQ', 'ReadPosEndDist']
    metrics = ['FS', 'MFE',
               'NBQ', 'ReadPosEndDist']
    format_metrics = ["AD", "PL", "QUAL"]
    with open(tp_vcf) as in_handle:
        df_tp = read_vcf_metrics(in_handle, metrics, format_metrics, 1)
    with open(fp_vcf) as in_handle:
        df_fp = read_vcf_metrics(in_handle, metrics, format_metrics, -1)
    df = pandas.concat([df_tp, df_fp])
    for val, name in [(0, "snp"), (1, "indel")]:
        explore_ml_decisiontree(df[df["indel"] == val],
                                metrics + format_metrics, out_decision % name)
    print df_tp.describe()
    print df_fp.describe()

def explore_ml_decisiontree(df, metrics, out_decision):
    store_dir = os.path.join(os.getcwd(), "ramp")
    if os.path.exists(store_dir):
        shutil.rmtree(store_dir)
    os.makedirs(store_dir)
    context = ramp.DataContext(store=os.path.join(os.getcwd(), "ramp"),
                               data=df)
    config = ramp.Configuration(target="target", metrics=[ramp.metrics.AUC()])
    factory = ramp.ConfigFactory(config,
                                 features=[[ramp.BaseFeature(x) for x in metrics]],
                                 model=[sklearn.tree.DecisionTreeClassifier(max_depth=3)])
    for x in factory:
        ramp.models.fit(x, context)
        sklearn.tree.export_graphviz(x.model, out_file=out_decision,
                                     feature_names=metrics)

def read_vcf_metrics(in_handle, metrics, format_metrics, target):
    d = {"target": [],
         "indel": []}
    for x in metrics + format_metrics:
        d[x] = []
    for rec in vcf.VCFReader(in_handle):
        for x in metrics:
            d[x].append(rec.INFO.get(x, None))
        d["target"] = target
        d["AD"].append(_calc_ad(rec.samples[0].data))
        d["PL"].append(_calc_pl(rec.samples[0].data))
        d["QUAL"].append(rec.QUAL)
        d["indel"].append(int(rec.is_indel))
    return pandas.DataFrame(d)

def _calc_ad(data):
    if hasattr(data, "AD"):
        if data.GT == "0":
            want, other = data.AD
        elif data.GT == "1":
            other, want = data.AD
        else:
            raise ValueError
        return 1.0 - (want / float(want + other))

def _calc_pl(data):
    if hasattr(data, "PL"):
        if data.GT == "0":
            return data.PL[-1] / 10.0
        else:
            return data.PL[0] / 10.0

if __name__ == "__main__":
    main()
