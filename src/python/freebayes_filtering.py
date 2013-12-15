#!/usr/bin/env python
"""Explore approaches for improved filtering of FreeBayes variant calls.

The rough filters used (DP < 5 and QUAL < 20.0) leave ~3200 concordant variables
filtered while only removing ~850. This tries to identify approaches that can
improve on this.

The improved approach identifies an additional ~2100 concordant variants while
adding only ~400 additional discordants. It will hopefully also capture additional
discordants that slip through the current filter since it includes strand bias
calculations.

Usage:
  freebayes_filtering.py <concordant> <discordant>
"""
import collections
import os
import subprocess
import sys

import matplotlib.pyplot as plt
import pandas as pd
import ramp
import sklearn
import vcf

def main(concordant_vcf, discordant_vcf):
    with open(concordant_vcf) as in_handle:
        df_tp = vcf_to_stats(in_handle, "tp")
    with open(discordant_vcf) as in_handle:
        df_fp = vcf_to_stats(in_handle, "fp")
    df = pd.concat([df_tp, df_fp], keys=["tp", "fp"])
    for name, gt in [("hom", "1/1"), ("het", "0/1")]:
        cur_df = df[df["GT"] == gt]
        check_filtering(cur_df, hom_filter if name == "hom" else het_filter)
        explore_ml_decisiontree(cur_df, ["DP", "QR_QA", "QUAL"], 2, "freebayes-filter-tree-%s.graphviz" % name)
        plt.figure()
        cur_df.boxplot(by="target")
        plt.savefig("metrics-%s.pdf" % name)

        plt.figure()
        cur_df["QR_QA"].hist(by=cur_df["target"])
        plt.savefig("strandbias-%s.pdf" % name)

        plt.figure()
        cur_df["AD"].hist(by=cur_df["target"])
        plt.savefig("ad-%s.pdf" % name)

        if name == "hom":
            plot_dp_qual_box(cur_df, name)
        else:
            plot_qual_scatter(cur_df, name)

# ## Filters

def het_filter(x):
    return ((x["DP"] < 4) or
            (x["DP"] < 13 and x["QUAL"] < 20 and x["QR_QA"] > 10 and x["AD"] > 0.1))

def hom_filter(x):
    return ((x["DP"] < 4 and x["QUAL"] < 50) or
            (x["DP"] < 13 and (x["QR_QA"] > -90 or x["AD"] > 0.1)))

def check_filtering(df, fn):
    print fn
    for calltype in ["tp", "fp"]:
        print "before", calltype, len(df[df["target"] == calltype])
    f_df = df[df.apply(fn, axis=1)]
    for calltype in ["tp", "fp"]:
        print "after", calltype, len(f_df[f_df["target"] == calltype])

# ## Plotting

def plot_dp_qual_box(df, name):
    """Boxplot of filtering metrics for true/false positives, split by depth.
    """
    for calltype in ["tp", "fp"]:
        plt.figure()
        c_df = df[df["target"] == calltype]
        c_df.boxplot(by="DP")
        plt.savefig("box-%s-%s.pdf" % (name, calltype))

def plot_qual_scatter(cur_df, name):
    x, y = "QR_QA", "AD"
    cur_df = cur_df[cur_df["DP"] > 3]
    cur_df = cur_df[cur_df["QUAL"] < 20]
    plt.figure()
    for color, calltype in [("b", "tp"), ("r", "fp")]:
        c_df = cur_df[cur_df["target"] == calltype]
        plt.scatter(c_df[x], c_df[y], # s=c_df["DP"],
                    c=color, label=calltype)
    plt.legend()
    plt.savefig("scatter-%s.pdf" % (name))

# ## Machine learning

def explore_ml_decisiontree(df, metrics, depth, out_decision):
    """Try to identify a decision tree for differentiating TPs/FPs.
    """
    context = ramp.DataContext(data=df)
    config = ramp.Configuration(target="target", metrics=[ramp.metrics.AUC()])
    factory = ramp.ConfigFactory(config,
                                 features=[[ramp.BaseFeature(x) for x in metrics]],
                                 model=[sklearn.tree.DecisionTreeClassifier(max_depth=depth,
                                                                            criterion="gini")])
    for x in factory:
        ramp.models.fit(x, context)
        out_file = sklearn.tree.export_graphviz(x.model, out_file=out_decision,
                                                feature_names=metrics)
        out_file.close()
    out_pdf = "%s.pdf" % os.path.splitext(out_decision)[0]
    subprocess.check_call(["dot", "-T", "pdf", "-o", out_pdf, out_decision])

# ## Calculate metrics

def strand_bias(data):
    """Approach to assess strand bias by looking at reference/alt qualities
    https://groups.google.com/d/msg/freebayes/fX4TOAqXJrA/VTNf-xXKSB8J
    """
    altc = float(sum(data.AO) if isinstance(data.AO, list) else data.AO)
    refc = float(data.DP - altc)
    if refc > 0:
        qr = sum(data.QR) if isinstance(data.QR, list) else data.QR / float(refc)
    else:
        qr = 0
    if altc > 0:
        qa = sum(data.QA) if isinstance(data.QA, list) else data.QA / float(altc)
    else:
        qa = 0
    return (qr - qa ) / float(max([qr, qa])) * 100.0

def percent_ad_deviation(data):
    altc = float(sum(data.AO) if isinstance(data.AO, list) else data.AO)
    refc = float(data.DP - altc)
    expected = 0.5 if data.GT == "0/1" else 1.0
    return abs(expected - (float(altc) / float(altc + refc)))

def vcf_to_stats(in_handle, target):
    d = collections.defaultdict(list)
    for rec in vcf.VCFReader(in_handle):
        data = rec.samples[0].data
        d["target"] = target
        d["DP"].append(data.DP)
        d["QUAL"].append(rec.QUAL)
        d["GT"].append(data.GT)
        d["AD"].append(percent_ad_deviation(data))
        d["QR_QA"].append(strand_bias(data))
    return pd.DataFrame(d)

if __name__ == "__main__":
    main(*sys.argv[1:])
