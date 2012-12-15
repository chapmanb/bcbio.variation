#!/usr/bin/env python
"""Explore classification metrics for validated true positives and negatives.

This identifies useful classification metrics to distinguish
true and false variants.
"""
import os
import sys
import math
import shutil
import itertools

import yaml
import vcf
import numpy as np
import pandas
import ramp
from ramp.estimators.sk import BinaryProbabilities
import sklearn

base_dir = "/usr/local/projects/xprize/quaternary"
prep_dir = os.path.join(base_dir, "NA12878_fosmid/work/prep")
mp_tps_vcf = os.path.join(prep_dir, "NA12878-allfos-nomnp-fullcombine-wrefs-cleaned-annotated-ffilter-nosv-nofilter-tps.vcf")
mp_fps_vcf = os.path.join(prep_dir, "NA12878-allfos-nomnp-fullcombine-wrefs-cleaned-annotated-ffilter-nosv-nofilter-fps.vcf")

default_config = {"sanger": {"tp": "NA12878-validate-tps.vcf",
                             "fp": "NA12878-validate-fps.vcf"},
                  "full": {"tp": mp_tps_vcf,
                           "fp": mp_fps_vcf}}

def main(config_file = None):
    if config_file:
        with open(config_file) as in_handle:
            config = yaml.load(in_handle)
    else:
        config = default_config
    metrics = ['Entropy', 'FS', 'GC', 'HRun', 'HaplotypeScore', 'MFE',
               'MQ', 'NBQ', 'ReadPosEndDist']
    format_metrics = ["AD", "PL", "QUAL", "DP"]
    #sanger_decisiontree(config["sanger"]["tp"], config["sanger"]["fp"], metrics, format_metrics)
    ml_params(config["full"]["tp"], config["full"]["fp"], metrics, format_metrics)

# ## Machine learning parameters

def ml_params(tp_vcf, fp_vcf, metrics, format_metrics):
    """Explore machine learning parameters to identify approaches to help separate data.
    """
    metrics = ['Entropy', 'FS', 'MFE',
               'MQ', 'NBQ', 'ReadPosEndDist']
    exploring = False
    with open(tp_vcf) as in_handle:
        df_tp = read_vcf_metrics(in_handle, metrics, format_metrics, 1,
                                 exploring)
    with open(fp_vcf) as in_handle:
        df_fp = read_vcf_metrics(in_handle, metrics, format_metrics, -1,
                                 exploring)
    df = pandas.concat([df_tp, df_fp], keys=["tp", "fp"])
    df = df.fillna({"NBQ": df["NBQ"].mean(), "PL" : df["PL"].mean(),
                    "AD" : df["AD"].mean(), "FS": 0.0, "DP": df["DP"].mean()})
    df = normalize_inputs(df, metrics + format_metrics)
    for val, name in [(0, "snp"), (1, "indel")]:
        print "--->", name
        linear_metric_explore(df[df["indel"] == val], metrics + format_metrics)
        #ml_param_explore(df[df["indel"] == val], metrics + format_metrics,
        #                 exploring)

def normalize_inputs(df, metrics):
    """Normalize all inputs around mean and standard deviation.
    """
    for m in metrics:
        mean = np.mean(df[m])
        stdev = np.std(df[m])
        def std_normalize(x):
            return (x - mean) / stdev
        #df[m] = df[m].map(std_normalize)
        xmin = min(df[m])
        xmax = max(df[m])
        def minmax_normalize(x):
            return (x - xmin) / (xmax - xmin)
        df[m] = df[m].map(minmax_normalize)
    return df

def ml_param_explore(df, metrics, test_all=False):
    """Explore classification approaches and parameters, leveraging ramp to compare multiple models.
    """
    print df.describe()
    context = ramp.DataContext(data=df)
    config = ramp.Configuration(target="target", metrics=[ramp.metrics.AUC(), ramp.metrics.F1(),
                                                          ramp.metrics.HingeLoss()])
    #rf_params = _find_rf_params(df, metrics)
    models = [sklearn.ensemble.RandomForestClassifier(n_estimators=50,
                                                      max_features=int(math.ceil(math.sqrt(len(metrics))))),
              sklearn.linear_model.LogisticRegression()]
    if test_all:
        svm_tester = "linear"
        if svm_tester == "linear":
            linear_params = _find_svm_rbf_params(df, metrics, "linear")
            models.append(sklearn.svm.SVC(kernel="linear", C=linear_params["C"]))
        else:
            rbf_params = _find_svm_rbf_params(df, metrics, "rbf")
            models.append(sklearn.svm.SVC(kernel="rbf", C=rbf_params["C"], gamma=rbf_params["gamma"]))

    factory = ramp.ConfigFactory(config,
                                 features=[[ramp.BaseFeature(x) for x in metrics]],
                                 model = models)
    for x in factory:
        ramp.models.cv(x, context, folds=5, repeat=2,
                       print_results=True)

def _get_metrics_groups(metrics):
    # keeps: AD -- snps, MQ -- indels
    to_remove = [['Entropy', 'MFE', 'NBQ'],
                 ['DP', 'Entropy', 'MFE', 'NBQ'],
                 ['PL', 'QUAL'],
                 ['Entropy', 'MFE', 'NBQ', 'PL', 'QUAL'],
                 ['DP', 'Entropy', 'MFE', 'NBQ', 'PL', 'QUAL'],
                 ['DP', 'Entropy', 'MFE', 'NBQ', 'QUAL'],
                 ['DP', 'MFE', 'NBQ', 'QUAL'],
                 ['MFE', 'NBQ', 'QUAL'],
                 ]
    for rem in to_remove:
        new = metrics[:]
        for x in rem:
            new.remove(x)
        yield new

def linear_metric_explore(df, metrics):
    """Explore different combinations of metrics with a linear classifier.
    """
    print df.describe()
    context = ramp.DataContext(data=df)
    config = ramp.Configuration(target="target", metrics=[ramp.metrics.AUC()])
    models = [sklearn.svm.SVC(kernel="linear", C=100.0)]
    for sub_metrics in [metrics] + list(_get_metrics_groups(metrics)):
        print "==>", sub_metrics
        factory = ramp.ConfigFactory(config, model=models,
                                     features=[[ramp.BaseFeature(x) for x in sub_metrics]])
        for x in factory:
            ramp.models.cv(x, context, folds=5, repeat=2,
                           print_results=True)

def _find_rf_params(df, metrics):
    """Perform a grid search to find best parameters for random forest.
    """
    context = ramp.DataContext(data=df)
    config = ramp.Configuration(target="target",
                                features=[ramp.BaseFeature(x) for x in metrics])
    x, y = ramp.models.get_xy(config, context)

    n = len(metrics)
    param_grid = dict(max_features=range(int(math.ceil(math.sqrt(n))), n+1, 3),
                      n_estimators=range(20, 101, 20))
    grid = sklearn.grid_search.GridSearchCV(sklearn.ensemble.RandomForestClassifier(),
                                            param_grid=param_grid,
                                            cv=sklearn.cross_validation.StratifiedKFold(y=y, k=3))
    grid.fit(x, y)
    print grid.best_estimator_
    out = {}
    for attr in param_grid.keys():
        out[attr] = getattr(grid.best_estimator_, attr)
    return out

def _find_svm_rbf_params(df, metrics, kernel):
    """Perform a grid search to find best parameters for a SVM RBF kernel.
    """
    context = ramp.DataContext(data=df)
    config = ramp.Configuration(target="target",
                                features=[ramp.BaseFeature(x) for x in metrics])
    x, y = ramp.models.get_xy(config, context)

    if kernel == "linear":
        param_grid = dict(C=10.0 ** np.arange(-2, 5))
    else:
        param_grid = dict(gamma=10.0 ** np.arange(-5, 4),
                          C=10.0 ** np.arange(-2, 9))
    grid = sklearn.grid_search.GridSearchCV(sklearn.svm.SVC(kernel=kernel),
                                            param_grid=param_grid,
                                            cv=sklearn.cross_validation.StratifiedKFold(y=y, k=3),
                                            verbose=True)
    grid.fit(x, y)
    print grid.best_estimator_
    out = {}
    for attr in param_grid.keys():
        out[attr] = getattr(grid.best_estimator_, attr)
    return out

# ## Decision tree visualization

def sanger_decisiontree(tp_vcf, fp_vcf, metrics, format_metrics):
    """Explore sanger true/false calls with a decision tree.

    This helps identify primary metrics helping to discriminate the data.
    """
    out_decision = "decision-tree-%s.graphviz"
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
    context = ramp.DataContext(data=df)
    config = ramp.Configuration(target="target", metrics=[ramp.metrics.AUC()])
    factory = ramp.ConfigFactory(config,
                                 features=[[ramp.BaseFeature(x) for x in metrics]],
                                 model=[sklearn.tree.DecisionTreeClassifier(max_depth=3)])
    for x in factory:
        ramp.models.fit(x, context)
        sklearn.tree.export_graphviz(x.model, out_file=out_decision,
                                     feature_names=metrics)

# ## Parse VCFs into pandas data frames

def read_vcf_metrics(in_handle, metrics, format_metrics, target,
                     use_subset = False):
    d = {"target": [],
         "indel": []}
    for x in metrics + format_metrics:
        d[x] = []
    if use_subset:
        recs = itertools.islice(vcf.VCFReader(in_handle), 10000)
    else:
        recs = vcf.VCFReader(in_handle)
    for rec in recs:
        for x in metrics:
            d[x].append(rec.INFO.get(x, None))
        d["target"] = target
        d["AD"].append(_calc_ad(rec.samples[0].data))
        d["PL"].append(_calc_pl(rec.samples[0].data))
        format_dp = _calc_dp(rec.samples[0].data)
        d["DP"].append(format_dp)
        d["QUAL"].append(rec.QUAL)
        d["indel"].append(int(rec.is_indel))
    return pandas.DataFrame(d)

def _calc_dp(data):
    if hasattr(data, "DP"):
        return data.DP

def _calc_ad(data):
    if hasattr(data, "AD"):
        if data.GT in ["0", "0/0"]:
            target = 1.0
            want, other = data.AD
        elif data.GT in ["1", "1/1"]:
            target = 1.0
            other, want = data.AD
        elif data.GT in ["0/1"]:
            target = 0.5
            want, other = data.AD
        else:
            print "Need to handle", data.AD, data.GT
            return None
            raise ValueError
        if want + other > 0:
            return target - (want / float(want + other))

def _calc_pl(data):
    if hasattr(data, "PL"):
        if data.GT == "0":
            return data.PL[-1] / 10.0
        elif data.GT == "1":
            return data.PL[0] / 10.0
        else:
            return min(x for x in data.PL if x > 0) / 10.0

if __name__ == "__main__":
    main(*sys.argv[1:])
