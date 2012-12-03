#!/usr/bin/env python
"""Explore classification metrics for validated true positives and negatives.

This identifies useful classification metrics to distinguish
true and false variants.
"""
import os
import math
import shutil
import itertools

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

def main():
    tp_vcf = "NA12878-validate-tps.vcf"
    fp_vcf = "NA12878-validate-fps.vcf"
    metrics = ['DP', 'Entropy', 'FS', 'GC', 'HRun', 'HaplotypeScore', 'MFE',
               'MQ', 'NBQ', 'ReadPosEndDist']
    format_metrics = ["AD", "PL", "QUAL"]
    sanger_decisiontree(tp_vcf, fp_vcf, metrics, format_metrics)
    ml_params(mp_tps_vcf, mp_fps_vcf, metrics, format_metrics)

# ## Machine learning parameters

def ml_params(tp_vcf, fp_vcf, metrics, format_metrics):
    """Explore machine learning parameters to identify approaches to help separate data.
    """
    metrics = ['Entropy', 'FS', 'HRun', 'MFE',
               'MQ', 'NBQ', 'ReadPosEndDist']
    exploring = True
    with open(tp_vcf) as in_handle:
        df_tp = read_vcf_metrics(in_handle, metrics, format_metrics, 1,
                                 exploring)
    with open(fp_vcf) as in_handle:
        df_fp = read_vcf_metrics(in_handle, metrics, format_metrics, -1,
                                 exploring)
    df = pandas.concat([df_tp, df_fp], keys=["tp", "fp"])
    if exploring:
        df = df.fillna({"NBQ": df["NBQ"].mean(), "PL" : df["PL"].mean(),
                        "AD" : df["AD"].mean(), "FS": 0.0})
    for val, name in [(0, "snp"), (1, "indel")]:
        print "--->", name
        ml_param_explore(df[df["indel"] == val], metrics + format_metrics,
                         exploring)

def ml_param_explore(df, metrics, test_all=False):
    """Explore classification approaches and parameters, leveraging ramp to compare multiple models.
    """
    print df.describe()
    context = ramp.DataContext(data=df)
    config = ramp.Configuration(target="target", metrics=[ramp.metrics.AUC(), ramp.metrics.F1(),
                                                          ramp.metrics.HingeLoss()])
    #rf_params = _find_rf_params(df, metrics)
    models = [sklearn.ensemble.RandomForestClassifier(n_estimators=50,
                                                      max_features=int(math.ceil(math.sqrt(len(metrics)))))]
    if test_all:
        rbf_params = _find_svm_rbf_params(df, metrics)
        models += [sklearn.svm.SVC(kernel="linear"),
                   sklearn.linear_model.LogisticRegression(),
                   sklearn.svm.SVC(kernel="rbf", C=rbf_params["C"], gamma=rbf_params["gamma"])]
    factory = ramp.ConfigFactory(config,
                                 features=[[ramp.BaseFeature(x) for x in metrics]],
                                 model = models)
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

def _find_svm_rbf_params(df, metrics):
    """Perform a grid search to find best parameters for a SVM RBF kernel.
    """
    context = ramp.DataContext(data=df)
    config = ramp.Configuration(target="target",
                                features=[ramp.BaseFeature(x) for x in metrics])
    x, y = ramp.models.get_xy(config, context)

    param_grid = dict(gamma=10.0 ** np.arange(-5, 4),
                      C=10.0 ** np.arange(-2, 9))
    grid = sklearn.grid_search.GridSearchCV(sklearn.svm.SVC(kernel="rbf"),
                                            param_grid=param_grid,
                                            cv=sklearn.cross_validation.StratifiedKFold(y=y, k=3))
    grid.fit(x, y)
    print grid.best_estimator_
    return {"C": grid.best_estimator_.C, "gamma": grid.best_estimator_.gamma}

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
        recs = itertools.islice(vcf.VCFReader(in_handle), 5000)
    else:
        recs = vcf.VCFReader(in_handle)
    for rec in recs:
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
