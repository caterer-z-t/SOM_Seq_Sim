import pandas as pd
import numpy as np
import statsmodels.api as sm
import statsmodels.formula.api as smf


def glmm_interact(
    dataset,
    cluster,
    contrast1,
    contrast2,
    random_effects=None,
    fixed_effects=None,
    verbose=False,
    save_models=False,
    save_model_dir=None,
    save_name=None,
):
    # Convert cluster to string
    cluster = dataset[cluster].astype(str)
    dataset = dataset.copy()

    # Generate design matrix from cluster assignments
    designmat = pd.get_dummies(cluster, prefix="cluster", drop_first=False)
    dataset = pd.concat([designmat, dataset], axis=1)

    # Create output list to hold results
    cluster_names = designmat.columns
    cluster_models = {}

    # Create model formulas
    if fixed_effects is not None and random_effects is not None:
        model_rhs = " + ".join(
            [*fixed_effects, f"(1|{random_effects})", contrast1, contrast2]
        )
    elif fixed_effects is not None:
        model_rhs = " + ".join([*fixed_effects, contrast1, contrast2])
        if verbose:
            print("No random effects specified")
            raise ValueError("No random effects specified")
    elif random_effects is not None:
        model_rhs = " + ".join([f"(1|{random_effects})", contrast1, contrast2])
    else:
        model_rhs = " + ".join([contrast1, contrast2])
        if verbose:
            print("No random or fixed effects specified")
            raise ValueError("No random or fixed effects specified")

    if verbose:
        print(f"Using full model: cluster ~ {model_rhs}")

    # Run nested mixed-effects models for each cluster
    for test_cluster in cluster_names:
        if verbose:
            print(f"Creating logistic mixed models for {test_cluster}")

        null_formula = f"{test_cluster} ~ 1 + {model_rhs}"
        full_formula = f"{test_cluster} ~ {contrast1}:{contrast2} + {model_rhs}"

        # Run null and full mixed-effects models
        null_model = smf.mixedlm(
            null_formula, dataset, groups=dataset[random_effects]
        ).fit()
        full_model = smf.mixedlm(
            full_formula, dataset, groups=dataset[random_effects]
        ).fit()

        # Likelihood Ratio Test
        model_lrt = sm.stats.anova_lm(null_model, full_model)

        # Calculate confidence intervals for the contrast term
        contrast_lvl2 = f"{contrast1}_{dataset[contrast1].cat.categories[1]}:{contrast2}_{dataset[contrast2].cat.categories[1]}"
        contrast_ci = full_model.conf_int()  # Adjust as necessary

        # Save model objects to list
        cluster_models[test_cluster] = {
            "null_model": null_model,
            "full_model": full_model,
            "model_lrt": model_lrt,
            "confint": contrast_ci,
        }

    # Organize results into output dataframe
    output = pd.DataFrame({"cluster": cluster_names, "size": designmat.sum().values})
    output["model.pvalue"] = [
        model["model_lrt"].iloc[1]["PR(>F)"] for model in cluster_models.values()
    ]
    output[f"{contrast_lvl2}.OR"] = [
        np.exp(model["full_model"].params[contrast_lvl2])
        for model in cluster_models.values()
    ]
    output[f"{contrast_lvl2}.OR.95pct.ci.lower"] = [
        np.exp(model["confint"].loc[contrast_lvl2][0])
        for model in cluster_models.values()
    ]
    output[f"{contrast_lvl2}.OR.95pct.ci.upper"] = [
        np.exp(model["confint"].loc[contrast_lvl2][1])
        for model in cluster_models.values()
    ]

    return output
