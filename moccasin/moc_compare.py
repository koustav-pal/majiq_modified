import numpy as np
import pandas as pd
import time

from numba import cuda, float32, jit

def fitmodel_adjust_counts_log_linear_default(Y_input, modelmat_df, adjustment_modelmat_df, build_model_samples ):
    """Return adjusted counts (numpy array)

    For each junction, fit a linear model on the input (counts or PSI). Compute adjusted counts
    by removing the terms labeled as confounders and re-applying the additive
    residuals from the full model. Return adjusted counts.

    In the full_correction version: the init_modelmat should include only the
    confounders. This function adds an intercept to the confounders, and the
    intercept serves as the only non-confounder. Conceptually this is a more
    aggressive adjustment since it possibly throws away more variation from the
    main effects of interest (which are not modeled here, but were possibly used
    to derive RUV factors).
    """

    ts = time.time()

    Y = Y_input.copy()

    ##### PART 1: BUILD THE MODEL

    M = modelmat_df.values

    # Only build the model with the build_model_samples, which is already a boolean numpy array
    M_build = M[ build_model_samples, : ]
    Y_build = Y[ : , build_model_samples ]

    n_samples_build = Y_build.shape[1]

    controls_designated = int(n_samples_build < Y.shape[1]) # the build_model option was used, i.e. a subset of samples were designated as controls


    te = time.time()
    print(3, te-ts)


    #### OLS fit method
    # For each junction, exclude samples with missing values and try to model & adjust

    ts = time.time()
    Mterm = np.linalg.inv( M_build.T @ M_build ) @ M_build.T # default Mterm, to be used for juncs with no missing values

    did_adjust = np.zeros(shape=Y.shape, dtype=np.int8)

    adjusted_arr = Y.copy() # default new values to old values

    n_juncs = Y.shape[0]

    te = time.time()
    print(4, te-ts)

    ts = time.time()


    n_juncs_i = np.zeros((n_juncs,), dtype=np.float32)

    build_samples_i = np.ones((n_juncs, n_samples_build), dtype=np.float32)


    threadsperblock = 32
    blockspergrid = (n_juncs_i.size + (threadsperblock - 1)) // threadsperblock
    # test[blockspergrid, threadsperblock]

    fitmodel_kerm[blockspergrid, threadsperblock](n_juncs_i, Y, M, Y_build, M_build, [n_samples_build], adjustment_modelmat_df.values, adjusted_arr, Mterm, [controls_designated], build_samples_i)


    te = time.time()
    print(5, te-ts)

    return (adjusted_arr, did_adjust)

@cuda.jit
def fitmodel_kerm(n_juncs_i, Y, M, Y_build, M_build, n_samples_build, adjustment_modelmat_values, adjusted_arr, Mterm, controls_designated, build_samples_i):

    """
    :param i: iterator / mat position
    :param Y: np array, const
    :param M: np array, const
    :param Y_build: np array, const
    :param M_build: np array, const
    :param n_samples_build: int, const
    :param adjustment_modelmat_values: np array, const
    :param adjusted_arr:  np array, output writable
    :param Mterm: np array, const
    :param controls_designated: bool, const
    :return:
    """

    i = cuda.grid(1)

    n_samples_build = n_samples_build[0]
    controls_designated = controls_designated[0]

    ## Fit the model:
    build_samples_i_sum = 0

    #build_samples_i = Y_build[i, :] < 0 # check for missing values, represented by psi < 0
    for _i, v in enumerate(Y_build[i]):
        if v < 0:
            build_samples_i[i, _i] = 1
            build_samples_i_sum += 1
        # else:
        #     build_samples_i[i, _i] = 0

    control_has_missing_value = ( build_samples_i_sum < n_samples_build ) and controls_designated
    if control_has_missing_value:
        # Design decision: if the build_model option was used, i.e. control samples were designated, and a control has a missing value, then do not adjust the lsv
        # Log this here: did_adjust is already zeros by default, done and logged
        return

    Y_build_i = Y_build[ i , build_samples_i ]

    if build_samples_i_sum == n_samples_build:
        # No missing values for this junction! We already have this Mterm.
        Mterm_i = Mterm
    else:
        M_build_i = M_build[ build_samples_i, : ]
        try:
            Mterm_i = np.linalg.inv( M_build_i.T @ M_build_i ) @ M_build_i.T
        except np.linalg.LinAlgError:
            # Model matrix was no longer full-rank after removing samples with missing values. Cannot adjust this junction.
            # Log this here: did_adjust is already zeros by default, done and logged
            return

    params_i = Y_build_i @ Mterm_i.T
    # Done fitting model

    ## Adjust values:
    #n_samples_adjust = Y.shape[1]


    #adjust_samples_i = [Y[i, k] >= 0 for k in range(n_samples_adjust)] # check for missing values, represented by psi < 0
    # re-use build_samples_i instead of making new adjustname array.
    for _i, v in enumerate(Y[i]):
        if v >= 0:
            build_samples_i[i, _i] = 1
        else:
            build_samples_i[i, _i] = 0

    fit_predictions_i = ( M[ build_samples_i, : ] @ params_i.T ).T

    mu_no_confounders_arr_i = ( adjustment_modelmat_values[ build_samples_i, : ] @ params_i.T ).T

    resid_arr_i = Y[i, build_samples_i] - fit_predictions_i

    adjusted_arr[i, build_samples_i] = resid_arr_i + mu_no_confounders_arr_i  # update to new values for those samples and juncs which can be adjusted

    #did_adjust[i, :] += np.asarray([1 if v else 0 for v in adjust_samples_i], dtype=np.int8)

def LSV_sum_reexpand(df, LSVs):
    """Return a numpy array

    Assumes the input dataframe has an "LSV" column.
    Groupby this column and sum to aggregate, then re-expand to the original
    dimensions so that each LSV record gets its sum. Return as a numpy array.

    Utility function used to 1) compute total reads over junctions for each LSV,
    and 2) compute the sum of adjusted PSI in order to renormalize PSI for each LSV
    """

    df["LSV"] = LSVs

    perlsv_df = df.groupby("LSV", sort=False).sum()

    reexp_arr = perlsv_df.loc[ LSVs ].reset_index().drop("LSV", axis=1).values

    return reexp_arr

def fitmodel_adjustcounts_cij( coverage_df, linreg_modelmat_df, adjustment_modelmat_df, build_model_samples, LSVs ):

    ts = time.time()
    nreads_df = coverage_df.copy()

    target_cov_arr = LSV_sum_reexpand(nreads_df, LSVs)
    del nreads_df
    te = time.time()
    print(1, te-ts)

    ts = time.time()
    norm_arr = LSV_sum_reexpand(coverage_df.copy(), LSVs)
    norm_arr[norm_arr == 0] = 1  # specifically avoid 0/0
    norm_arr[norm_arr < 0] = 1  # ensure psi value is negative in the case of missing values, represented as -1 up to now
    psi_arr = coverage_df.values / norm_arr
    te = time.time()
    print(2, te-ts)


    (adjusted_counts_arr, did_adjust_counts) = fitmodel_adjust_counts_log_linear_default(psi_arr,
                                                                                 linreg_modelmat_df,
                                                                                 adjustment_modelmat_df,
                                                                                 build_model_samples)



    ts = time.time()
    ##### clip to 0 and re-normalize so sum_psi=1 per lsv in each sample
    adjusted_counts_arr = np.clip(adjusted_counts_arr, a_min=0, a_max=None)

    te = time.time()
    print(6, te-ts)

    ## ENDELSE

    ts = time.time()

    # Convert adjusted counts to psi and then get convert psi to the target coverage
    norm_arr = LSV_sum_reexpand(pd.DataFrame( adjusted_counts_arr ), LSVs)
    norm_arr[ norm_arr==0 ] = 1 # specifically avoid 0/0
    adjusted_psi_arr = adjusted_counts_arr / norm_arr

    ##### Convert psi to the target coverage
    coverage_arr = adjusted_psi_arr * target_cov_arr

    te = time.time()
    print(7, te-ts)

    return coverage_arr


@cuda.jit
def matmul(A, B, C):
    """Perform square matrix multiplication of C = A * B
    """
    i, j = cuda.grid(2)
    if i < C.shape[0] and j < C.shape[1]:
        tmp = 0.
        for k in range(A.shape[1]):
            tmp += A[i, k] * B[k, j]
        C[i, j] = tmp

@cuda.jit
def test(mat1, mat2, mat3):
    """Perform square matrix multiplication of C = A * B
    """
    i = cuda.grid(1)



    # for j, v in enumerate(mat3):
    #     if v < 0.5:
    #         mat2[i, j] = 0
    mat2 = np.linalg.inv(mat3)


import jax.numpy as jnp
from jax import grad, jit, vmap

@jit
def selu(x, alpha=1.67, lmbda=1.05):
    return lmbda * jnp.where(x > 0, x, alpha * jnp.exp(x) - alpha)

def selu_nj(x, alpha=1.67, lmbda=1.05):
    return lmbda * np.where(x > 0, x, alpha * np.exp(x) - alpha)

t_s = time.time()
for i in [x / 1000000.0 for x in range(0, 1000000)]:
    selu(i).block_until_ready()
t_e = time.time()

print(t_e - t_s)
assert False


if __name__ == "__main__":

    mat1 = np.ones((10, ), dtype=np.float64)

    mat2 = np.ones((10, 10), dtype=np.float64)
    mat3 = np.ones((10, 10), dtype=np.float64)

    mat3[2] = 0.3

    #mat2[0, mat3 < 0.5] = 0


    #mat3 = np.zeros((10, ), dtype=np.float32)

    threadsperblock = 32
    blockspergrid = (mat1.size + (threadsperblock - 1)) // threadsperblock
    test[blockspergrid, threadsperblock](mat1, mat2, mat3)

    print(mat1)
    print(mat2)
    assert False


    if True:
        coverage_df = pd.read_pickle("/home/pjewell/PycharmProjects/moccasin/testing/coverage_df.pkl")
        linreg_modelmat_df = pd.read_pickle("/home/pjewell/PycharmProjects/moccasin/testing/linreg_modelmat_df.pkl")
        adjustment_modelmat_df = pd.read_pickle("/home/pjewell/PycharmProjects/moccasin/testing/adjustment_modelmat_df.pkl")
        build_model_samples = np.load("/home/pjewell/PycharmProjects/moccasin/testing/build_model_samples.npy")
        LSVs = np.load("/home/pjewell/PycharmProjects/moccasin/testing/LSVs.npy")


        output_coverage_arr_old = np.load("/home/pjewell/PycharmProjects/moccasin/testing/output_coverage_arr.npy")

        start_t = time.time()
        print("Calculation start~")
        output_coverage_arr = fitmodel_adjustcounts_cij(coverage_df, linreg_modelmat_df, adjustment_modelmat_df, build_model_samples, LSVs)
        end_t = time.time()
        print(f"Calculation end, time was: {(end_t - start_t)} s")


        print(f"Equiv1? {np.array_equiv(output_coverage_arr, output_coverage_arr_old)}")






    # mat1 = np.ones((10, 10), dtype=np.float32)
    # mat2 = np.ones((10, 10), dtype=np.float32)
    # mat2[2:6, 3:8] = 2
    # mat3 = np.empty((10, 10), dtype=np.float32)
    #
    # threadsperblock = 32
    # blockspergrid = (mat1.size + (threadsperblock - 1)) // threadsperblock
    # matmul[blockspergrid, threadsperblock](mat1, mat2, mat3)

