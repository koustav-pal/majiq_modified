
from moccasin import MoccasinModel
import numpy as np
import pandas as pd

def load_base_lsv_order(path):
    """
    Find Key (lsv_id, start, end) to position-order-index for the file we will align everything to
    (usually just the first file)

    """

    np_arr = np.load(path)['junc_info']

    keys = {}
    for i, row in enumerate(np_arr):
        keys[(row[0], row[1], row[2])] = i

    return keys


def load_grouping_from_majiq(path):
    np_arr = np.load(path)['junc_info']
    pd_df = pd.DataFrame(np_arr)
    grouping = pd_df.groupby("f0", sort=False).head(1).drop(["f0", "f1", "f2", "f3", "f4"], axis=1).index.to_numpy()

    return grouping

def load_coverage(path):
    npf = np.load(path)
    base_coverage = npf['coverage']
    return base_coverage

def load_coverage_with_order(path, base_order):
    npf = np.load(path)
    base_coverage = npf['coverage']
    base_info = npf['junc_info']
    final_coverage = np.zeros(base_coverage.shape)

    mapping = []
    for i, row in enumerate(base_info):
        key = (row[0], row[1], row[2],)
        mapping.append(base_order[key])

    final_coverage[mapping, :] = base_coverage[:]


    return final_coverage, mapping

def revert_coverage_with_order(coverage, mapping):
    # using the mapping obtained earlier, revert the order back to how it was originally
    return coverage[mapping, :]



# read from majiq files to get coverage

base_order = load_base_lsv_order('/home/pjewell/scf_test_cases/4majiqfile/workshop_Cer3.majiq')
grouping = load_grouping_from_majiq('/home/pjewell/scf_test_cases/4majiqfile/workshop_Cer3.majiq')


workshop_Cer3_f_coverage, workshop_Cer3_f_order = load_coverage_with_order('/home/pjewell/scf_test_cases/4majiqfile/workshop_Cer3.majiq', base_order)
workshop_Adr2_f_coverage, workshop_Adr2_f_order = load_coverage_with_order('/home/pjewell/scf_test_cases/4majiqfile/workshop_Adr2.majiq', base_order)


model_defs = {'workshop_Cer3': {"intercept": 1, "da_batch": 0},
              'workshop_Adr2': {"intercept": 1, "da_batch": 1}}

model = MoccasinModel(model_defs, confounding_factors={'da_batch': 0}, max_num_processes=12)

output1 = model.model_with_bootstraps(grouping, {
    'workshop_Cer3': workshop_Cer3_f_coverage,
    'workshop_Adr2': workshop_Adr2_f_coverage,
              })

output1['workshop_Adr2'] = revert_coverage_with_order(output1['workshop_Adr2'], workshop_Adr2_f_order)
output1['workshop_Cer3'] = revert_coverage_with_order(output1['workshop_Cer3'], workshop_Cer3_f_order)

# with just one bootstrap test, as 1d array

output2 = {
    'workshop_Cer3' : np.zeros(workshop_Cer3_f_coverage.shape),
    'workshop_Adr2' : np.zeros(workshop_Adr2_f_coverage.shape),
}

for i in range(workshop_Cer3_f_coverage.shape[1]):
    output2_tmp = model.model(grouping, {
        'workshop_Cer3': workshop_Cer3_f_coverage[:, i],
        'workshop_Adr2': workshop_Adr2_f_coverage[:, i],
                  })
    output2['workshop_Cer3'][:, i] = output2_tmp['workshop_Cer3']
    output2['workshop_Adr2'][:, i] = output2_tmp['workshop_Adr2']

output2['workshop_Adr2'] = revert_coverage_with_order(output2['workshop_Adr2'], workshop_Adr2_f_order)
output2['workshop_Cer3'] = revert_coverage_with_order(output2['workshop_Cer3'], workshop_Cer3_f_order)


# compare to default moccasin adjustment (test case)
moccasin_cer3_coverage = np.load("/home/pjewell/scf_test_cases/4majiqfile/output_quantifiable_lsvs_squash/workshop_Cer3.scf_adjusted.majiq")['coverage']
moccasin_adr2_coverage = np.load("/home/pjewell/scf_test_cases/4majiqfile/output_quantifiable_lsvs_squash/workshop_Adr2.scf_adjusted.majiq")['coverage']

print(np.allclose(output1['workshop_Adr2'], moccasin_adr2_coverage))
print(np.allclose(output1['workshop_Cer3'], moccasin_cer3_coverage))

print(np.allclose(output2['workshop_Adr2'], moccasin_adr2_coverage))
print(np.allclose(output2['workshop_Cer3'], moccasin_cer3_coverage))

# print("final_out", output['workshop_Adr2'][:20, 0])
# print("final_out m", moccasin_cer3_coverage[:20, 0])


