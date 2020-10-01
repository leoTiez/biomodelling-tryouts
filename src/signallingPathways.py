#!/usr/bin/python3
# ##################################################
# Implementation of the investigative framework of gene and signalling pathways after the paper
# Kholodenko, Boris N., et al.
# "Untangling the wires: a strategy to trace functional interactions in signaling and gene networks."
# Proceedings of the National Academy of Sciences 99.20 (2002): 12841-12846.
# ##################################################

import numpy as np
from data.dataObjects import DUMMY_DATA_1, DUMMY_DATA_2


def calc_local_resp(global_resp):
    """
    Calculate the local response matrix as a function of the global response matrix
    :param global_resp: Matrix with the global response to perturbations of particular species
    :return: Local response matrix. The j-th row shows the effect on the species i, meaning the chnange
    in activation brought about by the changed activity of module j.
    """
    inv_glob_resp = np.linalg.inv(global_resp)
    local_resp = -1./np.diag(inv_glob_resp) * inv_glob_resp
    return local_resp


def create_global_resp_vec(before_vec, after_mat):
    """
    Row i represents the perturbation of the species, ie. the after matrix's element ij
    is the perturbation of the i-th species and j-th species' corresponding response
    :param before_vec: vector (or even matrix) with the stable state concentrations before the perturbation
    :param after_mat: matrix with the concentrations after the respective perturbations
    :return: global response
    """
    global_resp = 2 * (after_mat / before_vec - 1) / (after_mat / before_vec + 1)
    return global_resp


def main():
    loc_resp_trial1 = calc_local_resp(DUMMY_DATA_1)
    loc_resp_trial2 = calc_local_resp(DUMMY_DATA_2)
    loc_resp_mean = np.mean([loc_resp_trial1, loc_resp_trial2], axis=0)
    print(loc_resp_mean)


if __name__ == '__main__':
    main()
