#!/usr/bin/python3
# ##################################################
# Implementation of the investigative framework of the DNA repair pathway
# Luijsterburg, Martijn S., et al.
# "Stochastic and reversible assembly of a multiprotein DNA repair complex ensures
# accurate target site recognition and efficient repair."
# Journal of Cell Biology 189.3 (2010): 445-463.
# ##################################################
import numpy as np

from data.dataObjects import XPC, TFIIH, XPG, XPA, XPF, RPA, PCNA
from data.dataObjects import DAMAGED, PARTIALLY, FULLY, INCISED, RESYNTH, RECHROM


def binary_to_int(binary_str):
    return int(binary_str, 2)


def int_to_binary(integer, num_digits):
    binary = str(bin(integer))[2:]
    return str(0) * (num_digits - len(binary)) + binary


class Protein:
    def __init__(self, init_conc, on_dict, off_dict, name=""):
        self.conc = init_conc
        self.on = on_dict
        self.off = off_dict
        self.name = name

    def protein_conc_update(self, protein_update):
        self.conc += protein_update


class Intermediates:
    def __init__(self, enzyme_dict, mask, num_proteins=7, name=""):
        self.conc = np.zeros(2**num_proteins)
        self.enzyme_dict = enzyme_dict
        self.mask = np.asarray(mask)
        self.name = name
        self.trans_inter = {}

    def set_transition_inter(self, trans_inter_dict):
        self.trans_inter = trans_inter_dict

    def intermediate_conc_update(self, proteins, dt=1e-4):
        """
        Update of the intermediate concentrations using the Euler approximation
        """
        protein_bound_on = np.zeros(len(proteins))
        protein_bound_off = np.zeros(len(proteins))
        protein_unbound_on = np.zeros(len(proteins))
        protein_unbound_off = np.zeros(len(proteins))
        for num, (protein, m) in enumerate(zip(proteins, self.mask)):
            if not m:
                continue

            protein_bound_on[num] = protein.on[self.name] * protein.conc
            protein_unbound_on[num] = -protein_bound_on[num]
            protein_unbound_off[num] = protein.off[self.name]
            protein_bound_off[num] = -protein_unbound_off[num]

        protein_off = [protein_unbound_off, protein_bound_off]
        protein_on = [protein_unbound_on, protein_bound_on]
        protein_change = np.zeros(len(proteins))
        for num_state in range(self.conc.size):
            binary = int_to_binary(num_state, len(proteins))
            if np.sum(np.asarray(list(binary)).astype("int")[~self.mask]) > 0:
                continue
            enzyme = 0
            if binary in self.enzyme_dict.keys():
                for inter, factor in zip(self.trans_inter[binary], self.enzyme_dict[binary]["factor"]):
                    enzyme += inter.conc[binary_to_int(binary)] * factor
            for num_protein, b in enumerate(binary):
                b = int(b)
                self.conc[num_state] += (protein_on[b][num_protein] * self.conc[num_state]**(1 - b)
                                         + protein_off[b][num_protein] * self.conc[num_state]**b
                                         + enzyme) * dt
                if b == 0:
                    protein_change[num_protein] -= proteins[num_protein].on[self.name] * proteins[num_protein].conc \
                                                   * self.conc[num_state]
                else:
                    protein_change[num_protein] += proteins[num_protein].off[self.name] * self.conc[num_state]

        return protein_change * dt


class RepairPathway:
    def __init__(self, proteins, intermediates):
        self.proteins = proteins
        self.intermediates = intermediates

    def run(self, num_steps=1):
        for _ in range(num_steps):
            protein_diff = np.zeros(len(self.proteins))
            for inter in self.intermediates:
                protein_diff += inter.intermediate_conc_update(self.proteins)
            for protein, p_diff in zip(self.proteins, protein_diff):
                protein.protein_conc_update(protein_diff)


def repair_pathway_factory(proteins, inters, init_damage=0.14):
    protein_objs = []
    intermediate_objs = []
    for protein in proteins:
        protein_objs.append(
            Protein(
                init_conc=protein["init_conc"],
                on_dict=protein["on"],
                off_dict=protein["off"],
                name=protein["name"]
            )
        )

    for intermediate in inters:
        intermediate_objs.append(
            Intermediates(
                enzyme_dict=intermediate["enzyme_dict"],
                name=intermediate["name"],
                mask=intermediate["mask"]
            )
        )
        if intermediate["name"] == "damaged":
            intermediate_objs[-1].conc[0] = init_damage

    for intermediate in intermediate_objs:
        for state, enzyme in intermediate.enzyme_dict.items():
            intermediate.trans_inter[state] = np.asarray(intermediate_objs)[enzyme["inter"]]

    return RepairPathway(protein_objs, intermediate_objs)


def main():
    proteins = [XPC, TFIIH, XPG, XPA, XPF, RPA, PCNA]
    intermediates = [DAMAGED, PARTIALLY, FULLY, INCISED, RESYNTH, RECHROM]
    repair_pathway = repair_pathway_factory(proteins, intermediates)
    repair_pathway.run()


if __name__ == '__main__':
    main()

