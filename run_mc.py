## Functions for doing the monte-carlo runs

import numpy as np
from all_data import *
import matplotlib.pyplot as plt
## Create the relevant objects

print "Creating Data Object PDZ_Data"
PDZ_Data = Data()
PDZ_Data.load_data()

print "PDZ_Data ready!"

def eval_score(domain, sequence, pos = 0):
    """
    Function which evaluates the score using Stiffler model
    """
    score = 0.0
    for i in range(5):
        score += domain.thetas[i,sequence[i]]
    return score - domain.thresholds[pos]

### Some utility functions for doing math
def sigmoid(x, a=1):
    return 1.0/(1+np.exp(-1.0*a*x))

def log_modified(x):
    if x > 0:
        return np.log(1+np.exp(-x))
    else:
        return -x+ np.log(1+np.exp(x))

def calc_log_proba_mod(peptide, domain, sequence):
    """
    Function which computes the log of the updated probability.
    For numerical stability, the sum of the logs is computed as
    log(exp(logA)+exp(logB)) which is just log(A+B)
    """
    ix = PDZ_Data.domain_names.index(domain.name)
    alpha = PDZ_Data.fp_interaction_matrix[peptide.name][ix]

    score = eval_score(domain, sequence,0)
    z_1 = log_modified(score)
    z_2 = log_modified(-1.0*score)
    if alpha > 0:
        a = peptide.posterior_matrix[1,1]
        x = np.log(a) -z_1
        b = peptide.posterior_matrix[1,0]
        y = np.log(b) - z_2
        result = np.logaddexp(x,y)
    else:
        a = peptide.posterior_matrix[0,1]
        x = np.log(a) - z_1
        b = peptide.posterior_matrix[0,0]
        y = np.log(b) - z_2
        result = np.logaddexp(x,y)
    return result*-1.0

def eval_log_energy(peptide,sequence):
    en = 0.0
    for domain in PDZ_Data.domains:
        en += calc_log_proba_mod(peptide, domain, sequence)
    return en

def calc_energy_ground():
    for peptide in PDZ_Data.peptides:
        peptide.energy_ground = eval_log_energy(peptide, PDZ_Data.convert2int(peptide.sequence_bis))


def pair_quantities():
    """
    Function which computes relevant quantities for each peptide-domain pair
    """
    for peptide in PDZ_Data.peptides:
        peptide.domain_data = {}
        for i in range(len(PDZ_Data.domain_names)):
            quant_list = []
            alpha = PDZ_Data.fp_interaction_matrix[peptide.name][i]
            if alpha > 0:
                quant_list.append(1.0)
            else:
                quant_list.append(alpha)
            energy = eval_score(domain, convert2int(peptide.sequence_bis),2)
            quant_list.append(energy) ## Score Calculated from Stiffler Model
            ## P(y_manip=1|seq) = P(y_manip=1|y_model = -1)*P(mod=-1|seq) + P(y_manip=1|y_model=1)*P(mod=1|seq)
            manip_1 = peptide.posterior_matrix[1,0]*sigmoid(energy,-1) + peptide.posterior_matrix[1,1]*sigmoid(energy)
            ## P(y_manip=-1|seq) = P(y_manip=-1|y_model = -1)*P(mod=-1|seq) + P(y_manip=-1|y_model=1)*P(mod=1|seq)
            manip_0 = peptide.posterior_matrix[0,0]*sigmoid(energy,-1) + peptide.posterior_matrix[0,1]*sigmoid(energy)
            quant_list.append(manip_1)
            quant_list.append(manip_0)
        peptide.domain_data[domain.name] = quant_list

def run_mc(nb_runs, peptide,beta = 1.0, nb_cycles=10, plot=False, verbose=True, print_reject = False):
    sims = []
    print "Name of Peptide {}".format(peptide.name)
    print "Base Energy {}".format(peptide.energy_ground)
    print "Base Sequence {}".format(peptide.sequence_bis)
    base_seq = PDZ_Data.convert2int(peptide.sequence_bis)
    peptide.mutations = []
    for j in range (nb_cycles):
        print "\n Cycle number : {}\n".format(j+1)
        sim_results = []
        mutated_sequences = []
        mutated_energies = []
        for_plot = []
        sequences_accepted = []
        mut_seq = base_seq
        mut_energy = peptide.energy_ground
        for i in range(nb_runs):
            y = np.random.randint(5)
            z = np.random.randint(19)
            ## Remove if the amino acid change is the same as before
            if z >= mut_seq[y]:
                z = z+1
            temp_seq = mut_seq[:]
            #print "Last sequence seen {}".format(convert2seq(mut_seq))
            temp_seq[y] = z
            #print "Sequence after mutation {}\n".format(convert2seq(temp_seq))
            temp_energy = eval_log_energy(peptide, temp_seq)
            ratio = np.exp(-1.0*beta*(temp_energy-mut_energy))
            prob_trans = min(1, ratio)
            x = np.random.uniform()
            if x < prob_trans:
                mut_energy = temp_energy
                mut_seq = temp_seq
                if verbose:
                    print "Run number: {}\n".format(i)
                    print "Uniform {} Ratio {} Prob_Trans {} ".format(x,ratio,prob_trans)
                    print "Accepted {} {} {} {} \n".format(temp_seq, temp_energy, PDZ_Data.convert2seq(temp_seq), PDZ_Data.convert2seq(mut_seq))
                sim_results.append({'Sequence': temp_seq, 'Energy': temp_energy, 'Status': 'Accepted'})
                for_plot.append(temp_energy)
            else:
                if verbose:
                    if print_reject:
                        print "Run number: {}\n".format(i)
                        print "Uniform {} Ratio {} Prob_Trans {} ".format(x,ratio,prob_trans)
                        print "Rejected {} {} {} {}\n".format(temp_seq, temp_energy, PDZ_Data.convert2seq(temp_seq), PDZ_Data.convert2seq(mut_seq))
                sim_results.append({'Sequence': temp_seq, 'Energy': temp_energy, 'Status': 'Rejected'})

                ##print "Rejected {} {}".format(temp_seq, temp_energy)
            mutated_sequences.append(temp_seq)
            mutated_energies.append(temp_energy)
        peptide.mutations.append(sequences_accepted)
        print "Lowest Energy {} Sequence {}\n".format(np.min(mutated_energies), PDZ_Data.convert2seq(mutated_sequences[np.argmin(mutated_energies)]))
        if plot == True:
            plt.figure(j)
            plt.axhline(y = peptide.energy_ground, hold = None, c = 'r', linewidth=0.5)
            if len(for_plot) == 1:
                plt.axhline(y = for_plot[0], hold=None, c = 'b', linewidth = 1.5)
            plt.plot(for_plot)
            plt.show()
        sims.append({'Results' : sim_results, 'Mutated sequences': mutated_sequences, 'Mutated Energies': mutated_energies})
    peptide.sims = sims
    print " Completed run for {}\n".format(peptide.name)

### Functions to be used once the simulations have been run
def compute_unique_mut(peptide):
    unique = []
    for run in peptide.mutations:
        for item in run:
            if item not in unique:
                unique.append(item)
    return unique

def compute_freq_matrix(pep):
    freq_matrix = np.zeros((5,20))
    accepted = []
    for i in range(len(pep.sims)):
        temp = pep.sims[i]['Results']
        for j in range(len(temp)):
            if temp[j]['Status'] == 'Accepted':
                accepted.append(temp[j]['Sequence'])
    for mut in accepted:
        for i in range(5):
            freq_matrix[i,mut[i]] += 1
    freq_matrix_normalized = freq_matrix.astype('float') / freq_matrix.sum(axis=1)[:, np.newaxis]
    return [freq_matrix, freq_matrix_normalized]

def compute_proba_pos(pep,pos):
    return 1-compute_freq_matrix(pep)[1][pos,convert2int(pep.sequence_bis)[pos]]

def plot_freq_matrix(pep, normalized = True):
    if normalized:
        fm_normalized = compute_freq_matrix(pep)[1]
    else:
        fm_normalized = compute_freq_matrix(pep)[0]
    print "Sequence of Peptide {}".format(pep.sequence_bis)
    plt.figure(num=None, figsize=(8, 6), dpi=80, facecolor='w', edgecolor='k')
    plt.title("Frequency matrix for Peptide {}".format(pep.name))
    plt.xticks(range(len(PDZ_Data.aminoacids)), PDZ_Data.aminoacids, fontsize=12)
    plt.imshow(fm_normalized,interpolation='nearest', cmap = plt.cm.Blues)
    plt.colorbar()
    plt.show()

def print_seq(pep):
    """
    Returns the sequence with the names of the amino acids
    """
    for acid in pep.sequence_bis:
        print "{} {} ".format(acid,acid_dict[acid])

def compute_entropy_sequence(peptide):
    test_matrix = compute_freq_matrix(peptide)[1]
    entropy_sequence = []
    x, y = test_matrix.shape
    for i in range(x):
        w = 0.0
        for j in range(y):
            if test_matrix[i][j] == 0:
                w +=0
            else:
                w += test_matrix[i][j] * np.log(test_matrix[i][j])
        w = -1.0*w
        entropy_sequence.append(w)
    peptide.entropy_sequence = entropy_sequence
    return entropy_sequence
