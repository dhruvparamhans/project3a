## Functions for doing the monte-carlo runs

import numpy as np
from all_data import *


## Create the relevant objects

PDZ_Data = Data()
PDZ_Data.load_data()

def eval_score(domain, sequence, pos = 0):
    score = 0.0
    for i in range(5):
        score += domain.thetas[i,sequence[i]]
    return score - domain.thresholds[pos]

def sigmoid(x, a=1):
    return 1.0/(1+np.exp(-1.0*a*x))

def log_modified(x):
    if x > 0:
        return np.log(1+np.exp(-x))
    else:
        return -x+ np.log(1+np.exp(x))

def calc_log_proba_mod(peptide, domain, sequence):
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

def pair_quantities():
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
        temp1 = pep.sims[i]['Results']
        for j in range(len(temp1)):
            if temp[j]['Status'] == 'Accepted':
                accepted.append(temp1[j]['Sequence'])
    for mut in accepted:
        for i in range(5):
            freq_matrix[i,mut[i]] += 1
    freq_matrix_normalized = freq_matrix.astype('float') / freq_matrix.sum(axis=1)[:, np.newaxis]
    return [freq_matrix, freq_matrix_normalized]

def compute_proba_pos(pep,pos):
    return 1-compute_freq_matrix(pep)[1][pos,convert2int(pep.sequence_bis)[pos]]
