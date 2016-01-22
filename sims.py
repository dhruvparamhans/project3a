### Simulations done over all peptides

from run_mc import *

nb_cycles = 1
beta = 1.0
nb_runs = 1000

test_peptide = PDZ_Data.peptides[11]
calc_energy_ground()

run_mc(nb_runs, test_peptide, beta, nb_cycles = nb_cycles, plot = True, verbose = False)

plot_freq_matrix(test_peptide)

print compute_entropy_sequence(test_peptide)

mi = compute_mutual_information(test_peptide)

plt.imshow(mi, interpolation='nearest', cmap = plt.cm.Blues, aspect = 'auto')
plt.colorbar()
plt.show()
