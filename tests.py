### Testing whether code functions well or not

from all_data import *
from run_mc import *

bn = 7.0
PDZ_Data.divide_peps()
calc_energy_ground()
test_peptide = PDZ_Data.peptide_dist[bn][0]
print test_peptide.name
print test_peptide.sequence_bis
print test_peptide.energy_ground

for peptide in PDZ_Data.peptide_dist[bn]:
    print "{} {}".format(peptide.name, peptide.sequence_bis)

run_mc(100, test_peptide, beta = 1.01, nb_cycles = 5, plot = True)

print compute_entropy_sequence(test_peptide)
