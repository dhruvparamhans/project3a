"""
This script functions and classes that I need to declare every time I run the simulations. It is more convenient to use the python notebook for demonstration purposes, and keep the code here
"""


ROOT = 'E:\\Ecole\\Year 3\\Projet 3A'
DATA = ROOT+'\Data_PDZ'
import pandas as pd
import numpy as np

class Domain:
	"""
	Holds stuff relevant to the domains
	"""

	def __init__(self, name):
		self.name = name
		self.thresholds = None
		self.thetas = None


class Peptide:
	"""
	Peptide stuff
	"""

	def __init__(self, name):
		self.name = name
		self.sequence = None
		self.sequence_bis = None ## Last five amino acids
		self.energy_ground = None ## Useful for previous simulations
		self.score = None ## Score calculated from Stiffler model
		self.y_manip_bind = 0.0
		self.y_manip_nbind = 0.0
		self.y_model_bind = 0.0
		self.y_model_nbind = 0.0
		self.posterior_matrix = None

class Data:
	"""
	Constructs Peptide and Domain objects for each peptide and domain. Also includes other relevant data such as the names of all the domains, peptides and aminoacids
	"""

	def __init__(self):
		## PDZ Domains
		temp_df = pd.read_excel(DATA+'\\theta_data.xlsx')
		self.aminoacids = [acid.encode('utf-8') for acid in list(temp_df.columns[:20])]
		self.df = temp_df.T
		self.domains = [Domain(domain.encode('utf-8')) for domain in list(self.df.columns)]
		self.domain_names = [domain.name for domain in self.domains]
		### Peptide sequences
		self.pep_seqs = []
		self.pep_names = []
		with open(DATA+'\\peptides.free') as f:
			for line in f:
				x = line.split()
				self.pep_seqs.append(x[1])
				self.pep_names.append(x[0])
		self.peptides = [Peptide(name) for name in self.pep_names]

		## Interaction: Which peptides bind to which domains
		self.fp_interaction_matrix = pd.read_excel(DATA+"\\fp_interaction_matrix.xlsx")
		for column in self.fp_interaction_matrix.columns:
			self.fp_interaction_matrix.loc[self.fp_interaction_matrix[column] == 0.0, column] = -1.0
		self.fp_interaction_matrix = self.fp_interaction_matrix.rename(columns = lambda x: str(x).replace(" ", ""))

		## Classification matrix
		self.class_matrix = np.zeros((2,2))
		self.class_matrix[0,0] = 0.85
		self.class_matrix[0,1] = 0.04
		self.class_matrix[1,0] = 0.15
		self.class_matrix[1,1] = 0.96

	def create_domains(self):
		for domain in self.domains:
			domain.thetas = self.df[domain.name][:100]
			domain.thresholds = np.asarray(self.df[domain.name][100:])
			domain.thetas = np.asarray(domain.thetas)
			domain.thetas = domain.thetas.reshape(5,20)

	def convert2seq(self,seq_int):
		"""
		Function to convert integer sequences into letter sequences
		Useful for visualization purposes
		"""
		return [self.aminoacids[i] for i in seq_int]

	def convert2int(self,seq_pep):
		"""
		Function to convert letter sequences into integer sequences
		Useful when doing Monte-carlo runs
		"""
		return [self.aminoacids/index(pep) for pep in seq_pep]

	def create_peptide(self):
		for i in range(len(self.pep_seqs)):
			self.peptides[i].sequence = self.pep_seqs[i]
			self.peptides[i].sequence_bis = list(self.pep_seqs[i])[5:]

	def calc_y_manip(self):
		"""
		Calculate the binding probabilities using the updated probability rule.
		The updated rule models the false positive and true positive rates as
		probabilities and integrates them into the model.
		"""
		bind = 0.0
		nbind = 0.0
		for peptide in self.peptides:
			peptide.y_manip_bind = 0.0
			peptide.y_manip_nbind = 0.0
			for i in range(len(self.domain_names)):
				alpha = self.fp_interaction_matrix[peptide.name][i]
				if alpha > 0:
					peptide.y_manip_bind +=1
				else:
					peptide.y_manip_nbind +=1
			bind += peptide.y_manip_bind
			nbind += peptide.y_manip_nbind

			peptide.y_manip_bind /=74
			peptide.y_manip_nbind /=74

			peptide.y_model_bind = self.class_matrix[1,1]*peptide.y_manip_bind + self.class_matrix[1,0]*peptide.y_manip_nbind
			peptide.y_model_nbind = self.class_matrix[0,0]*peptide.y_manip_nbind + self.class_matrix[0,1]*peptide.y_manip_bind

			peptide.posterior_matrix = np.zeros((2,2))

			peptide.posterior_matrix[0,0] = self.class_matrix[0,0]*peptide.y_manip_nbind / peptide.y_model_nbind
			peptide.posterior_matrix[1,0] = self.class_matrix[0,1]*peptide.y_manip_bind / peptide.y_model_nbind

			peptide.posterior_matrix[0,1] = self.class_matrix[1,0]*peptide.y_manip_nbind / peptide.y_model_bind
			peptide.posterior_matrix[1,1] = self.class_matrix[1,1]*peptide.y_manip_bind / peptide.y_model_bind

		self.y_manip_bind = bind/np.size(self.fp_interaction_matrix)
		self.y_manip_nbind = nbind/np.size(self.fp_interaction_matrix)

		self.y_model_bind = self.class_matrix[1,1]*self.y_manip_bind + self.class_matrix[1,0]*self.y_manip_nbind

		self.y_model_nbind = self.class_matrix[0,0]*self.y_manip_nbind + self.class_matrix[0,1]*self.y_manip_bind


		## Posterior probability matrix

		self.posterior_matrix = np.zeros((2,2))
		## P(y_manip|y_model) = P(y_model|y_manip)*P(y_manip) / P(y_model)
		## P(manip=-1|model=-1) = P(model=-1|manip=-1) * P(manip=-1) /P(model=-1)
		self.posterior_matrix[0,0] = self.class_matrix[0,0]*self.y_manip_nbind / self.y_model_nbind
		## P(manip=1|model=-1) = P(model=-1|manip=1) * P(manip=1) /P(model=-1)
		self.posterior_matrix[1,0] = self.class_matrix[0,1]*self.y_manip_bind / self.y_model_nbind
		## P(manip=-1|model=1) = P(model=1|manip=-1) * P(manip=-1) /P(model=1)
		self.posterior_matrix[0,1] = self.class_matrix[1,0]*self.y_manip_nbind / self.y_model_bind
		## P(manip=1|model=1) = P(model=1|manip=1) * P(manip=1) /P(model=1)
		self.posterior_matrix[1,1] = self.class_matrix[1,1]*self.y_manip_bind / self.y_model_bind

	def divide_peps(self):
		"""
		Function which divides the peptides according to their binding number.
		The binding number for any peptide is the number of domains that it
		binds to according to fp_interaction_matrix.

		Once the simulations have been run, it will be useful to see the behavior
		of peptides as a function of their binding number 
		"""
		self.binds = [peptide.y_manip_bind*74 for peptide in self.peptides]
		from collections import OrderedDict
		self.peptide_dist = OrderedDict()
		for value in self.binds:
			self.peptide_dist[value] = []
		for peptide in self.peptides:
			x = peptide.y_manip_bind*74
			self.peptide_dist[x].append(peptide)
		return self.peptide_dist

	def load_data(self):
		"""
		Function which loads the relevant data to start performing analysis
		"""
		self.create_domains()
		self.create_peptide()
		self.calc_y_manip()
		self.divide_peps()
