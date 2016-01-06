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

	def __init___(self, name):
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

class Data:
	"""
	Constructs Peptide and Domain objects for each peptide and domain. Also includes other relevant data such as the names of all the domains, peptides and aminoacids
	"""

	def __init__(self):
		## PDZ Domains 
		temp_df = pd.read_excel('Data_PDZ/theta_data.xlsx')
		self.aminoacids = [acid.encode('utf-8') for acid in list(temp_df.columns[:20])]
		self.df = temp_df.T
		self.domains = [Domain(domain.encode('utf-8')) for domain in list(self.df.columns)]
		self.domain_names = [domain.name for domain in self.domains]
		### Peptide sequences
		self.pep_seqs = []
		self.pep_names = []
		with open('Data_PDZ/peptides.free') as f:
			for line in f:
				x = line.split()
				self.pep_seqs.append(x[1])
				self.pep_names.append(x[0])
		self.peptides = [Peptide(name) for name in self.pep_names]

	def create_domains(self):
		for domain in self.domains:
			domain.thetas = self.df[domain.name][:100]
			domain.thresholds = np.asarray(self.df[domain.name][100:])
			domain.thetas = np.asarray(domain.thetas)
			domain.thetas = domain.thetas.reshape(5,20)

	def create_peptide(self):
		for i in range(len(self.pep_seqs)):
			self.peptides[i].sequence = self.pep_seqs[i]
			self.peptides[i].sequence_bis = list(self.pep_seqs[i])[5:]

	def load_data(self):
		self.create_domains()
		self.create_peptide()


