""" 
Data model for performing Lasso Minimization to integrate
PDZ Domain sequences into the model. 

I could have extended the previous data class but then 
lets concentrate on getting things right without worrying
too much about optimized code
"""


ROOT = 'E:\\Ecole\\Year 3\\Projet 3A'
DATA = ROOT+'\Data_PDZ'
import pandas as pd
import numpy as np

class Domain: 
	"""
	Domain object. Similar to the one used in all_data.py 
	Will eventually contain sequences as well
	""" 

	def __init__(self, name, seq_length = 50):
		self.name = name 
		self.thetas = None 
		self.thresholds = None 
		self.sequence = '' ## Not yet defined 
		self.seqlen = 0 ## Not yet defined
		self.sequence_list = list()
		self.sequence_number = list()
		self.seq_length = seq_length

	def remove_unwanted(self):
		while '\n' in self.sequence_list:
			self.sequence_list.remove('\n')

	def compute_one_hot(self):
		sigma = np.zeros((self.seq_length,20))
		for i in range(self.seq_length):
			sigma[i,self.sequence_number[i]] +=1
		self.sigma = sigma

class Data:

	def __init__(self):
		temp_df = pd.read_excel(DATA+'\\theta_data_removed.xlsx')
		self.aminoacids = [acid.encode('utf-8') for acid in list(temp_df.columns[:20])]
		self.theta_df = temp_df.T
		self.domains = [Domain(domain.encode('utf-8')) for domain in list(self.theta_df.columns)]
		self.domain_names = [domain.name for domain in self.domains]

		seq_df = pd.read_excel(DATA+'\\seq_pdz.xlsx')
		self.seq_df = seq_df.T

		self.encoding = 'utf-8'

	def create_domains(self):
		for domain in self.domains:
			domain.thetas = self.theta_df[domain.name][:100]
			domain.thetas = np.asarray(domain.thetas)
			domain.thetas = domain.thetas.reshape(5,20)

		for i in range(len(self.domain_names)):
			self.domains[i].sequence = self.seq_df[i]['Sequence'].encode(self.encoding)
			self.domains[i].sequence_list = list(self.domains[i].sequence)
			self.domains[i].remove_unwanted()
			self.domains[i].seqlen = len(self.domains[i].sequence_list)
			self.domains[i].sequence_number = self.convert2int(self.domains[i].sequence_list)
			self.domains[i].compute_one_hot()

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
		return [self.aminoacids.index(pep) for pep in seq_pep]

