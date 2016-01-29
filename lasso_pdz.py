""" 
This script contains code for performing lasso minimization 
to integrate PDZ Domain sequences in our model
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

	def __init__(self):
		self.name = name 
		self.thetas = None 
		self.thresholds = None 
		self.sequence = None ## Not yet defined 
		self.seqlen = None ## Not yet defined

	def create_domains(self):
		temp_df = pd.read_excel(DATA+'\\theta_data.xlsx')
		self.df = temp_df.T
		self.thetas = self.df[self.name][:100]
		self.thresholds = np.asarray(self.df[self.name][100:])
		self.thetas = np.asarray(self.thetas)
		self.thetas = self.thetas.reshape(5,20)


