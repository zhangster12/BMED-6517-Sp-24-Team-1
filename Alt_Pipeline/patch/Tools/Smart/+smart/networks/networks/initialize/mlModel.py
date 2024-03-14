# -------------------------------------------------------
# mlModel Object for Neural Networks
# Created By: Jonathan Zia
# Date Created: Wednesday, March 27, 2019
# Inan Research Lab @ Gatech
# -------------------------------------------------------

from keras.models import Model, load_model
from keras.layers import Input, Dense, Conv1D, MaxPooling1D, UpSampling1D, Flatten, Reshape, Lambda
from keras import backend as K
from keras.losses import mse
import os

class mlModel():

	"""
	Class containing general methods, member variables, and functions
	for neural network models
	"""

	def __init__(self, layers, units, activation, numIn, numOut, root, load=False, save_path='', load_path='', encoder=False):

		"""
		Specify model architecture and initialize Keras model

		ARGUMENTS:
		- layers 		Int 	Number of hidden layers
		- units 		[Int]	Number of units per hidden layer
		- activation 	[Str]	Activation function for each hidden layer and output layer
		- numIn 		Int 	Number of input features
		- numOut 		Int 	Number of output features
		- root 			Str 	Root directory
		- load 			Bool 	Load model from .h5 file?
		- save_path 	Str 	Save path for results/output from root
		- load_path 	Str 	Load path for .h5 file from root (including filename and extension)
		- encoder 		Int 	Layer number for encoder output layer (if indicated)
		"""

		# -------------------------------------------------------
		# Member Variables
		# -------------------------------------------------------

		# Architecture
		self.layers = layers
		self.units = units
		self.activation = activation
		self.inputs = numIn
		self.outputs = numOut
		self.encoder = encoder

		# Save/Load preferences
		self.root = root
		self.load = load
		self.save_path = os.path.join(self.root, save_path)
		self.load_path = os.path.join(self.root, load_path)

		# Filepaths
		self.root_dir = os.path.join(self.save_path, 'results')				# Save directory
		self.save_path_full = os.path.join(self.root_dir, 'full.h5')		# Saving full model
		if self.encoder != False:
			self.save_path_enc = os.path.join(self.root_dir, 'encoder.h5')	# Saving encoder only
		self.load_path_full = os.path.join(self.root, 'full.h5')			# Load path for full model
		if self.encoder != False:
			self.load_path_enc = os.path.join(self.root, 'encoder.h5')		# Load path for encoder model

		# Special flags
		self.expand = False 	# Indicating whether to expand input dimensions

	def construct(self):

		"""
		Construct Keras model after initialization
		"""

		# -------------------------------------------------------
		# Initialize Network
		# -------------------------------------------------------

		# Initialize list of layer objects
		layer = []

		# Define input layer
		layer.append(Input(shape=(self.inputs,)))

		# Define hidden layers
		for i in range(0, self.layers):
			layer.append(Dense(self.units[i], activation=self.activation[i])(layer[i]))

		# Define output layer
		output = Dense(self.outputs, activation=self.activation[-1])(layer[-1])

		# -------------------------------------------------------
		# Define the Model
		# -------------------------------------------------------

		# Load the full model, if indicated
		if self.load:
			network = load_model(self.load_path_full)
			# Load the encoder model, if indicated
			if self.encoder != False:
				encoder = load_model(self.load_path_enc)
		else:
			# Initialize the full model
			network = Model(inputs=layer[0], outputs=output)
			# Initialize the encoder, if indicated
			if self.encoder != False:
				encoder = Model(inputs=layer[0], outputs=layer[self.encoder])

		# -------------------------------------------------------
		# Return the Model
		# -------------------------------------------------------

		if self.encoder != False:
			return network, encoder
		else:
			return network

	def construct_conv(self, kernel=3, pool=2):

		"""
		Construct a convolutional autoencoder after initialization

		ARGUMENTS:
		- kernel 	Int 	Kernel length
		"""

		# Indicate expanded input dimensions
		self.expand = True

		# -------------------------------------------------------
		# Initialize Network
		# -------------------------------------------------------

		# Initialize list of layer objects
		layer = []

		# Define input layer
		layer.append(Input(shape=(self.inputs, 1)))
		inputs = layer[-1]
		
		# Define encoder layers
		# Hidden layers are composed of a convolutional layer and max pooling layer
		for i in range(0, self.layers):
			layer.append(Conv1D(self.units[i], kernel, activation='relu', padding='same')(layer[-1]))
			layer.append(MaxPooling1D(pool_size=pool)(layer[-1]))

		# Define latent layer
		layer.append(Dense(self.units[-1], activation=self.activation[0])(Flatten()(layer[-1])))
		latent = layer[-1]

		# Re-expand latent layer
		layer.append(Dense(int(layer[-2].shape[1]*layer[-2].shape[2]), activation=self.activation[0])(latent))
		
		# Expand dimensions of latent layer
		layer[-1] = Reshape((int(int(layer[-3].shape[1]*layer[-3].shape[2])/pool), pool))(layer[-1])

		# Define decoder layers
		for i in range(self.layers-1, -1, -1):
			layer.append(Conv1D(self.units[i], kernel, activation='relu', padding='same')(layer[-1]))
			layer.append(UpSampling1D(size=pool)(layer[-1]))

		# Define output layer
		layer.append(Conv1D(1, kernel, activation='relu', padding='same')(layer[-1]))
		output = layer[-1]

		# -------------------------------------------------------
		# Define the Model
		# -------------------------------------------------------

		# Load the full model, if indicated
		if self.load:
			network = load_model(self.load_path_full)
			# Load the encoder model, if indicated
			if self.encoder != False:
				encoder = load_model(self.load_path_enc)
		else:
			# Initialize the full model
			network = Model(inputs=inputs, outputs=output)
			# Initialize the encoder, if indicated
			if self.encoder != False:
				encoder = Model(inputs=inputs, outputs=latent)

		# -------------------------------------------------------
		# Return the Model
		# -------------------------------------------------------

		return network, encoder