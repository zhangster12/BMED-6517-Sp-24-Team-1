# -------------------------------------------------------
# Testing mlModel Objects
# Created By: Jonathan Zia
# Date Created: Wednesday, March 27, 2019
# Inan Research Lab @ Gatech
# -------------------------------------------------------

import networks.pipeline.pywrapper as py
from keras.models import load_model
import numpy as np
import os

def test(model, dataset, testCells, batch=32, encoder=False,
	fieldName='data', predName='predictions.txt', tarName='targets.txt'):


	"""
	Method for testing mlModel objects

	ARGUMENTS
	- model 		mlModel 	Trained mlModel object
	- dataset 		String 		Path to testing dataset including filename and extension
	- testCells 	[Int]		List of cells on which to test
	- batch 		Int 		Batch size for predictions
	- encoder 		Bool 		Generate subspace representation from encoder model?
	- fieldName 	String 		Name of data field in .mat file for processing
	- predName 		String 		Name of file to store predictions ( + extension )
	- tarName 		String 		Name of file to store targets ( + extension )
	"""

	# -------------------------------------------------------
	# Get Testing Data
	# -------------------------------------------------------

	# Get testing inputs
	data, _ = py.import_cell(dataset, fieldName, testCells, [])
	data = np.transpose(data)
	# Expand data dimensions if necessary
	if model.expand:
		data = np.expand_dims(data, axis=-1)

	# Get testing targets
	targets = data

	# -------------------------------------------------------
	# Load Trained Model
	# -------------------------------------------------------
	if encoder == False:
		network = load_model(model.save_path_full)
	else:
		network = load_model(model.save_path_enc)

	# -------------------------------------------------------
	# Get Predictions
	# -------------------------------------------------------
	predictions = network.predict(data, batch_size=batch)

	# -------------------------------------------------------
	# Write Predictions and Targets to Files
	# -------------------------------------------------------

	# Write predictions
	with open(os.path.join(model.root_dir, predName), 'w') as file_object:
		np.savetxt(file_object, np.squeeze(predictions))

	# Write targets
	with open(os.path.join(model.root_dir, tarName), 'w') as file_object:
		np.savetxt(file_object, np.squeeze(targets))