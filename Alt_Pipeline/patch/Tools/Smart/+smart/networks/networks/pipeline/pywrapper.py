# -------------------------------------------------------
# Python Wrapper for Importing MATLAB Data
# Created by: Jonathan Zia
# Date Created: Monday, March 25, 2019
# Inan Research Lab @ Gatech
# -------------------------------------------------------

import scipy.io as sio
import numpy as np
import h5py

# -------------------------------------------------------
# Loading MATLAB Data from Cell Array
# -------------------------------------------------------
def import_cell(filePath, fieldName, trainCells, testCells, sequential=False, HDF5=False, returnValid=True, trainIdx=[]):

	"""
	Importing cell array from MATLAB into Python

	Arguments:
	- filePath 		String 	Path to MATLAB file (including filename)
	- fieldName 	String 	Name of cell array in file
	- trainCells 	[Int] 	List of cells for training set
	- testCells 	[Int] 	List of cells for testing set
	- sequential 	Bool 	Is the training data sequential?
	- HDF5			Bool 	Is the data in HDF5 format?
	- trainIdx 		[Int] 	List of training indices to return
	"""

	# Load cell array into placeholder
	if not HDF5:
		cell = sio.loadmat(filePath)
	else:
		cell = h5py.File(filePath)

	# Return the field name only
	data = cell[fieldName]

	# Condense cell array into list of numpy arrays
	if data.size > 1:
		data = np.squeeze(data)
	
	# If there is only one training/testing index, return the single array
	if len(trainCells) == 1:
		training, testing = data, data
		return training, testing

	output = []
	for i in range(data.shape[0]):
		if HDF5:
			temp = cell[data[i]].value
			temp = np.moveaxis(temp, [0, 1, 2, 3], [3, 2, 1, 0])
			output.append(temp)
		else:
			output.append(data[i])
	
	# Separate the list into training and testing sets
	# Initialize placeholder for training data
	if len(trainCells) > 1:
		training = output[trainCells[0]]
	else:
		training = []
	# Initialize placeholder for testing data
	if len(testCells) > 1:
		testing = output[testCells[0]]
	else:
		testing = []

	# Populate placeholder for training data

	# Set axis based on type of data
	if sequential:
		axis = 0
	else:
		axis = 1

	# Populate placeholders
	if len(trainCells) > 1:
		for cell in trainCells[1:len(trainCells)]:
			training = np.concatenate((training, output[cell]), axis=axis)
	# Populate placeholder for testing data
	if len(testCells) > 1:
		for cell in testCells[1:len(testCells)]:
			testing = np.concatenate((testing, output[cell]), axis=axis)


	# Remove training data if indicated
	if len(trainIdx) > 0:
		training = training[trainIdx, :, :]


	# Return the output arrays
	return training, testing