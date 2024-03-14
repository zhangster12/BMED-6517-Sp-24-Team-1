# -------------------------------------------------------
# Training mlModel Objects
# Created By: Jonathan Zia
# Date Created: Wednesday, March 27, 2019
# Inan Research Lab @ Gatech
# -------------------------------------------------------

from keras.callbacks import ReduceLROnPlateau, EarlyStopping
from keras.callbacks import ModelCheckpoint, TensorBoard
from keras.utils.vis_utils import plot_model
from keras import optimizers, regularizers
import networks.pipeline.pywrapper as py
#import matplotlib.pyplot as plt
import numpy as np
import os

def train(model, network, dataset, epoch, lr, batch, decay, d_patience, s_patience, trainCells, validCells, graphics=False, encoder=False):

	"""
	Method for training mlModel objects

	ARGUMENTS
	- model 		mlModel 	Neural network object
	- network 		mlModel 	Constructed neural network
	- dataset 		Str 		Path to dataset including filename and extension
	- epoch 		Int 		Maximum number of epochs
	- lr 			Double 		Initial learning rate
	- batch 		Int 		Batch size
	- decay 		Double 		Learning rate decay (newLR = decay * oldLR)
	- d_patience 	Int 		Number of epochs before decaying learning rate
	- s_patience 	Int 		Number of epochs before training termination
	- trainCells 	[Int]		List of cells on which to train
	- validCells 	[Int]		List of cells on which to validate
	- graphics 		Bool 		Plot results?
	- encoder 		mlModel 	Constructed encoder
	"""

	# -------------------------------------------------------
	# Compile model and Deine Optimizer
	# -------------------------------------------------------

	# Defining Adam optimizer
	optimizer = optimizers.Adam(lr=lr)

	# Compile model
	network.compile(loss='mean_squared_error', optimizer=optimizer, metrics=['accuracy'])

	# Print a model summary
	print(network.summary())

	# -------------------------------------------------------
	# Specify the Input and Validation Data
	# -------------------------------------------------------

	# Get training and validation data
	x_train, x_valid = py.import_cell(dataset, 'az_rec', trainCells, validCells)
	y_train, y_valid = py.import_cell(dataset, 'az_rec', trainCells, validCells)
	x_train, x_valid, y_train, y_valid = np.transpose(x_train), np.transpose(x_valid), np.transpose(y_train), np.transpose(y_valid)

	# Expand data dimensions if necessary
	if model.expand:
		x_train, x_valid = np.expand_dims(x_train, axis=-1), np.expand_dims(x_valid, axis=-1)
		y_train, y_valid = np.expand_dims(y_train, axis=-1), np.expand_dims(y_valid, axis=-1)

	# -------------------------------------------------------
	# Run the Model
	# -------------------------------------------------------

	# Obtain callbacks during training. Callbacks include:
	# ReduceLROnPlateau 	Learning rate decay
	# EarlyStopping 		Early stopping
	# ModelCheckpoint 		Model checkpoints (optimal full model)
	# TensorBoard 			Writing TensorBoard output file

	# Declare callbacks
	lr_decay = ReduceLROnPlateau(factor=decay, patience=d_patience)
	early_stop = EarlyStopping(patience=s_patience, restore_best_weights=True)
	checkpoint = ModelCheckpoint(model.save_path_full, save_best_only=True)
	tensorboard = TensorBoard(log_dir=os.path.join(model.root_dir, 'logs'))
	callbacks = [lr_decay, early_stop, checkpoint, tensorboard]

	# Run model
	history = network.fit(x=x_train, y=y_train, batch_size=batch, 
		epochs=epoch, validation_data=(x_valid, y_valid), callbacks=callbacks)

	# Save the model (optimal, since weights are restored)
	network.save(model.save_path_full)
	if model.encoder != False:
		encoder.save(model.save_path_enc)

	# -------------------------------------------------------
	# Record Training Statistics
	# -------------------------------------------------------

	# Write training and validation losses to a file
	with open(os.path.join(model.root_dir, 'train_loss.txt'), 'w') as file_object:
		np.savetxt(file_object, history.history['loss'])
	with open(os.path.join(model.root_dir, 'valid_loss.txt'), 'w') as file_object:
		np.savetxt(file_object, history.history['val_loss'])

	# -------------------------------------------------------
	# Visualize Training
	# -------------------------------------------------------

	plot_model(network, to_file=os.path.join(model.root_dir,"graph.png"), 
		show_shapes=True, show_layer_names=True)

	# Plot the loss
	if graphics:
		plt.plot(history.history['loss'])
		plt.plot(history.history['val_loss'])
		plt.title('Model Loss')
		plt.ylabel('Loss')
		plt.xlabel('Epoch')
		plt.legend(['train','test'], loc='upper right')
		# Save plot loss to file and display the plot
		plt.savefig(os.path.join(model.root_dir, 'loss.png'))
		plt.show()