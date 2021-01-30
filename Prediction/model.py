import os
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2' 
import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import layers
import numpy as np
import matplotlib.pyplot as plt
import math

data_input,data_output = [],[]

n_prediction = 20
shift = 50

def predictRandom(model,size,test_input,test_output):
	
	sample_in,sample_out = [],[]
	index = np.random.random(size)
	
	for i in range(size):
		index[i] = int(index[i]*100000)%len(test_input)
		sample_in.append(test_input[int(index[i])])
		sample_out.append(test_output[int(index[i])])

	sample_in = np.array(sample_in).reshape(size,len(test_input[0]))
	sample_out = np.array(sample_out).reshape(size,len(test_output[0]))

	prediction = model.predict(sample_in)
	for i in range(size):
		printLine_test(int(index[i]),test_input,test_output)
		print("Prediction : ",prediction[i])

	print(sample_out)
	print(prediction)

def printLine_test(i,test_input,test_output):
	print("Matrix n : ",test_input[i][0])
	print("Matrix nnz : ",test_input[i][1])
	print("Matrix sparsity : ",test_input[i][2])
	print("Matrix trace real : ",test_input[i][3])
	print("Matrix trace imag : ",test_input[i][4])
	print("Matrix amplitude mean : ",test_input[i][5])
	print("Matrix amplitude variance : ",test_input[i][6])
	print("Matrix min mod : ",test_input[i][7])
	print("Matrix max mod : ",test_input[i][8])
	print("Matrix diag min mod : ",test_input[i][9])
	print("Matrix diag max mod : ",test_input[i][10])
	print("Matrix diag mean : ",test_input[i][11])
	print("Matrix diag variance : ",test_input[i][12])
	print("Matrix diag gravity : ",test_input[i][13])
	print("Matrix left bandwidth : ",test_input[i][14])
	print("Matrix right bandwidth : ",test_input[i][15])

	if(test_input[i][16] == float(0)):
		print("Matrix solver : GMRES")
	else:
		print("Matrix solver : BICG")

	print("RHS n : ",test_input[i][17])
	print("RHS mean : ",test_input[i][18])
	print("RHS var : ",test_input[i][19])
	print("RHS gravity : ",test_input[i][20])
	print("RHS min mod : ",test_input[i][21])
	print("RHS max mod : ",test_input[i][22])
	#print("Converged : ",data_output[i][0])
	print("Iterations : ",test_output[i][0])
	#print("Elapsed Time : ",data_output[i][2])


def printLine(i):
	print("Matrix n : ",data_input[i][0])
	print("Matrix nnz : ",data_input[i][1])
	print("Matrix sparsity : ",data_input[i][2])
	print("Matrix trace real : ",data_input[i][3])
	print("Matrix trace imag : ",data_input[i][4])
	print("Matrix amplitude mean : ",data_input[i][5])
	print("Matrix amplitude variance : ",data_input[i][6])
	print("Matrix min mod : ",data_input[i][7])
	print("Matrix max mod : ",data_input[i][8])
	print("Matrix diag min mod : ",data_input[i][9])
	print("Matrix diag max mod : ",data_input[i][10])
	print("Matrix diag mean : ",data_input[i][11])
	print("Matrix diag variance : ",data_input[i][12])
	print("Matrix diag gravity : ",data_input[i][13])
	print("Matrix left bandwidth : ",data_input[i][14])
	print("Matrix right bandwidth : ",data_input[i][15])

	if(data_input[i][16] == float(0)):
		print("Matrix solver : GMRES")
	else:
		print("Matrix solver : BICG")

	print("RHS n : ",data_input[i][17])
	print("RHS mean : ",data_input[i][18])
	print("RHS var : ",data_input[i][19])
	print("RHS gravity : ",data_input[i][20])
	print("RHS min mod : ",data_input[i][21])
	print("RHS max mod : ",data_input[i][22])
	#print("Converged : ",data_output[i][0])
	print("Iterations : ",data_output[i][0])
	#print("Elapsed Time : ",data_output[i][2])


with open('/home/rustin/Projets/JSC_PI/v1.0/shuffled_sample.txt','r') as f:
    
    data_input_line,data_output_line = [],[]
    lines = f.readlines()
    count = 0
    for line in lines:
        count = count + 1
        atom = line.split(",")
        if(count > 1):
        	for i in range(len(atom)):
        		atom[i] = float(atom[i])

        	data_input_line  = atom[0:len(atom)-3]
        	data_output_line = atom[len(atom)-3	:len(atom)]

        	if(data_output_line[0] == 1):
        		data_input.append(data_input_line)
        		data_output.append(data_output_line[1:2])

        	input_size  = len(atom)-3
        	output_size = 1
        else:
        	header_input  = atom[0:len(atom)-3]
        	header_output = atom[len(atom)-3:len(atom)]



printLine(10)

train_test_ratio = 0.90
train_input  = np.array(data_input[0:int((count-1)*train_test_ratio)])
train_output = np.array(data_output[0:int((count-1)*train_test_ratio)])

print(len(train_input))

test_input  = np.array(data_input[int((count-1)*train_test_ratio):count-1])
test_output = np.array(data_output[int((count-1)*train_test_ratio):count-1])

model = keras.Sequential(name="Prediction_Model")

n1,n2,n3,n4 = 15,15,3,output_size
b_size,vs,epoch = 550,0.9,170

model.add(layers.Dense(n1, activation="relu", input_shape=(input_size,)))
model.add(layers.Dense(n2, activation="sigmoid"))
model.add(layers.Dense(n2, activation="relu"))
model.add(layers.Dense(n2, activation="relu"))
model.add(layers.Dense(n2, activation="relu"))
model.add(layers.Dense(n3, activation="relu"))
model.add(layers.Dense(n4, activation="relu"))


model.compile(
	optimizer="adam",
	loss=tf.keras.losses.MeanSquaredError())

metrics = model.fit(train_input,train_output,batch_size=b_size,validation_split=vs,epochs=epoch,shuffle=True)
results = model.evaluate(test_input,test_output,batch_size=64)

model.summary()

res = model.predict(test_input[0+shift:n_prediction+shift].reshape((n_prediction,23)))

print(b_size,vs,epoch)

print("Layers : ",n1,n2,n3,n4)

predictRandom(model,20,test_input,test_output)

val_loss = metrics.history['val_loss']