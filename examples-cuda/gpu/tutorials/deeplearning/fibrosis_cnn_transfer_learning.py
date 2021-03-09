#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan  4 11:26:48 2017

@author: vkola
"""
from keras.preprocessing.image import ImageDataGenerator
from keras.models import Model
from keras.layers import Dense
from keras.layers import GlobalAveragePooling2D
from keras.utils import np_utils
from keras.callbacks import EarlyStopping, LearningRateScheduler
from keras.optimizers import SGD
import scipy.io
import numpy as np
import random
from sklearn.cross_validation import train_test_split
from sklearn.metrics import confusion_matrix, roc_curve, auc
from keras.applications.vgg16 import VGG16

random.seed(18)

def get_mat_data():
    matData = scipy.io.loadmat('fibrosis_creatinine_data.mat')
    xT = matData.get('xTrain')
    yT = matData.get('yTrain')
    return xT, yT
       

def generate_results(y_test, y_score):
    fpr, tpr, _ = roc_curve(y_test, y_score)
    roc_auc = auc(fpr, tpr)
    print('AUC: %f' % roc_auc)

    
def scheduler(epoch):
    if epoch == 5:
        model.lr.set_value(.005)
    return model.lr.get_value()

    
batch_size = 64
nb_classes = 2 
nb_epoch = 100 
data_augmentation =True

xT, yT = get_mat_data()
X_train, X_test, y_train, y_test = train_test_split(xT, yT, test_size=0.3, random_state=18)
print(X_train.shape, X_test.shape, y_train.shape, y_test.shape)

Y_train = np_utils.to_categorical(y_train, nb_classes)
Y_test = np_utils.to_categorical(y_test, nb_classes)

vgg16_model = VGG16(weights="imagenet", include_top=True)
base_model = Model(input=vgg16_model.input, 
                   output=vgg16_model.get_layer("block5_pool").output)

x = base_model.output
x = GlobalAveragePooling2D()(x)

x = Dense(1024, activation='relu')(x)

predictions = Dense(nb_classes, activation='softmax')(x)

for layer in base_model.layers[0:14]:
    layer.trainable = False

model = Model(input=base_model.input, output=predictions)

sgd = SGD(lr=1e-6, decay=1e-6, momentum=0.9, nesterov=True)

model.compile(loss='binary_crossentropy',
              optimizer='sgd',
              metrics=['accuracy'])

X_train = X_train.astype('float32')
X_test = X_test.astype('float32')
X_train /= 255
X_test /= 255

if not data_augmentation:
    print('Not using data augmentation.')
    early_stopping = EarlyStopping(monitor='val_loss', patience=5, verbose=1, mode='auto')
    history = model.fit(X_train, Y_train,
              batch_size=batch_size,
              nb_epoch=nb_epoch,
              validation_data=(X_test, Y_test),
              callbacks = [early_stopping],
              shuffle=True)    
else:
    print('Using real-time data augmentation.')
    datagen = ImageDataGenerator(
        rotation_range=180,  
        width_shift_range=0.1,   
        height_shift_range=0.1)  
    datagen.fit(X_train)
    early_stopping = EarlyStopping(monitor='val_loss', patience=5, verbose=1, mode='auto')
    change_lr = LearningRateScheduler(scheduler)
    history = model.fit_generator(datagen.flow(X_train, Y_train,
                        batch_size=batch_size),
                        samples_per_epoch=1000*X_train.shape[0],
                        nb_epoch=nb_epoch,
                        callbacks = [early_stopping],
                        validation_data=(X_test, Y_test))

    
# Post processing results
y_scores = model.predict(X_test)
y_preds = np.round(y_scores)
cm = confusion_matrix(y_test, y_preds[:,1])
print(cm)    

generate_results(y_test, y_scores[:, 1])
