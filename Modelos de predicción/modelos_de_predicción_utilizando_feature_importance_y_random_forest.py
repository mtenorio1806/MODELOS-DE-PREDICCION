# -*- coding: utf-8 -*-
"""MODELOS DE PREDICCIÓN UTILIZANDO FEATURE IMPORTANCE Y RANDOM FOREST

Automatically generated by Colaboratory.

Original file is located at
    https://colab.research.google.com/drive/1MNgM-V84tSaxrlZtulooQmDrhZf1iAG-
"""

#MODELOS DE PREDICCIÓN UTILIZANDO FEATURE IMPORTANCE Y RANDOM FOREST

  #IMPORTAR LIBRERIAS NECESARIAS

import pandas as pd
import numpy as np
import sklearn.metrics as metrics
from sklearn.metrics import accuracy_score, confusion_matrix, plot_confusion_matrix, classification_report
from sklearn.model_selection import train_test_split
from sklearn.metrics import scorer
from sklearn.model_selection import StratifiedKFold, KFold
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.model_selection import GridSearchCV
from pandas import ExcelWriter
from sklearn.ensemble import RandomForestClassifier

                              #CARGAR DATOS

data = pd.read_csv('/content/drive/My Drive/MATRICES PR/CONJUNTOS DE DATOS_MODELO PR/ALELOS AL 95/TGCEUCAST_95_TEST_PR.csv')
#data.describe()
#print(data)

                       #PREPROCESAMIENTO DE LOS DATOS

  #Eliminar valores perdidos

data_na = data.dropna(axis=0)
#print(data_na)

  #Convertir las características a variables binarias 

lista = data_na.drop('TGC_EUCAST', axis=1)
datadum = pd.get_dummies(lista.astype(str))
#print(datadum)

  #Eliminar columnas duplicadas
  
for col in datadum.columns: 
  if col[-2:] == '_0':
    #print(col)
    del(datadum[col])

#print(datadum.columns)
#print(datadum.shape)

                                      #DIVISIÓN DE DATOS


  #Definir variables independientes (carateristicas) y dependientes (fenotipo)

X = datadum.values

y= data_na['TGC_EUCAST'].values 

  #Dividir el conjunto de datos en entrenamiento (80%) y prueba (20%) de forma estratificada segun el fenotipo

X_train, X_test, y_train, y_test = train_test_split(X, y, test_size= 0.2, stratify= y, random_state=42)
print("Valor a predecir: ", y_test)
print('X_train shape :', X_train.shape)


                                #ENTRENAR MODELO

  #Inicialiar variables

skf = StratifiedKFold(n_splits=10)
Ramfor = RandomForestClassifier(random_state=42)
Ramfor.fit(X_train, y_train)

 #Selección de características mediante Feature Importance

characteristic = datadum.columns
importances = Ramfor.feature_importances_
#print(importances)
# get importance
characteristics_importances = [(characteristic, round(importance, 2)) for characteristic, importance in zip(characteristic, importances)]
characteristics_importances = sorted(characteristics_importances, key = lambda x: x[1], reverse = True)
#[print('Variable: {:20} Importance: {}'.format(*pair)) for pair in characteristics_importances];
print(characteristics_importances)

my_df = pd.DataFrame(datadum.columns, [importances])
#print(my_df)
writer = ExcelWriter('/content/drive/My Drive/MATRICES PR/CONJUNTOS DE DATOS_MODELO PR/FEATURES/RF_FEATIMPO.xls')
my_df.to_excel(writer, 'Hoja de datos')
writer.save()

for dat in characteristics_importances:
  if dat[1] < 0.01:
    del(datadum[dat[0]])
print(dat)
print(X_train.shape)


                            #AJUSTE DE HIPERPARÁMETROS

  #Definir el rango de valores e hiperparámetros para la optimización

parameters = {'max_depth' : [2, 4, 6], 'min_samples_leaf': [1, 2, 4], 'min_samples_split' : [2, 5, 10], 'n_estimators': [50, 100, 200, 400], 'max_features': ['sqrt', 'auto']}
 
 #Entrenar el modelo para obtener los mejores hiperparámetros

model = GridSearchCV(Ramfor, parameters, cv=skf)
model.fit(X_train, y_train)
print('Parametros: ', model.best_params_)
print('Mejor puntaje', model.best_score_)
print(model.best_estimator_)

y_predicted = model.predict(X_test)

print('Accuracy test: ', accuracy_score(y_test, y_predicted))
print('Accuracy training: ', accuracy_score(y_train, model.predict(X_train)))

  #Matriz de confusion y reporte de clasificacion training

cm_train =confusion_matrix(y_train, model.predict(X_train))
print(cm_train)

print('Class report train: ', classification_report(y_train, model.predict(X_train)))

  #Matriz de confusion y reporte de clasificacion test

cm =confusion_matrix(y_test, y_predicted)
print(cm)

print(classification_report(y_test, y_predicted))

ax= plt.subplot()
sns.heatmap(cm, annot=True, ax = ax); #annot=True to annotate cells

# labels, title and ticks
ax.set_xlabel('Predicted');ax.set_ylabel('True'); 
ax.set_title('Confusion Matrix'); 
ax.xaxis.set_ticklabels(['Sensible','Resistente']); ax.yaxis.set_ticklabels(['Sensible', 'Resistente'])