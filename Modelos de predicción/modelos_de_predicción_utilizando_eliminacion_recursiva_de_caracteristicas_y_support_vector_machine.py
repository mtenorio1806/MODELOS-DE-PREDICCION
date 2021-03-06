# -*- coding: utf-8 -*-
"""MODELOS DE PREDICCIÓN UTILIZANDO ELIMINACION RECURSIVA DE CARACTERISTICAS Y SUPPORT VECTOR MACHINE

Automatically generated by Colaboratory.

Original file is located at
    https://colab.research.google.com/drive/1fpghZXa-oXw9tDBPWsVrZbU5cXx0D0yt
"""

#MODELOS DE PREDICCIÓN UTILIZANDO ELIMINACION RECURSIVA DE CARACTERISTICAS Y SUPPORT VECTOR MACHINE

                         #IMPORTAR LIBRERIAS NECESARIAS

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import accuracy_score, confusion_matrix, plot_confusion_matrix, classification_report
from sklearn.model_selection import train_test_split
from sklearn.model_selection import StratifiedKFold, KFold
from sklearn.model_selection import GridSearchCV
from pandas import ExcelWriter
from sklearn.feature_selection import RFECV
from sklearn.svm import LinearSVC
from sklearn.svm import SVC

                         #CARGAR CONJUNTO DE DATOS

data = pd.read_csv('/content/drive/MyDrive/MATRICES PR/CONJUNTOS DE DATOS_MODELO PR/ALELOS 97/MER_97_TEST_PR.csv')
#data.describe()
#print(data)

                       #PREPROCESAMIENTO DE LOS DATOS

  #Eliminar valores perdidos

data_na = data.dropna(axis=0)
print(data_na)


  #Convertir las características a variables binarias

lista = data_na.drop('MER', axis=1)
print ('lista: ', lista)
datadum = pd.get_dummies(lista.astype(str))
print('datadum: ', datadum)

  #Eliminar columnas duplicadas 

for col in datadum.columns: 
  if col[-2:] == '_0':
    print(col)
    del(datadum[col])

print('datadum columns',datadum.columns)
print('datadum shape:', datadum.shape)

                        #DIVISIÓN DE DATOS

  #Definir variables independientes (carateristicas) y dependientes (fenotipo)

X = datadum.values

y= data_na['MER'].values 

  #Dividir el conjunto de datos en entrenamiento (80%) y prueba (20%) de forma estratificada segun el fenotipo

X_train, X_test, y_train, y_test = train_test_split(X, y, test_size= 0.2, stratify= y, random_state=42)
print("Conjunto de prueba: ", y_test)





                                #ENTRENAR MODELO


#Inicialiar variables

skf = StratifiedKFold(n_splits=10)
svm = SVC(kernel='linear', probability=True, C=0.01)  #Utilizar valor del parametro C obtenido de la curva de validación

svm.fit(X_train, y_train)
y_pred = svm.predict(X_test)
print(accuracy_score(y_test, y_pred))
print(accuracy_score(y_train, svm.predict(X_train)))

 #Selección de características mediante eliminación recursiva de características.

rfecv = RFECV(estimator=svm, cv=skf, scoring='accuracy')
rfecv.fit(X_train, y_train)
y_rfecv = rfecv.predict(X_test)

print('Optimal features :', rfecv.n_features_)
print('Best features: ', datadum.columns[rfecv.support_])

feature_impo = pd.DataFrame(datadum.columns[rfecv.support_])
writer = ExcelWriter('/content/drive/My Drive/MATRICES PR/CONJUNTOS DE DATOS_MODELO PR/FEATURES/FEATURES_SVM_REFCV.xls')
feature_impo.to_excel(writer, 'Hoja de datos')
writer.save()

  #Mantener las caracteristicas seleccionadas como importantes

X_train_rfecv = rfecv.transform(X_train)
X_test_rfecv = rfecv.transform(X_test)
print('X transformado', X_train_rfecv.shape)

                            #AJUSTE DE HIPERPARÁMETROS

  #Definir el rango de valores e hiperparámetros para la optimización

parameters = {'C':[0.0001, 0.001, 0.01, 0.1, 1, 10, 100], 'kernel': ['linear', 'poly', 'rbf']}

  #Entrenar el modelo para obtener los mejores hiperparámetros

model = GridSearchCV(svm, parameters, cv=skf)
model.fit(X_train_rfecv, y_train)
print('Parametros: ', model.best_params_)
print('Mejor puntaje', model.best_score_)
print(model.best_estimator_)

y_predicted = model.predict(X_test_rfecv)

print('Accuracy test', accuracy_score(y_test, y_predicted))
print('Accuracy training', accuracy_score(y_train, model.predict(X_train_rfecv)))

# Matriz de confusion y reporte de clasificacion del entrenamiento

cm_train =confusion_matrix(y_train, model.predict(X_train_rfecv))
print(cm_train)

print('Class report train: ', classification_report(y_train, model.predict(X_train_rfecv)))

# Matriz de confusion y reporte de clasificacion de prueba

cm =confusion_matrix(y_test, y_predicted)
print(cm)

print('Class report test: ', classification_report(y_test, y_predicted))

ax= plt.subplot()
sns.heatmap(cm, annot=True, ax = ax); #annot=True to annotate cells

# labels, title and ticks
ax.set_xlabel('Predicted');ax.set_ylabel('True'); 
ax.set_title('Confusion Matrix'); 
ax.xaxis.set_ticklabels(['Sensible','Resistente']); ax.yaxis.set_ticklabels(['Sensible', 'Resistente'])
