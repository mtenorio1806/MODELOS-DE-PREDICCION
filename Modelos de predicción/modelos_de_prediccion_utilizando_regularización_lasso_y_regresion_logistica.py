# -*- coding: utf-8 -*-
"""MODELOS DE PREDICCION UTILIZANDO REGULARIZACIÓN LASSO Y REGRESION LOGISTICA

Automatically generated by Colaboratory.

Original file is located at
    https://colab.research.google.com/drive/1QBYheoXuXFHUPagbQLE5j5FBbQ9NblwK
"""

#MODELOS DE PREDICCION UTILIZANDO REGULARIZACIÓN LASSO Y REGRESION LOGISTICA

  #IMPORTAR LIBRERIAS NECESARIAS

import pandas as pd
import numpy as np
import sklearn.metrics as metrics
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score, confusion_matrix, plot_confusion_matrix, classification_report
from sklearn.model_selection import train_test_split
from sklearn.model_selection import StratifiedKFold, KFold
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.model_selection import GridSearchCV
from sklearn.feature_selection import SelectFromModel
from pandas import ExcelWriter


                              #CARGAR DATOS

data = pd.read_csv('/content/drive/MyDrive/MATRICES PR/CONJUNTOS DE DATOS_MODELO PR/ALELOS AL 100/TGCEUCAST_100_PR_TEST.csv')
#data.describe()
#print(data)

                       #PREPROCESAMIENTO DE LOS DATOS

#Eliminar missing values

data_na = data.dropna(axis=0)
print(data_na)

#Convertir las características a variables binarias

lista = data_na.drop('TGC_EUCAST', axis=1)
datadum = pd.get_dummies(lista.astype(str))
print(datadum)

#Eliminar columnas duplicadas

for col in datadum.columns: 
  if col[-2:] == '_0':
    print(col)
    del(datadum[col])

print(datadum.columns)
print(datadum.shape)

                       #DIVISIÓN DE DATOS
#Definir variables independientes (carateristicas) y dependientes (fenotipo)

X = datadum.values

y= data_na['TGC_EUCAST'].values 

#Dividir el conjunto de datos en entrenamiento (80%) y prueba (20%) de forma estratificada segun el fenotipo

X_train, X_test, y_train, y_test = train_test_split(X, y, test_size= 0.2, stratify= y, random_state=42)
print("Conjunto de prueba: ", y_test)


                                #ENTRENAR MODELO

#Selección de características mediante regularización L1

skf = StratifiedKFold(n_splits=10)
lr = LogisticRegression(penalty='l1', solver='saga', C=0.1) #Utilizar valor del parametro C obtenido de la curva de validación
lr.fit(X_train, y_train)

sel_lasso = SelectFromModel(lr)
sel_lasso.fit(X_train, y_train)
criterio = sel_lasso.get_support()
#print(criterio_feat)
my_df = pd.DataFrame(datadum.columns, [criterio])
#print(my_df)

writer = ExcelWriter('/content/drive/My Drive/MATRICES PR/CONJUNTOS DE DATOS_MODELO PR/FEATURES_lasso.xls')
my_df.to_excel(writer, 'Hoja de datos')
writer.save()

X_train_lasso = sel_lasso.transform(X_train)
X_test_lasso = sel_lasso.transform(X_test)

                            #AJUSTE DE HIPERPARÁMETROS

  #Definir el rango de valores e hiperparámetros para la optimización

parameters = {'C':[0.0001, 0.001, 0.01, 0.1, 1, 10, 100]}

 #Entrenar el modelo para obtener los mejores hiperparámetros

lr = LogisticRegression(solver='saga')

model_lasso = GridSearchCV(lr, parameters, cv=skf)
model_lasso.fit(X_train_lasso, y_train)
print('Parametros: ', model_lasso.best_params_)
print('Mejor puntaje', model_lasso.best_score_)
print(model_lasso.best_estimator_)

y_pred_lasso = model_lasso.predict(X_test_lasso)
print('Accuracy test', accuracy_score(y_test, y_pred_lasso))
print('Accuracy training', accuracy_score(y_train, model_lasso.predict(X_train_lasso)))

  #Matriz de confusion y reporte de clasificacion de training

cm_train =confusion_matrix(y_train, model_lasso.predict(X_train_lasso))
print(cm_train)

print('Class report train: ', classification_report(y_train, model_lasso.predict(X_train_lasso)))

  #Matriz de confusion y reporte de clasificacion de test

cm =confusion_matrix(y_test, y_pred_lasso)
print(cm)
print(classification_report(y_test, y_pred_lasso))

ax= plt.subplot()
sns.heatmap(cm, annot=True, ax = ax); #annot=True to annotate cells

# labels, title and ticks
ax.set_xlabel('Predicted');ax.set_ylabel('True'); 
ax.set_title('Confusion Matrix'); 
ax.xaxis.set_ticklabels(['Sensible', 'Resistente']); ax.yaxis.set_ticklabels(['Sensible','Resistente'])