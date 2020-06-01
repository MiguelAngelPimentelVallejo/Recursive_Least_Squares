#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""  This script is to do a recursive least square method comparation """

__author__ = '{Miguel Angel Pimentel Vallejo}'
__email__ = '{miguel.pimentel@umich.mx}'
__date__= '{18/may/2020}'

import numpy as np
import matplotlib.pyplot as plt

#Function with the algorithm of recursive least squares
def R_Lst_Sqr(phi,y,omega_past,P_past):

    for i in range( omega_past.shape[0] ):
        phi = np.concatenate((phi,[phi[0]**i]),axis=0)

    phi = np.array([phi[1:,0]]).T

    K = P_past.dot( phi.dot( np.linalg.inv( np.identity( phi.shape[1] ) + phi.T.dot( P_past.dot( phi ) ) ) ) ) 
    P = P_past.dot( np.identity( phi.shape[0] ) - K.dot( phi.T ) )
    omega = omega_past +  K.dot( y - phi.T.dot( omega_past ) )
    yhat = phi.T.dot(omega_past)

    #print( y - phi.T.dot( omega_past ))

    return yhat,P,omega
        

#Generate al polynomial signal with noise
x = np.linspace(0,50,10000)
coef = np.array([[1.1,0.45,0.1]])
y = coef[0,0] + coef[0,1]*x + coef[0,2]*x**2 + np.random.normal(0.1, 0.5, len(x))

"Zero grade approximation"

#Select the grade of the aproximation
grade = 1

#Caculate the intial values for P matrix and omega vector
zeta = -1000 #0.01
P = zeta*np.identity( grade )
omega = np.zeros((grade,1))

#Arrays to save the results
yhat_v = np.array([])
omega_v = np.empty([0,grade])

for i in range(grade,x.shape[0]):
    yhat,P,omega = R_Lst_Sqr(np.array([[x[i]]]),np.array([[y[i]]]),omega,P)
    yhat_v = np.concatenate((yhat_v,yhat[0]))
    omega_v = np.concatenate((omega_v,[omega.flatten()]),axis=0)

#Plot results
plt.figure()
plt.title("Recursive Least Squares")
plt.xlabel("x")
plt.ylabel("magnitude")
plt.plot(x,y,label = "$y$")
plt.plot(x[grade:],yhat_v,'--', label = "$\^y$" )
plt.grid(True)
plt.legend()

plt.figure()
plt.title("Dynamic's Coeficients")
plt.xlabel("x")
plt.ylabel("magnitude")
labels = ['$\\^\\theta_1$','$\\theta_1$']
omega_l =plt.plot(x,(np.ones(len(x))*coef.T).T)
omegahat_l = plt.plot(x[grade:],omega_v, '--')
plt.grid(True)
plt.legend(omega_l + omegahat_l, labels, loc = 4)

plt.figure()
plt.title("Error")
plt.xlabel("x")
plt.ylabel("magnitude")
plt.plot(x[grade:],y[grade:]-yhat_v,label="$e$" )
plt.grid(True)
plt.legend()


"First grade approximation"

#Select the grade of the aproximation
grade = 2

#Caculate the intial values for P matrix and omega vector
zeta = 0.01
P = zeta*np.identity( grade )
omega = np.zeros((grade,1))

#Arrays to save the results
yhat_v = np.array([])
omega_v = np.empty([0,grade])

for i in range(grade,x.shape[0]):
    yhat,P,omega = R_Lst_Sqr(np.array([[x[i]]]),np.array([[y[i]]]),omega,P)
    yhat_v = np.concatenate((yhat_v,yhat[0]))
    omega_v = np.concatenate((omega_v,[omega.flatten()]),axis=0)

#Plot results
plt.figure()
plt.title("Recursive Least Squares")
plt.xlabel("x")
plt.ylabel("magnitude")
plt.plot(x,y,label = "$y$")
plt.plot(x[grade:],yhat_v,'--', label = "$\^y$" )
plt.grid(True)
plt.legend()

plt.figure()
plt.title("Dynamic's Coeficients")
plt.xlabel("x")
plt.ylabel("magnitude")
labels = ['$\\^\\theta_1$','$\\^\\theta_2$','$\\theta_1$','$\\theta_2$']
omega_l =plt.plot(x,(np.ones(len(x))*coef.T).T)
omegahat_l = plt.plot(x[grade:],omega_v, '--')
plt.grid(True)
plt.legend(omega_l + omegahat_l, labels, loc = 4)

plt.figure()
plt.title("Error")
plt.xlabel("x")
plt.ylabel("magnitude")
plt.plot(x[grade:],y[grade:]-yhat_v,label="$e$" )
plt.grid(True)
plt.legend()

"Second grade approximation"

#Select the grade of the aproximation
grade = 3

#Caculate the intial values for P matrix and omega vector
zeta = 0.01
P = zeta*np.identity( grade )
omega = np.zeros((grade,1))

#Arrays to save the results
yhat_v = np.array([])
omega_v = np.empty([0,grade])

for i in range(grade,x.shape[0]):
    yhat,P,omega = R_Lst_Sqr(np.array([[x[i]]]),np.array([[y[i]]]),omega,P)
    yhat_v = np.concatenate((yhat_v,yhat[0]))
    omega_v = np.concatenate((omega_v,[omega.flatten()]),axis=0)

#Plot results
plt.figure()
plt.title("Recursive Least Squares")
plt.xlabel("x")
plt.ylabel("magnitude")
plt.plot(x,y,label = "$y$")
plt.plot(x[grade:],yhat_v,'--', label = "$\^y$" )
plt.grid(True)
plt.legend()

plt.figure()
plt.title("Dynamic's Coeficients")
plt.xlabel("x")
plt.ylabel("magnitude")
labels = ['$\\^\\theta_1$','$\\^\\theta_2$','$\\^\\theta_3$','$\\theta_1$','$\\theta_2$','$\\theta_3$']
omega_l =plt.plot(x,(np.ones(len(x))*coef.T).T)
omegahat_l = plt.plot(x[grade:],omega_v, '--')
plt.grid(True)
plt.legend(omega_l + omegahat_l, labels, loc = 4)

plt.figure()
plt.title("Error")
plt.xlabel("x")
plt.ylabel("magnitude")
plt.plot(x[grade:],y[grade:]-yhat_v,label="$e$" )
plt.grid(True)
plt.legend()


"Third grade approximation"

#Select the grade of the aproximation
grade = 4

#Caculate the intial values for P matrix and omega vector
zeta = 0.01
P = zeta*np.identity( grade )
omega = np.zeros((grade,1))

#Arrays to save the results
yhat_v = np.array([])
omega_v = np.empty([0,grade])

for i in range(grade,x.shape[0]):
    yhat,P,omega = R_Lst_Sqr(np.array([[x[i]]]),np.array([[y[i]]]),omega,P)
    yhat_v = np.concatenate((yhat_v,yhat[0]))
    omega_v = np.concatenate((omega_v,[omega.flatten()]),axis=0)

#Plot results
plt.figure()
plt.title("Recursive Least Squares")
plt.xlabel("x")
plt.ylabel("magnitude")
plt.plot(x,y,label = "$y$")
plt.plot(x[grade:],yhat_v,'--', label = "$\^y$" )
plt.grid(True)
plt.legend()

plt.figure()
plt.title("Dynamic's Coeficients")
plt.xlabel("x")
plt.ylabel("magnitude")
labels = ['$\\^\\theta_1$','$\\^\\theta_2$','$\\^\\theta_3$','$\\^\\theta_4$','$\\theta_1$','$\\theta_2$','$\\theta_3$','$\\theta_4$']
omega_l =plt.plot(x,(np.ones(len(x))*coef.T).T)
omegahat_l = plt.plot(x[grade:],omega_v, '--')
plt.grid(True)
plt.legend(omega_l + omegahat_l, labels, loc = 4)

plt.figure()
plt.title("Error")
plt.xlabel("x")
plt.ylabel("magnitude")
plt.plot(x[grade:],y[grade:]-yhat_v,label="$e$" )
plt.grid(True)
plt.legend()

plt.show()