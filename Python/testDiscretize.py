from Discretize import tauchenhussey,tauchen
import numpy as np
import math

mu=0.0
rho=0.95
sigma=0.22
numMarkovState=5

w=0.5+rho/4.0
sigmaZ = sigma/math.sqrt(1-rho**2)
baseSigma=w*sigma +(1-w)*sigmaZ

stateList1,transitionMat1=tauchenhussey(numMarkovState,mu,rho,sigma, baseSigma)
stateList2,transitionMat2=tauchen(numMarkovState,mu,rho,sigma)