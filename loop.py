# Performs a circular orbital fit to data
# looping over one of the 
import numpy as n
import matplotlib.pyplot as plt
from scipy.optimize import leastsq
#import cataclysmic as cv
from matplotlib import mpl,pyplot
import matplotlib.cm as cm
import lmfit as lm
import os as os
import param as op
import sys
import time as ti
plt.ion() # Activate update
reload(op)

def loop(variable,initial,delta,num):
    import param as op
    reload(op)
    runner=initial+n.arange(num)*delta
    chisqr=[]
    for kk in runner: 
        if os.path.isfile(op.file) == False:
            print 'Error: No input file'
            print 'Check that file is in current folder'
            sys.exit()

        data=n.loadtxt(op.file)

        hjd,vel,err=[],[],[]
        if len(data[0])==3:
            for i in data:
                hjd.append(i[0]),vel.append(i[1]),err.append(i[2])
        if len(data[0])==2:
            for i in data:
                hjd.append(i[0]),vel.append(i[1])
        hjd,vel=n.array(hjd),n.array(vel)

        # %%%%%%%%%%%%%%%%%%%%%%%% ORBIT PARAMETERS

        params = lm.Parameters()
        if variable == 'porb':
            params.add('porb', value= kk, vary=False)
        else:
            params.add('porb', value= op.porb, vary=op.fix_porb)
        if variable == 'hjd0':
            params.add('hjd0', value= kk, vary=False)
        else:
            params.add('hjd0',value=op.hjd0, vary=op.fix_hjd0)
        if variable == 'gama':
            params.add('gama', value= kk, vary=False)
        else:
            params.add('gama',value=op.gama, vary=op.fix_gama)   
        if variable == 'k1':
            params.add('k1', value= kk, vary=False)
        else:
            params.add('k1',value=op.k1, vary=op.fix_k1)
        
        params.add('phaseoff', value= op.phaseoff, vary=op.fix_phaseoff) 


        def res_sin3(pars, x, data=None,sigma=op.sigma+n.zeros(len(hjd))):
            porb = pars['porb'].value
            hjd0 = pars['hjd0'].value
            gama = pars['gama'].value
            k1 = pars['k1'].value
            phaseoff=pars['phaseoff'].value
            model=gama+k1*n.sin(2*n.pi*((x-hjd0)/porb+phaseoff))
            if data is None:
                return model
            if sigma is  None:
                return (model - data)
            return (model - data)/sigma

        results=lm.minimize(res_sin3,params, args=(hjd,vel))
        print 'Looping '+str(variable)+': '+str(kk)+'\n'
        lm.report_errors(params,show_correl=False)
        
        print 'DoF = ',results.nfree
        print 'Chi-squared = ',results.chisqr
        print '--------------'
        chisqr.append(results.chisqr)
    plt.figure(2,facecolor='w')
    plt.clf()
    plt.plot(runner,chisqr,'rs',markersize=6)
    plt.plot(runner,chisqr,'k-')
    plt.xlabel(variable)
    plt.ylabel('$\chi^2$')
    plt.show()
loop(op.variable,op.initial,op.delta,op.num)
