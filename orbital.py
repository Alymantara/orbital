# Orbital fitting program
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
params.add('porb', value= op.porb, vary=op.fix_porb)
params.add('hjd0',value=op.hjd0, vary=op.fix_hjd0)
params.add('gama', value= op.gama, vary=op.fix_gama)
params.add('k1', value= op.k1, vary=op.fix_k1)
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
print 'Best fit for '+op.object
print '--------------'
lm.report_errors(params,show_correl=False)
print '--------------'
print 'DoF = ',results.nfree
print 'Chi-squared = ',results.chisqr
phase1=[]
for i in hjd:
	tmp=(i-params['hjd0'].value)/params['porb'].value-n.fix((i-params['hjd0'].value)/params['porb'].value)
	if tmp < 0.0:
		tmp = tmp+1.
	phase1.append(tmp)
phase=n.array(phase1)
    
fig=plt.figure(1,facecolor='w')
plt.clf() 
ax1=plt.subplot2grid((4, 1), (0, 0),rowspan=3)
#ffpar=cv.orbital(phase,vel,[70.0,30.0,saved[3][0],porb])
'''plt.step(phase,vel,'kd',where='mid')
plt.step([i+1.0 for i in phase],vel,'kd',where='mid')
plt.step([i-1.0 for i in phase],vel,'k.',where='mid')
plt.step([i+1.0 for i in phase],vel,'k.',where='mid')'''
plt.errorbar(phase+params['phaseoff'].value,vel,yerr=params['k1'].stderr,marker='.',mec='k',ecolor='k',mfc='k',linestyle='None')
plt.errorbar([i+1.0+params['phaseoff'].value for i in phase],vel,yerr=params['k1'].stderr,marker='.',mec='k',ecolor='k',mfc='k',linestyle='None')
plt.errorbar([i-1.0+params['phaseoff'].value for i in phase],vel,yerr=params['k1'].stderr,marker='.',mec='k',ecolor='k',mfc='k',linestyle='None')

#plt.plot(n.arange(0.,1.,0.01),ffpar[0][0]+ffpar[0][1]*n.sin((n.arange(0.,1.,0.01)-0.05)*2*n.pi))
plt.plot([i+params['phaseoff'].value for i in n.arange(-1,3,0.01)],params['gama'].value+params['k1'].value*n.sin(2.0*n.pi*n.arange(-1,3,0.01)),'b-')
plt.axis([0,2,min(vel)*1.3,max(vel)*1.3])
plt.ylabel('Radial Velocity, km s$^{-1}$')

ax1.yaxis.set_major_locator(plt.MaxNLocator(9))
ax1.xaxis.set_major_locator(plt.MaxNLocator(6))
plt.minorticks_on()
plt.setp(ax1.get_xticklabels(), visible=False)
import datetime
plt.suptitle('Object: '+op.object+'\n Filename: '+op.file+'  '+datetime.datetime.fromtimestamp(ti.time()).strftime('%Y-%m-%d %H:%M:%S'))
#plt.plot(cv.res_sin3(params,n.arange(0,1.9,0.01)),n.arange(0,1.9,0.01),'r-')
size=12
plt.xticks(size=size)
plt.yticks(size=size)

ax2=plt.subplot2grid((4, 1), (3, 0))
res=vel-res_sin3(params,hjd)

plt.errorbar(phase+params['phaseoff'].value,res,yerr=params['k1'].stderr,marker='.',mec='k',ecolor='k',mfc='k',linestyle='None')
plt.errorbar([i+1.0+params['phaseoff'].value for i in phase],res,yerr=params['k1'].stderr,marker='.',mec='k',ecolor='k',mfc='k',linestyle='None')
plt.errorbar([i-1.0+params['phaseoff'].value for i in phase],res,yerr=params['k1'].stderr,marker='.',mec='k',ecolor='k',mfc='k',linestyle='None')
ax2.yaxis.set_major_locator(plt.MaxNLocator(5))
ax2.xaxis.set_major_locator(plt.MaxNLocator(4))
plt.minorticks_on()
plt.xlabel('Orbital Phase, $\phi$')
plt.axhline(y=0.0,linestyle='--',color='k')
plt.ylabel('Residuals')
plt.subplots_adjust(hspace = .001)
plt.axis([0,2,min(res)*1.6,max(res)*1.6])

plt.savefig(op.psname+'.'+op.output,dpi=72,facecolor='w',format=op.output)

