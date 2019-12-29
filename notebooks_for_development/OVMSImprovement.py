#import urllib2
import time
import math
import datetime
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

FileName=['1','2']
#root='.'

nmfile='./data/PLC_2019_04_20T04_33.txt'



#load data into an array
#fn=root+nmfile
fn=nmfile

first=1
array=np.loadtxt(fn,delimiter=',')
fullarray=array

#extract columns
JDtime=fullarray[:,8]
time=(JDtime-JDtime[0])*3600*24

Closed=fullarray[:,20]
Ksp=fullarray[:,28]
Ksp0=fullarray[:,29]
#Kph=fullarray[:,29]
Kunph=fullarray[:,31] #listed as 32  (-1)
#Hsp=fullarray[:,33]
#Hsp0=fullarray[:,34]
#Hph=fullarray[:,35]
Hunph=fullarray[:,37]
tip=fullarray[:,64] / 10.0 #65
tilt=fullarray[:,72] / 10.0  #73
SNR=fullarray[:,75]  #76
PL=fullarray[:,81]
tipcmd=fullarray[:,82]
tiltcmd=fullarray[:,83]
v1=fullarray[:,84]
v2=fullarray[:,85]
v3=fullarray[:,86]
v4=fullarray[:,87]
v5=fullarray[:,88]
ovmsgain=fullarray[:,88] #89 OVMS gain
ovmspl=fullarray[:,89]

#plt.scatter(time,ovmsgain)
#plt.show()
#ovms_status = pd.read_csv(fn)
#print(ovms_status.keys())

print(np.where(ovmsgain == 0))

v8=fullarray[:,91]

pathK=(Kunph-Ksp)*2.2/360
pathH=(Hunph)*1.65/360
diffHK=pathK-pathH

rms = np.std(pathK )
print('Overall RMS  '+ str(rms) +' microns')
numlines =len(pathK)

startopenindex = int(25.0/60.0*numlines)
endopenindex = int(33.0/60.0*numlines)
startclosedindex = int(35.0/60.0 * numlines)

pathKopen=pathK[startopenindex:endopenindex]
pathKclosed=pathK[startclosedindex:numlines]

rmsopen = np.std(pathKopen)
rmsclosed = np.std(pathKclosed)

print("Open loop RMS is "+str(rmsopen)+" microns")
print("Closed loop RMS is "+str(rmsclosed)+" microns")




plt.figure(1,figsize=(20,12))
plt.rcParams.update({'font.size': 15})
plt.plot(np.subtract(time,30), pathK, label='Pathlength Error (um)')
#plt.plot(time, ovmspl, label='OVMS OPD')
plt.axvline(x=np.subtract(time[30863],30), linestyle="--", color="k") # where OVMS is turned on
#plt.plot(time, ovmsgain, label='OVMS gain')
#plt.plot(time, PL,label='Commanded Pathlength')
#plt.plot(time, pathH, label='H band Pathlength Error (um)')
#plt.plot(time, diffHK,label='difference')         
#plt.plot(time, Ksp/100, label='K setpoint')
#plt.title(nmfile);
#plt.plot(time, Closed, '+', label='Closed')
#plt.plot(time,tip-10.0,label='Tip Error/10 -5')
#plt.plot(time,tilt-20.0,label='Tilt Error/10 -20')
#plt.plot(time,tipcmd*10.0,'o',label='Corrected Tip*10')
#plt.plot(time,tiltcmd*10.0,'o',label='Corrected Tilt*10')
#plt.plot(time, SNR*3.0, label='SNR*3')
plt.ylim((-3,3))
plt.xlim((0,8))
#plt.plot(time, Hunphase-Ksp, label='H phase ')
#plt.plot(time, Kavunph-Ksp[0], label='K unwr. av. ')
#plt.plot(time, Havunph-Ksp[0], label='H unwr. av. ')
#plt.plot(time, pseudoKav, label='pseudoK av')
#plt.plot(time, wt2, label='water term 2')
#plt.plot(time, watert, label='water term')

#plt.plot(time, Hsp, label='H setpoint')
plt.xlabel('Elapsed time (sec)')
plt.ylabel('Pathlength error ($\mu$m)')
plt.annotate("Phase loop on,\nOVMS off", [1.4,-2.5])
plt.annotate("Phase loop on,\nOVMS on", [5.4,-2.5])
#plt.legend(loc='lower right')
#plt.xlim((start,stop))
#plt.ylim((-360,720))
#plt.savefig('lowFreqFollowing1237.png')
#plt.savefig('OVMSonoff.png')
plt.show()


plt.figure(2, figsize=(20,12))
plt.plot(time,ovmsgain)
#plt.savefig('OVMSgain.png')
#plt.plot(time, v1, label='v1 ')
plt.show()


