#Script to visualize strong motions lcoated in a directory using obspy for comparison
import obspy
import os
import fnmatch
import glob
import dyna_convert
import obspy.signal.trigger as triggers
import pandas
import driven_harmonic_oscillator

#directory with subdirectories for gm
path = "/home/mhaas/PhD/IDA/RCFrame_ESDOF/"
#gm_folder = path+"gm/rexelite/C3"

def findFiles(path,pattern):
    '''
    Function to find files with pattern within directory 'path'
    return list of files
    '''
    #check if contains subdirectories
    files = glob.glob(path+'/*')
    for f in files:
        if os.path.isdir(f):
            files=files + glob.glob(f+'/'+pattern)
    #return glob.glob(path+pattern)
    return glob.fnmatch.filter(files,pattern)

#TODO: CHECK LATER WHICH ONE IS SELECTED GIVEN SOIL TYPE+T1
#find ground motions
matches = []
for root, dirnames, filenames in os.walk(path+'gm/rexelite'):
    for filename in fnmatch.filter(filenames, '*README*'):
        matches.append(os.path.join(root, filename))

#store them in a dictionary according to soil type
eigenperiods={}
for fn in matches:
    with open(fn,'r') as f:
        for line in f:
            #get soil type
            if 'site classification' in line:
                soil_type = line.split(':')[-1]
                soil_type = soil_type.strip()
                if soil_type not in eigenperiods.keys():
                    eigenperiods[soil_type]=[]
                    eigenperiods[soil_type+'dir']=[]
            #get maximum value for spectra match T1=max_val/2.
            if 'Period range max [s]' in line:
                period = float(line.split(':')[-1])/2.
                eigenperiods[soil_type].append(period)
                #append fn
                eigenperiods[soil_type+'dir'].append(os.path.dirname(fn))

#response spectrum

#collect pga,duration,soil_type,period
pga = []
duration = []
soil_type = []
period = []

#store each as stream
streams=[]
for s in ['B','C']:
    #all directories containing ground motions for soil type
    dirs = eigenperiods[s+'dir']
    #directories for soil type and periods
    for i,gm_folder in enumerate(dirs):
        #find all gm files in directory
        gmfiles = findFiles(gm_folder,'*ACC.ASC')
        #read the files
        #gms = [eqscale.read_gm(f) for f in gmfiles]
        for f in gmfiles:
            #get as obspy stream
            st = dyna_convert.read_ems_to_stream(f)
            streams.append(st)
            #get trace
            tr = st[0]
            #append pga, duration, soil type and period
            pga.append(tr.max())
            duration.append(tr.stats.endtime-tr.stats.starttime)
            soil_type.append(s)
            period.append(eigenperiods[s][i])

        #for st in streams:
        #    #only one component at a time (both horizontal in different streams)
        #    tr = st[0]
        #    ##get trigger points
        #    #sta = 0.5 #s
        #    #lta = 3 #s
        #    #delta = tr.stats['delta'] #s (sample)
        #    ##recursive sta/lta
        #    #trgs = triggers.classic_sta_lta(tr,int(sta/delta),int(lta/delta))
        #    #triggers.plot_trigger(tr,trgs,1.1,0.05)
        #    print tr.max(),tr.stats.endtime-tr.stats.starttime
        #    st.plot()

#make a pandas dataframe
stats = pandas.DataFrame({'soil':soil_type,'T0':period,'PGA':pga,'t':duration})

#get mean pga (abs) and mean duration per soil type and period
meanPGA  = []
minPGA  = []
maxPGA  = []
meanDur  = []
minDur  = []
maxDur  = []
soil_type = []
period = []
for s in stats.soil.unique():
    for T in stats.T0.unique():
        meanPGA.append(stats[(stats.soil==s)&(stats.T0==T)]['PGA'].abs().mean())
        minPGA.append(stats[(stats.soil==s)&(stats.T0==T)]['PGA'].abs().min())
        maxPGA.append(stats[(stats.soil==s)&(stats.T0==T)]['PGA'].abs().max())
        meanDur.append(stats[(stats.soil==s)&(stats.T0==T)]['t'].mean())
        minDur.append(stats[(stats.soil==s)&(stats.T0==T)]['t'].min())
        maxDur.append(stats[(stats.soil==s)&(stats.T0==T)]['t'].max())
        soil_type.append(s)
        period.append(T)

#make a pandas dataframe
stats_total = pandas.DataFrame({'soil':soil_type,'T0':period,'PGAmax':maxPGA,'PGAmean':meanPGA,'PGAmin':minPGA,'tmax':maxDur,'tmean':meanDur,'tmin':minDur})

#ICONS response
#find a particular trace
date = '19930714'
time = '123148'
component = 'HN2'
station = 'PAT1'

tr = [st for st in streams if st[0].stats['station']==station and st[0].stats['channel']==component and st[0].stats['dyna']['EVENT_DATE_YYYYMMDD']==date and st[0].stats['dyna']['EVENT_TIME_HHMMSS']==time]

#take only first if more than one times in different sets of ground motions
tr = tr[0][0]
dt = 0.005
#go to m/s2
tr = [a/100. for a in tr]

#ICONS ESDOF params
#mass [kg]
m = 91153.0644091
#damping ratio [-]
xi = 0.05
#stiffness [N/m]
k = 3090534.95577
#height for drift
H = 2.7
#
#   #   #test Chopra example 5.2
#   #   #expected ui: 0 0.0437 0.2326 0.6121 1.1116 1.1141 1.1143 1.1143 1.6213 1.9889 2.0947 1.9233 1.5593
#   #   #(succeeded)
#   #   m = 0.2533
#   #   k = 10.
#   #   xi = 0.05
#   #   dt = 0.1
#   #   tr = [0.,5.,8.6602,10.,8.6603,5.0,0.,0.,0.,0.,0.]
#
#
#   #   #test Clough & Penzien 1975 p. 104,106,109
#   #   # expected ui: 0,x,0.0002,x,0.0016,x,0.0054,x,0.0111,x,0.0169
#   #   #(succeeded)
#   #   m = 96.6/32.3 #W/g
#   #   k = 2700. #kips/ft**2
#   #   xi = 0.05
#   #   dt = 0.005
#   #   tr = [0.,19.32,38.64,57.96,77.28,96.60,77.28,57.96,38.64,19.32,0.]
#
#   #   #test Humar p.413-415 (checked only linear acceleration method)
#   #   # expected ui: 0.0000,0.0291,0.2119,0.5896,1.0532,1.3862,1.3644,0.8969,0.1678,-0.5389,-0.9785
#   #   k = 100. #kips/in
#   #   m = 2.533 #kip s2/in
#   #   omega = 6.283 #rad/s
#   #   c = 3.183 #kip s/in
#   #   xi = c/(2.*omega*m)
#   #   dt = 0.1
#   #   tr = [0.,50.,86.6,100.,86.6,50.,0.,0.,0.,0.,0.]
#
Newmark = driven_harmonic_oscillator.Newmark(m,xi,k)
Newmark.set_init_conditions(0,0,gamma=0.5,beta=1./6)
#force
p = [-m*a for a in tr]
deformation,velocity,acceleration = Newmark.groundmotion_response(p,dt)
#   #CDM
#   #CDM = driven_harmonic_oscillator.CentralDifferenceMethod(m,xi,k)
#   #deformation,velocity,acceleration = CDM.groundmotion_response(p,dt,0.,0.)
#   #
#   #       #instance of solver
#   #       Newmark = driven_harmonic_oscillator.Newmark(m,xi,k)
#   #
#   #       #initial conditions
#   #       u_0 = 0;
#   #       ut_0 = 0;
#   #
#   #       #linear acceleration method
#   #       gamma = 1./2
#   #       beta = 1./6
#   #
#   #
#   #       #solve ODE converting cm/s2 to m/s2
#   #       #tr = [a/100. for a in tr]
#   #       deformation,velocity,acceleration = Newmark.groundmotion_response(tr,dt,u_0,ut_0,gamma,beta)
#   #
#   #       #CDM
#   #       #CDM = driven_harmonic_oscillator.CentralDifferenceMethod(m,xi,k)
#   #       #deformation,velocity,acceleration = CDM.groundmotion_response(tr,dt,u_0,ut_0)
#   #
#   #
#   #       #sensitivity study
#   #       import matplotlib.pyplot as plt
#   #
#   #       ms = [0.001,0.005,0.01,0.02,0.05,0.07,0.1,1.0,2.0,5.0,7.,10.]
#   #       xi = 0.05
#   #       k = 2.0
#   #       out = []
#   #       for m in ms:
#   #           Newmark = driven_harmonic_oscillator.Newmark(m,xi,k)
#   #           deformation,velocity,acceleration = Newmark.groundmotion_response(tr,dt,u_0,ut_0,gamma,beta)
#   #           out.append(max([abs(d) for d in deformation]))
#   #       plt.plot(ms,out)
#   #       plt.show()
#   #       ##2) stiffness
#   #       m  = 1.
#   #       xi = 0.05
#   #       ks = [0.001,0.005,0.01,0.02,0.05,0.07,0.1,1.0,2.0,5.0,7.,10.]
#   #       out = []
#   #       for k in ks:
#   #           Newmark = driven_harmonic_oscillator.Newmark(m,xi,k)
#   #           deformation,velocity,acceleration = Newmark.groundmotion_response(tr,dt,u_0,ut_0,gamma,beta)
#   #           out.append(max([abs(d) for d in deformation]))
#   #       plt.plot(ks,out)
#   #       plt.show()
#   #
#   #       #damping
#   #       m  = 1.
#   #       xis = [0.01,0.02,0.03,0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0]
#   #       k = 7.0
#   #       out = []
#   #       for xi in xis:
#   #           Newmark = driven_harmonic_oscillator.Newmark(m,xi,k)
#   #           deformation,velocity,acceleration = Newmark.groundmotion_response(tr,dt,u_0,ut_0,gamma,beta)
#   #           out.append(max([abs(d) for d in deformation]))
#   #       plt.plot(xis,out)
#   #       plt.show()
#   #
#   #
#   #generate a plot
import matplotlib.pyplot as plt

times = [i*dt for i in range(len(tr))]
plt.plot(times,tr)
plt.show()

plt.plot(times,deformation)
plt.show()
#   #
#   #write ground motion file in [g]
#   #g = 9.81
#   #with open('gm.csv','w') as f:
#   #    for i in range(len(tr)):
#   #        f.write('{},{}\n'.format(times[i],tr[i]/9.81))
#   #OpenSees style gm m/s2
#   acc = tr
#   print 'timesteps:',len(acc)
#   with open('gm.acc','w') as outfile:
#       i=1
#       while i < len(acc)+1:
#           #first element
#           if i==1:
#               onerow = [acc[i-1]]
#           else:
#               onerow.append(acc[i-1])
#               #if full row print to file
#               if i%5==0:
#                   #do this in order to get two whitespace if positive and one if negative
#                   strings = ['{:8.7e}'.format(f) if f<0 else ' {:8.7e}'.format(f) for f in onerow]
#                   #write to file whitespace separated
#                   outfile.write(' {:9} {:9} {:9} {:9} {:9}\n'.format(*strings))
#                   #reset collector
#                   onerow = []
#           #next element
#           i+=1
#
#       #collect last elements (already collected above) and write row (no matter if 5 or less elements
#       #create strings
#       strings = [' {:8.7e}'.format(f) if f<0 else '  {:8.7e}'.format(f) for f in onerow]
#       #join it and write to file
#       outfile.write(''.join(strings))
#
