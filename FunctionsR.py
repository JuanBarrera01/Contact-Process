import numpy as np
import networkx as nx
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import scipy as sp
from numpy import linalg
from numpy.random import rand
from os import times
import matplotlib.colors as mcolors
import time
import json
import math
start_time = time.time()

#Creating a lattice with periodic boundaries
def create_lattice(size, perc, start):
    #latticetime=time.time()
    lattice=np.zeros((size**2,4), dtype=bool)
    #use order right, down, left, up
    for j in range(size**2):
        randomizer = np.random.rand(2)
        #bonds with right neighbour
        #randomizer removes bond with probability p
        if j % size != (size-1) and randomizer[0]>=perc:
            lattice[j, 0]=lattice[j+1, 2]=1
        #bonds horizontal periodic
        elif randomizer[0]>=perc:
            lattice[j, 0]=lattice[j-(size-1), 2]=1
        if j < (size*(size-1)) and randomizer[1]>=perc:
            lattice[j, 1]=lattice[j+size,3]=1
        elif randomizer[1]>=perc:
            lattice[j,1]=lattice[j-size*(size-1),3]=1
    #All infected matrix
    if start == 'all':
        inflist=np.arange(size**2)
        inflist = np.insert(inflist, 0, 2*size**2)
    else:
        inflist=np.array([2*size**2, np.random.randint(0,size**2)])
        
    #print ("Lattice built in %s"  %(time.time()-latticetime))
    return lattice, inflist
#Run on MC Step
def step(inflist, infrate, lattice, size, trace, t, wholetime, track, NR, NR2):
    check =False
    if (inflist.size > 1):          
        #choose random infected node
        if (inflist.size !=2):
            i=np.random.randint(1, inflist.size)
        else:
            i=1
        randomizer= np.random.rand(1)
        #determining infection or heal
        if  1/(1+infrate)>randomizer:
            NR-= np.sqrt((inflist [i]%size)**2+(inflist[i]//size)**2)
            NR2 -= (inflist [i]%size)**2+(inflist[i]//size)**2
            inflist[i]=inflist[-1]
            inflist=inflist[:-1]
        else:
            site=inflist[i]
            infection = np.random.randint(0,4)
            neighbours= [size*(site//size)+(site+1)%size,(site+size)%(size**2),size*(site//size)+(site-1)%size,(site-size)%(size**2)]
            if lattice[site][infection]== 1 and neighbours[infection] not in inflist:
                inflist=np.append(inflist, neighbours[infection])
                NR+= np.sqrt((neighbours[infection]%size)**2+(neighbours[infection]//size)**2)
                NR2+= (neighbours[infection]%size)**2+(neighbours[infection]//size)**2
        t+=1/inflist.size
        #Adding next measurement
        if t >=wholetime :
            trace[1][track] = (inflist.size-1)
            trace[2][track] = NR
            trace[3][track] = NR2
            track +=1
            if track < trace[0].size:
                wholetime = trace[0][track]
    else:
        check = True
        print ('Dead')
    return inflist, trace, t, check, wholetime, track, NR, NR2
#Run one Sim
def cycle(inflist, infrate, lattice, size, duration):

    maxlog = math.ceil(math.log(duration, 1.05))
    temp = np.arange(maxlog+1)
    temp =np.floor(1.05**temp)
    temp = np.insert(temp, 0, 0)
    temp [-1]=duration
    temp = np.unique(temp)
    trace = np.zeros((4, temp.size))
    #establish timescale
    trace[0]=temp  
    #establish infected size at beggining
    trace[1][0]= (inflist.size-1)
    #establish N<R> (only works for start one)
    trace[2][0]=math.sqrt((inflist[1]%size)**2+(inflist[1]//size)**2)
    NR=math.sqrt((inflist[1]%size)**2+(inflist[1]//size)**2)
    #establish N<R2>
    trace[3][0]=(inflist[1]%size)**2+(inflist[1]//size)**2
    NR2=(inflist[1]%size)**2+(inflist[1]//size)**2
    
    track =1
    wholetime = trace[0][track]
    t=0
    while t < duration:
        inflist, trace, t, stopper, wholetime, track, NR, NR2 =step(inflist, infrate, lattice, size , trace, t, wholetime, track, NR, NR2)
        if stopper == True:
            print('stopped {} at {}'.format(infrate, wholetime))
            break
    return  trace, inflist
#continue sims
def continues(inflist, infrate, lattice, size, duration,t, track, trace):

    wholetime = trace[0][track]
    NR = trace[2][track-1]
    NR2 = trace[3][track-1]
    while t < duration:
        inflist, trace, t, stopper, wholetime, track, NR, NR2 =step(inflist, infrate, lattice, size , trace, t, wholetime, track, NR, NR2)
        if stopper == True:
            print('stopped {} at {}'.format(infrate, wholetime))
            break
    return  trace, inflist
#store all sims
#store all sims
def simulation(size, perc, infrate, duration, numruns):
    s = np.zeros((numruns, duration+1))
    for trial in range(numruns):
        print(trial)
        lattice, inflist= create_lattice(size, perc)
        s[trial] = cycle(inflist, infrate, lattice, size, duration)
    return s
#store only surviving sims
def survivalsim(size, perc, infrate, duration, numsurv):
    s = np.zeros((numsurv, duration+1))
    survivors=0
    attempts=0
    while survivors < numsurv:
        print( "Attempt {} at {}".format(attempts,infrate))
        lattice, inflist = create_lattice(size, perc)
        trace = cycle(inflist, infrate, lattice, size, duration)
        #print (trace[-1])
        if trace[-1]!=0:
            s[survivors]=trace
            survivors+=1  
            print('survivor {} at {}'.format(survivors, infrate))
        attempts +=1  
    return s

#ignore everything from here:
#ploting funcitions:
def plotmany(manyruns, size, infrate, duration, perc):
    fig = plt.figure(figsize = (5,5))
    colors=list(mcolors.TABLEAU_COLORS)
    for trial in range (len(manyruns)):
        plt.plot( manyruns[trial] , c = colors[trial%10], marker = ".", linestyle = "None", markersize = 3)
    plt.xlabel("t", size = 12)
    plt.ylabel("density", size = 12)
    plt.xscale("log")
    plt.yscale("log")
    plt.legend(fontsize = 12, markerscale = 2)
    plt.grid(alpha = 0.5, linestyle = "--")
    plt.title("Size " + str(size) + "p= " + str(perc) + "lambda " + str(infrate))
    fig.savefig("Graphs/Multirun/" + f"Size  {size}  p={perc} lambda {infrate} steps {duration}  {len(manyruns)} runs.jpg", transparent = True, bbox_inches = 'tight', pad_inches = 0)
    fig.savefig("Graphs/Multirun/" + f"Size  {size}  p={perc} lambda {infrate} steps {duration}  {len(manyruns)} runs.pdf", transparent = True, bbox_inches = 'tight', pad_inches = 0)
def plotav(manyruns,size, infrate, duration, perc):
    average=np.average(manyruns, 0)   
    fig = plt.figure(figsize = (5,5))
    plt.plot( average , c = 'red' , marker = ".", linestyle = "None", label = infrate, markersize = 3)
    plt.xlabel("t", size = 12)
    plt.ylabel("density", size = 12)
    plt.xscale("log")
    plt.yscale("log")
    plt.legend(fontsize = 12, markerscale = 2)
    plt.grid(alpha = 0.5, linestyle = "--")
    plt.title("Size " + str(size) + "p= " + str(perc) + "lambda " + str(infrate) + " " + " average of " + str(len(manyruns)) + "  runs")
    fig.savefig(  f"Graphs/Averages/ Size {size}  p= {perc}  lambda {infrate} steps {duration} average {len(manyruns)} runs.jpg", transparent = True, bbox_inches = 'tight', pad_inches = 0)
    fig.savefig(  f"Graphs/Averages/ Size {size}  p= {perc}  lambda {infrate} steps {duration} average {len(manyruns)} runs.pdf", transparent = True, bbox_inches = 'tight', pad_inches = 0)
def plotavsurv(manyruns, size, infrate, duration,perc):
    surviving=0
    average = np.zeros(duration+1)
    for j in range(len(manyruns)):
        #only average surviving 
        if manyruns[j][-1]!=0:
            average+=manyruns[j]
            surviving+=1
            #print('survivor #' +str(surviving))
    if surviving != 0: average=average/surviving
    fig = plt.figure(figsize = (5,5))
    plt.plot( average , c = 'red' , marker = ".", linestyle = "None", label = infrate, markersize = 3)
    plt.xlabel("t", size = 12)
    plt.ylabel("density", size = 12)
    plt.xscale("log")
    plt.yscale("log")
    plt.legend(fontsize = 12, markerscale = 2)
    plt.grid(alpha = 0.5, linestyle = "--")
    plt.title("Size " + str(size) + "p= " + str(perc) + "lambda " + str(infrate) + " " + str(surviving) + " surviving runs")
    fig.savefig( f"Graphs/Survivorav/ Size {size} p= {perc} lambda {infrate} steps {duration} average {surviving} surviving runs.jpg", transparent = True, bbox_inches = 'tight', pad_inches = 0)
    fig.savefig( f"Graphs/Survivorav/ Size {size} p= {perc} lambda {infrate} steps {duration} average {surviving} surviving runs.pdf", transparent = True, bbox_inches = 'tight', pad_inches = 0)
def plotsurv(manyruns, size, infrate, duration, perc):
    colors=list(mcolors.TABLEAU_COLORS)
    surviving=0
    fig = plt.figure(figsize = (5,5))
    for j in range(len(manyruns)):
        #only plot surviving 
        if manyruns[j][-1]!=0:
            surviving+=1
            plt.plot( manyruns[j] , c = colors[j%10], marker = ".", linestyle = "None", markersize = 3)
    plt.xlabel("t", size = 12)
    plt.ylabel("density", size = 12)
    plt.xscale("log")
    plt.yscale("log")
    plt.legend(fontsize = 12, markerscale = 2)
    plt.grid(alpha = 0.5, linestyle = "--")
    plt.title("Size " + str(size) + "p= " + str(perc) + "lambda " + str(infrate) + " " + str(surviving) + " surviving runs")
    fig.savefig("Graphs/Survivors/" + "Size " + str(size) + " p= " + str(perc) + "lambda " + str(infrate) + "steps " + str(duration) + "  " + str(surviving) + " surviving runs.jpg", transparent = True, bbox_inches = 'tight', pad_inches = 0)
    fig.savefig( "Graphs/Survivors/"+ "Size " + str(size) + " p= " + str(perc) + "lambda " + str(infrate) + "steps " + str(duration) + "  " + str(surviving) + " surviving runs.pdf", transparent = True, bbox_inches = 'tight', pad_inches = 0)


def plotmanyone(manyruns, size, infrate, duration, perc):
    fig = plt.figure(figsize = (5,5))
    colors=list(mcolors.TABLEAU_COLORS)
    for trial in range (len(manyruns)):
        plt.plot( manyruns[trial] , c = colors[trial%10], marker = ".", linestyle = "None", markersize = 3)
    plt.xlabel("t", size = 12)
    plt.ylabel("density", size = 12)
    plt.xscale("log")
    plt.yscale("log")
    plt.legend(fontsize = 12, markerscale = 2)
    plt.grid(alpha = 0.5, linestyle = "--")
    plt.title("Size " + str(size) + "p= " + str(perc) + "lambda " + str(infrate))
    fig.savefig("Graphs/Multirunone/" + f"Size  {size}  p={perc} lambda {infrate} steps {duration}  {len(manyruns)} runs.jpg", transparent = True, bbox_inches = 'tight', pad_inches = 0)
    fig.savefig("Graphs/Multirunone/" + f"Size  {size}  p={perc} lambda {infrate} steps {duration}  {len(manyruns)} runs.pdf", transparent = True, bbox_inches = 'tight', pad_inches = 0)
def plotavone(manyruns,size, infrate, duration, perc):
    average=np.average(manyruns, 0)   
    logt = np.arange(len(average))
    logt = np.log10(logt)
    fig = plt.figure(figsize = (5,5))
    plt.plot(logt, average , c = 'red' , marker = ".", linestyle = "None", label = infrate, markersize = 3)
    plt.xlabel("log t", size = 12)
    plt.ylabel("density", size = 12)
    plt.xscale("log")
    plt.yscale("log")
    plt.legend(fontsize = 12, markerscale = 2)
    plt.grid(alpha = 0.5, linestyle = "--")
    plt.title("Size " + str(size) + "p= " + str(perc) + "lambda " + str(infrate) + " " + " average of " + str(len(manyruns)) + "  runs")
    fig.savefig(  f"Graphs/Averagesone/ Size {size}  p= {perc}  lambda {infrate} steps {duration} average {len(manyruns)} runs.jpg", transparent = True, bbox_inches = 'tight', pad_inches = 0)
    fig.savefig(  f"Graphs/Averagesone/ Size {size}  p= {perc}  lambda {infrate} steps {duration} average {len(manyruns)} runs.pdf", transparent = True, bbox_inches = 'tight', pad_inches = 0)
def plotavsurvone(manyruns, size, infrate, duration,perc):
    surviving=0
    average = np.zeros(duration+1)
    for j in range(len(manyruns)):
        #only average surviving 
        if manyruns[j][-1]!=0:
            average+=manyruns[j]
            surviving+=1
            #print('survivor #' +str(surviving))
    if surviving != 0: average=average/surviving
    fig = plt.figure(figsize = (5,5))
    plt.plot( average , c = 'red' , marker = ".", linestyle = "None", label = infrate, markersize = 3)
    plt.xlabel("t", size = 12)
    plt.ylabel("density", size = 12)
    plt.xscale("log")
    plt.yscale("log")
    plt.legend(fontsize = 12, markerscale = 2)
    plt.grid(alpha = 0.5, linestyle = "--")
    plt.title("Size " + str(size) + "p= " + str(perc) + "lambda " + str(infrate) + " " + str(surviving) + " surviving runs")
    fig.savefig( f"Graphs/Survivoravone/ Size {size} p= {perc} lambda {infrate} steps {duration} average {surviving} surviving runs.jpg", transparent = True, bbox_inches = 'tight', pad_inches = 0)
    fig.savefig( f"Graphs/Survivoravone/ Size {size} p= {perc} lambda {infrate} steps {duration} average {surviving} surviving runs.pdf", transparent = True, bbox_inches = 'tight', pad_inches = 0)
def plotsurvone(manyruns, size, infrate, duration, perc):
    colors=list(mcolors.TABLEAU_COLORS)
    surviving=0
    fig = plt.figure(figsize = (5,5))
    for j in range(len(manyruns)):
        #only plot surviving 
        if manyruns[j][-1]!=0:
            surviving+=1
            plt.plot( manyruns[j] , c = colors[j%10], marker = ".", linestyle = "None", markersize = 3)
    plt.xlabel("t", size = 12)
    plt.ylabel("density", size = 12)
    plt.xscale("log")
    plt.yscale("log")
    plt.legend(fontsize = 12, markerscale = 2)
    plt.grid(alpha = 0.5, linestyle = "--")
    plt.title("Size " + str(size) + "p= " + str(perc) + "lambda " + str(infrate) + " " + str(surviving) + " surviving runs")
    fig.savefig("Graphs/Survivorsone/" + "Size " + str(size) + " p= " + str(perc) + "lambda " + str(infrate) + "steps " + str(duration) + "  " + str(surviving) + " surviving runs.jpg", transparent = True, bbox_inches = 'tight', pad_inches = 0)
    fig.savefig( "Graphs/Survivorsone/"+ "Size " + str(size) + " p= " + str(perc) + "lambda " + str(infrate) + "steps " + str(duration) + "  " + str(surviving) + " surviving runs.pdf", transparent = True, bbox_inches = 'tight', pad_inches = 0)

def plotdensityandsize(manyruns, size, infrate, duration, perc):
    numsurv=np.average(manyruns, 0)
    probsurv=np.count_nonzero(manyruns, 0)/len(numsurv)
    fig = plt.figure(figsize = (5,5))
    plt.plot(probsurv, numsurv , c = 'red' , marker = ".", linestyle = "None", label = infrate, markersize = 3)
    plt.xlabel("Ps", size = 12)
    plt.ylabel("Ns", size = 12)
    plt.xscale("log")
    plt.yscale("log")
    plt.legend(fontsize = 12, markerscale = 2)
    plt.grid(alpha = 0.5, linestyle = "--")
    plt.title("Size " + str(size) + "p= " + str(perc) + "lambda " + str(infrate) + " " + " average of " + str(len(manyruns)) + "  runs")
    fig.savefig(  f"Graphs/Density-and-Probability/ Size {size}  p= {perc}  lambda {infrate} steps {duration} average {len(manyruns)} runs.jpg", transparent = True, bbox_inches = 'tight', pad_inches = 0)
    fig.savefig(  f"Graphs/Density-and-Probability/ Size {size}  p= {perc}  lambda {infrate} steps {duration} average {len(manyruns)} runs.pdf", transparent = True, bbox_inches = 'tight', pad_inches = 0)


#this one is not optimized, will worry later
def plotlambdas(lamplot, size, infrates, duration):
    fig = plt.figure(figsize = (5,5))
    colors=list(mcolors.TABLEAU_COLORS)
    for trial in range (infrates.size):
        plt.plot( lamplot[trial] , c = colors[trial%10], marker = ".", linestyle = "None", label = "λ = " + str(infrates[trial]),  markersize = 3)
    plt.xlabel("t", size = 12)
    plt.ylabel("density", size = 12)
    plt.xscale("log")
    plt.yscale("log")
    plt.legend(fontsize = 12, markerscale = 2)
    plt.grid(alpha = 0.5, linestyle = "--")
    plt.title("Size " + str(size))
    fig.savefig("Graphs/Lambdas/" + "Size" + str(size) + " Time " + str(duration) +str(infrates) + "Surv " +  " lambdas.jpg", transparent = True, bbox_inches = 'tight', pad_inches = 0 )
    fig.savefig( "Graphs/Lambdas/" + "Size" + str(size) + " Time " + str(duration) +str(infrates) + "Surv " + " lambdas.pdf", transparent = True, bbox_inches = 'tight', pad_inches = 0)

def plotlambdasone(lamplot, size, infrates, duration):
    fig = plt.figure(figsize = (5,5))
    colors=list(mcolors.TABLEAU_COLORS)
    for trial in range (infrates.size):
        plt.plot( lamplot[trial] , c = colors[trial%10], marker = ".", linestyle = "None", label = "λ = " + str(infrates[trial]),  markersize = 3)
    plt.xlabel("t", size = 12)
    plt.ylabel("density", size = 12)
    plt.xscale("log")
    plt.yscale("log")
    plt.legend(fontsize = 12, markerscale = 2)
    plt.grid(alpha = 0.5, linestyle = "--")
    plt.title("Size " + str(size))
    fig.savefig("Graphs/Lambdasone/" + "Size" + str(size) + " Time " + str(duration) +str(infrates) + "Surv " +  " lambdas.jpg", transparent = True, bbox_inches = 'tight', pad_inches = 0 )
    fig.savefig( "Graphs/Lambdasone/" + "Size" + str(size) + " Time " + str(duration) +str(infrates) + "Surv " + " lambdas.pdf", transparent = True, bbox_inches = 'tight', pad_inches = 0)
