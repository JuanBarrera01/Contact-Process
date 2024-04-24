import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import os
import numpy as np
import re
import math


dir_path= r'C:\Users\Juan Felipe Barrera\Documents\College Northwestern\Kovacs Lab\DataWQ24\__Condensed\NPY'
graphing =np.zeros((10,10))
fig = plt.figure(figsize = (8,8))
colors=list(mcolors.TABLEAU_COLORS)
colorcount=0

#choose lambda
for path in os.listdir(dir_path):
    print(path)
    begin=re.search('a-3.', path).span()
    end=re.search('-Start',path).span()
    lam = path[begin[0]+2:end[0]]
    if os.path.isdir(os.path.join(dir_path, path)):
        dir=os.path.join(dir_path, path)
        counter=0
        #cycle through all files for lambda
        for subpath in os.listdir(dir):
            if os.path.isfile(os.path.join(dir,subpath)):
                convert=np.load(os.path.join(dir,subpath), allow_pickle=True)
                length= int(convert.size/3)
                if counter==0:
                    graphing = np.zeros((3, length))
                    graphing[0]=convert[0]
                    #aligned time indices
                graphing[1]+=convert[1]
                graphing[2]+=convert[2]
                #sum how many survivors (NS)
                #for i in range(length):
                    #if convert[1,i]!=0:
                        #graphing [2,i]+=1
                        #sum if survived (PS)
                #sum how many survivors (NS)
                counter+=1
        print(counter)
        graphing[1]=graphing[1]/counter 
        #graphing 1 is NS
        graphing[2]=graphing[2]/counter
        #Graphing 2 is PS
        #divided both counters to get real NS and PS
        #Plotting Ns vs PS
        plt.plot(np.log(graphing[0]), graphing[1] , c = colors[colorcount%10] , marker = ".", linestyle = "None", label = lam, markersize = 3)
        colorcount+=1
plt.legend(fontsize = 12, markerscale = 2)
plt.grid(alpha = 0.5, linestyle = "--")
plt.xlabel("log(t)", size = 12)
plt.ylabel("Ns", size = 12)
plt.xscale("log")
plt.yscale("log")
fig.savefig(r"C:\Users\Juan Felipe Barrera\Documents\College Northwestern\Kovacs Lab\DataWQ24\__Condensed\NSvslogt04-14.jpg", transparent = True, bbox_inches = 'tight', pad_inches = 0)
fig.savefig(r"C:\Users\Juan Felipe Barrera\Documents\College Northwestern\Kovacs Lab\DataWQ24\__Condensed\NSvslogt04-14.pdf", transparent = True, bbox_inches = 'tight', pad_inches = 0)