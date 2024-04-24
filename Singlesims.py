from Functions import *
import os
import sys

n= int(sys.argv[1])
lam = float(sys.argv[2])
timesteps = int(sys.argv[3])
start = sys.argv[4]
id = sys.argv[5]
p = 0.5


lattice, inflist, = create_lattice(n, p, start)
trace, inflist = cycle(inflist, lam, lattice, n, timesteps)


dir_path = '__Stored-Data/Integer-Count/Size-{}-Time-{}-Lambda-{}-Start-{}-logscale'.format(n, timesteps, lam, start)
if inflist.size != 1:
    np.save('/home/jfl7897/contact-process/Optimize-Lattice/__Ongoing_Runs/__Size-{}-Time-{}-Lambda-{}/__Inflist/Size-{}-Time-{}-Lambda-{}-Start-{}-File-{}-inflist.npy'.format(n, timesteps, lam, n, timesteps, lam, start, int(id)), inflist)
    np.save('/home/jfl7897/contact-process/Optimize-Lattice/__Ongoing_Runs/__Size-{}-Time-{}-Lambda-{}/__Lattice/Size-{}-Time-{}-Lambda-{}-Start-{}-File-{}-lattice.npy'.format(n, timesteps, lam, n, timesteps, lam, start, int(id)), lattice)
# check whether directory already exists
if not os.path.exists(dir_path):
    os.mkdir(dir_path)

np.save('{}/Size-{}-Time-{}-Lambda-{}-Start-{}-File-{}-logscale.npy'.format(dir_path, n, timesteps, lam, start, int(id)), trace)




