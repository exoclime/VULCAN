from argparse import ArgumentParser
import numpy as np
import matplotlib.pylab as plt

# Get the file path from the command line argument
parser = ArgumentParser()
parser.add_argument("-f","--file", dest="filename",type=str,help="nk_file", metavar="FILE")
args = parser.parse_args()

fname = args.filename

data = np.loadtxt(fname,skiprows=5)

wl = data[:,0]
n = data[:,1]
k = data[:,2]

fig = plt.figure()

plt.plot(wl,n,c='red',label='n')
plt.plot(wl,k,c='blue',ls='dashed',label='k')

plt.legend()

#plt.xlim(9,12)

plt.xscale('log')
plt.yscale('log')

plt.show()

