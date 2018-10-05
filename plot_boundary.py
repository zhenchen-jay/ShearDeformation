import numpy
import matplotlib.pyplot as plt

def plot_boundary():
    boundary = numpy.loadtxt('/Users/chenzhen/UT/Research/Projects/ShearDeformation/build/Release/boundary.txt')
    boundary = numpy.delete(boundary, 2, axis=1)
    plt.figure()
    plt.plot(boundary[:,0],boundary[:,1])
    plt.show()


plot_boundary()


