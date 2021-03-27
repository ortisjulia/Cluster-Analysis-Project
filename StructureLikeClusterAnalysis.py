# -*- coding: utf-8 -*-

#!/usr/bin/python3
'''
    Title: StructureLikeClusterAnalysis.py
    Date: 2021-02-17
    Author: Julia Ortis Sunyer
    Description:
        This code allows the user to create a structure-like cluster analysis plot with error bars.
    List of functions:
        No user defined functions.
    List of non-standard modules:
        - Numpy: NumPy is the fundamental package for scientific computing in Python. 
          It is a Python library that provides a multidimensional array object, various 
          derived objects (such as masked arrays and matrices), and an assortment of routines
          for fast operations on arrays, including mathematical, logical, shape manipulation,
          sorting, selecting, I/O, discrete Fourier transforms, basic linear algebra, basic
          statistical operations, random simulation and much more.
        - Matplotlib: Matplotlib is a comprehensive library for creating static, animated,
          and interactive visualizations in Python.
    Procedure:
        First I imported all the modules necessary to run the script. The I specify the amount of
        columns I will have depending on the amount of taxa studied. Following that I input the data
        for cluster 1 (C1) and cluster 2 (C2). Then, I specify the locations and the width of the
        bars to create them. After doing this I specify the name of each bar, the X and Y axis,
        the title of the figure and the legend.
    Usage:
        StructureLikeClusterAnalysis.py
    
'''


#Plot results with standard error bars:
#Import the packages I will need
import numpy as np
import matplotlib.pyplot as plt

#The amount of columns I will have (depends on the amount of taxa)
N = 15

#C1 is the cluster 1 data
C1 = (0, 0, 0, 0, 0, 0.28, 0.30, 0.17, 0.22, 0.17, 1, 1, 1, 1, 1)
#C1 standard error
C1Std = (0, 0, 0, 0, 0, 0.050, 0.063, 0.050, 0.050, 0.055, 0, 0, 0, 0, 0)

#C2 is the cluster 2 data
C2 = (1, 1, 1, 1, 1, 0.72, 0.70, 0.83, 0.78, 0.83, 0, 0, 0, 0, 0)
#C2 standard error
C2Std = (0, 0, 0, 0, 0, 0.050, 0.063, 0.050, 0.050, 0.055, 0, 0, 0, 0, 0)

# the x locations for the taxa
ind = np.arange(N)

# the width of the bars
width = 0.35

#To create the plot bars for C1
p1 = plt.bar(ind, C1, width, yerr=C1Std, color='green')

#To create the plot bars for C2
p2 = plt.bar(ind, C2, width, bottom=C1, color='blue')

#To create the actual plot
plt.ylabel('Ancestry') #Y-axis label
plt.xlabel('Taxa') #X-axis label
plt.title('Structure-like cluster analysis') #Title of the plot
plt.xticks(ind, ('8N05240', '8N05890', '8N06612', '8N73248', '8N73604', 'K006', 'K010', 'K011', 'K015', 'K019', 'Lesina_280', 'Lesina_281', 'Lesina_282', 'Lesina_285', 'Lesina_286'), rotation='vertical') #Name of the bars in the X-axis
plt.yticks(np.arange(0, 1.5, 0.2)) #Numbers in the Y-axis
plt.legend((p1[0], p2[0]), ('Cluster 1', 'Cluster 2')) #Legend of the plot
plt.show() #Plot
