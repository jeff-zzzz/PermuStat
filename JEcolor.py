# -*- coding: utf-8 -*-
"""
Created on Fri Apr 24 11:07:29 2015
@author: jesong1126
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

cdict1 = {'red':   ((0.0, 0.0, 0.0),
                   (0.5, 0.0, 0.1),
                   (1.0, 1.0, 1.0)),

         'green': ((0.0, 0.0, 0.0),
                   (1.0, 0.0, 0.0)),

         'blue':  ((0.0, 0.0, 1.0),
                   (0.5, 0.1, 0.0),
                   (1.0, 0.0, 0.0))
        }


blue_red1 = LinearSegmentedColormap('BlueRed1', cdict1)

# Make some illustrative fake data:
x = np.arange(0, np.pi, 0.1)
y = np.arange(0, 2*np.pi, 0.1)
X, Y = np.meshgrid(x,y)
Z = np.cos(X) * np.sin(Y) * 10

# Make the figure:
plt.figure()
plt.imshow(Z, interpolation='nearest', cmap=blue_red1)
plt.colorbar()

plt.show()

"""

Example: suppose you want red to increase from 0 to 1 over the bottom
half, green to do the same over the middle half, and blue over the top
half.  Then you would use:

cdict = {'red':   ((0.0,  0.0, 0.0),
                   (0.5,  1.0, 1.0),
                   (1.0,  1.0, 1.0)),

         'green': ((0.0,  0.0, 0.0),
                   (0.25, 0.0, 0.0),
                   (0.75, 1.0, 1.0),
                   (1.0,  1.0, 1.0)),

         'blue':  ((0.0,  0.0, 0.0),
                   (0.5,  0.0, 0.0),
                   (1.0,  1.0, 1.0))}

If, as in this example, there are no discontinuities in the r, g, and b
components, then it is quite simple: the second and third element of
each tuple, above, is the same--call it "y".  The first element ("x")
defines interpolation intervals over the full range of 0 to 1, and it
must span that whole range.  In other words, the values of x divide the
0-to-1 range into a set of segments, and y gives the end-point color
values for each segment.

Now consider the green. cdict['green'] is saying that for
0 <= x <= 0.25, y is zero; no green.
0.25 < x <= 0.75, y varies linearly from 0 to 1.
x > 0.75, y remains at 1, full green.

If there are discontinuities, then it is a little more complicated.
Label the 3 elements in each row in the cdict entry for a given color as
(x, y0, y1).  Then for values of x between x[i] and x[i+1] the color
value is interpolated between y1[i] and y0[i+1].

Going back to the cookbook example, look at cdict['red']; because y0 !=
y1, it is saying that for x from 0 to 0.5, red increases from 0 to 1,
but then it jumps down, so that for x from 0.5 to 1, red increases from
0.7 to 1.  Green ramps from 0 to 1 as x goes from 0 to 0.5, then jumps
back to 0, and ramps back to 1 as x goes from 0.5 to 1.

row i:   x  y0  y1
                /
               /
row i+1: x  y0  y1

Above is an attempt to show that for x in the range x[i] to x[i+1], the
interpolation is between y1[i] and y0[i+1].  So, y0[0] and y1[-1] are
never used.

"""

##-----------------------------------------------------------------------------

#cmap = plt.get_cmap('BlueRed2')
#plt.imshow(Z, interpolation='nearest', cmap=cmap)
#
## Now we will set the third cmap as the default.  One would
## not normally do this in the middle of a script like this;
## it is done here just to illustrate the method.
#plt.rcParams['image.cmap'] = 'BlueRed3'
#
#plt.subplot(2,2,3)
#plt.imshow(Z, interpolation='nearest')
#plt.colorbar()
#plt.title("Alpha = 1")
#
## Or as yet another variation, we can replace the rcParams
## specification *before* the imshow with the following *after* imshow.
## This sets the new default *and* sets the colormap of the last
## image-like item plotted via pyplot, if any.
#
#plt.subplot(2,2,4)
## Draw a line with low zorder so it will be behind the image.
#plt.plot([0, 10*np.pi], [0, 20*np.pi], color='c', lw=20, zorder=-1)
#
#plt.imshow(Z, interpolation='nearest')
#plt.colorbar()

# Here it is: changing the colormap for the current image and its
# colorbar after they have been plotted.
#plt.set_cmap('BlueRedAlpha')
#plt.title("Varying alpha")
##
#plt.suptitle('Custom Blue-Red colormaps', fontsize=16)


##-----------------------------------------------------------------------------
##
#cmaps = [('Sequential',     ['Blues', 'BuGn', 'BuPu',
#                             'GnBu', 'Greens', 'Greys', 'Oranges', 'OrRd',
#                             'PuBu', 'PuBuGn', 'PuRd', 'Purples', 'RdPu',
#                             'Reds', 'YlGn', 'YlGnBu', 'YlOrBr', 'YlOrRd']),
#         ('Sequential (2)', ['afmhot', 'autumn', 'bone', 'cool', 'copper',
#                             'gist_heat', 'gray', 'hot', 'pink',
#                             'spring', 'summer', 'winter']),
#         ('Diverging',      ['BrBG', 'bwr', 'coolwarm', 'PiYG', 'PRGn', 'PuOr',
#                             'RdBu', 'RdGy', 'RdYlBu', 'RdYlGn', 'Spectral',
#                             'seismic']),
#         ('Qualitative',    ['Accent', 'Dark2', 'Paired', 'Pastel1',
#                             'Pastel2', 'Set1', 'Set2', 'Set3']),
#         ('Miscellaneous',  ['gist_earth', 'terrain', 'ocean', 'gist_stern',
#                             'brg', 'CMRmap', 'cubehelix',
#                             'gnuplot', 'gnuplot2', 'gist_ncar',
#                             'nipy_spectral', 'jet', 'rainbow',
#                             'gist_rainbow', 'hsv', 'flag', 'prism'])]
#
#nrows = max(len(cmap_list) for cmap_category, cmap_list in cmaps)
#gradient = np.linspace(0, 1, 256)
#gradient = np.vstack((gradient, gradient))
#
#def plot_color_gradients(cmap_category, cmap_list):
#    fig, axes = plt.subplots(nrows=nrows)
#    fig.subplots_adjust(top=0.95, bottom=0.01, left=0.2, right=0.99)
#    axes[0].set_title(cmap_category + ' colormaps', fontsize=14)
#
#    for ax, name in zip(axes, cmap_list):
#        ax.imshow(gradient, aspect='auto', cmap=plt.get_cmap(name))
#        pos = list(ax.get_position().bounds)
#        x_text = pos[0] - 0.01
#        y_text = pos[1] + pos[3]/2.
#        fig.text(x_text, y_text, name, va='center', ha='right', fontsize=10)
#
#    # Turn off *all* ticks & spines, not just the ones with colormaps.
#    for ax in axes:
#        ax.set_axis_off()
#
#for cmap_category, cmap_list in cmaps:
#    plot_color_gradients(cmap_category, cmap_list)
#
#plt.show()
#
#xi = np.array([0., 0.5, 1.0])
#yi = np.array([0., 0.5, 1.0])
#zi = np.array([[0., 1.0, 2.0],
#               [0., 1.0, 2.0],
#               [-0.1, 1.0, 2.0]])
#plt.figure()
#v = np.linspace(-.1, 2.0, 15, endpoint=True)
#plt.contour(xi, yi, zi, v, linewidths=0.5, colors='k')
#plt.contourf(xi, yi, zi, v, cmap=plt.cm.jet)
#x = plt.colorbar(ticks=v)
#print x
#plt.show() 
#   
