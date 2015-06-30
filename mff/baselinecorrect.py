# -*- coding: utf-8 -*-

#import config 
import numpy as np   
import matplotlib.pyplot as plt

from mff import read_mff_header, read_mff_data, Kmeans #, getEpochInfos, mff_getSummaryInfo
filePath ='/Users/jesong1126/Work/Data/VGT/VGT_8subj_bcr_blc_ave.mff'

hdr = read_mff_header.read_mff_header(filePath)
        
nC = hdr['nChans']
nSamples = hdr['nSamples']
nSamplesPre = hdr['nSamplesPre']
nTrials = hdr['nTrials']
srate = hdr['Fs']
summaryInfo = hdr['orig'] 
trialsName = summaryInfo['epochLabels']   
categoryName = list(set(trialsName))
nCategory = len(categoryName)
 
data = read_mff_data.read_mff_data(filePath, 'epoch', 1, hdr['nTrials'], hdr)    
baseline = (nSamplesPre * 1000 / srate)  
msSamples = np.arange(0, nSamples,1) * 1000/srate  - baseline 
xlimMin = msSamples[0]
xlimMax = msSamples[-1]


H = np.identity(nC) - np.ones((nC, nC))/nC 
s = np.zeros(data.shape)  
if nTrials > 1:
    for i in range(nTrials):
        s[:,:,i] = np.dot(H, data[:,:,i])
    s2 = s * s 
    gfp = s2.sum(axis=0)      
    sGFP = np.sqrt(gfp)
else :
    s = np.dot(H, data)  
    s2 = s * s 
    gfp = s2.sum(axis=0)      
    sGFP = np.sqrt(gfp)
    
plt.plot(sGFP)
plt.show()
      
Data = data[:,:,0]
nStable = 200 

GEV = np.zeros((20,)) 
for k in range(2,20):
    VGT106_Kmeans = Kmeans.Kmeans(Data, k, nStable) 
    GEV[k] = VGT106_Kmeans['GEV']


fig, axes = plt.subplots(2,1)#, sharex=True)
axes[0].scatter(range(2,20), GEV[2:]) #, ylabel= 'GEV') 
axes[0].set_ylabel('GEV')
axes[1].scatter(range(2,20), np.log(GEV[2:]/(1600-np.arange(2,20))))  
axes[1].set_xlabel('K: number of clusters')
axes[1].set_ylabel('log(CV)')
#plt.subplots_adjust(hspace=0)
plt.savefig('VGT106GEVKopt.png')
plt.show() 


kOpt = 5

VGT106_Kmeans = Kmeans.Kmeans(Data, kOpt, nStable) 
GEV = VGT106_Kmeans['GEV'] 
GEVs = VGT106_Kmeans['GEVs']
C_UTmaps = VGT106_Kmeans['C_UTmaps']
KmeanID = VGT106_Kmeans['KmeanID']
Tmaps = VGT106_Kmeans['Tmaps']

GEV 
GFP = np.std(Data, axis= 0) 

colorcache = {u'blue': (0.0, 0.0, 1.0),
'green': (0.0, 0.5019607843137255, 0.0), 
'purple': (0.5019607843137255, 0.0, 0.5019607843137255),
 'r': (1.0, 0.0, 0.0),   
'firebrick': (0.6980392156862745, 0.13333333333333333, 0.13333333333333333), 
'cyan': (0.0, 1.0, 1.0),
'yellow': (1.0, 1.0, 0.0),
u'c': (0.0, 0.75, 0.75), 
u'k': (0.0, 0.0, 0.0),
'y': (0.75, 0.75, 0),
'darkgoldenrod': (0.7215686274509804, 0.5254901960784314, 0.043137254901960784),
'magenta': (1.0, 0.0, 1.0), 
u'0.75': (0.75, 0.75, 0.75), 
'0.8': (0.8, 0.8, 0.8), 
u'w': (1.0, 1.0, 1.0),
u'#bcbcbc': (0.7372549019607844, 0.7372549019607844, 0.7372549019607844), 
u'white': (1.0, 1.0, 1.0), 
u'#ffed6f': (1.0, 0.9294117647058824, 0.43529411764705883), 
u'#467821': (0.27450980392156865, 0.47058823529411764, 0.12941176470588237), 
u'#eeeeee': (0.9333333333333333, 0.9333333333333333, 0.9333333333333333), 
u'#F0E442': (0.9411764705882353, 0.8941176470588236, 0.25882352941176473), 
u'0.50': (0.5, 0.5, 0.5), 
u'#E24A33': (0.8862745098039215, 0.2901960784313726, 0.2), 
u'#f0f0f0': (0.9411764705882353, 0.9411764705882353, 0.9411764705882353), 
u'0.40': (0.4, 0.4, 0.4), 
'#afeeee': (0.6862745098039216, 0.9333333333333333, 0.9333333333333333),  
'0.5': (0.5, 0.5, 0.5), 
u'#fc4f30': (0.9882352941176471, 0.30980392156862746, 0.18823529411764706),  
u'0.00': (0.0, 0.0, 0.0), 
u'#bfbbd9': (0.7490196078431373, 0.7333333333333333, 0.8509803921568627), 
u'#ccebc4': (0.8, 0.9215686274509803, 0.7686274509803922), 
u'#A60628': (0.6509803921568628, 0.023529411764705882, 0.1568627450980392), 
u'#988ED5': (0.596078431372549, 0.5568627450980392, 0.8352941176470589), 
u'#777777': (0.4666666666666667, 0.4666666666666667, 0.4666666666666667), 
u'#EEEEEE': (0.9333333333333333, 0.9333333333333333, 0.9333333333333333), 
u'#fdb462': (0.9921568627450981, 0.7058823529411765, 0.3843137254901961), 
 u'#FFB5B8': (1.0, 0.7098039215686275, 0.7215686274509804), 
 u'#30a2da': (0.18823529411764706, 0.6352941176470588, 0.8549019607843137), 
 u'#555555': (0.3333333333333333, 0.3333333333333333, 0.3333333333333333), 
 u'#7A68A6': (0.47843137254901963, 0.40784313725490196, 0.6509803921568628), 
 u'#8b8b8b': (0.5450980392156862, 0.5450980392156862, 0.5450980392156862), 
 u'gray': (0.5019607843137255, 0.5019607843137255, 0.5019607843137255), 
 u'#8dd3c7': (0.5529411764705883, 0.8274509803921568, 0.7803921568627451), 
 u'#bc82bd': (0.7372549019607844, 0.5098039215686274, 0.7411764705882353), 
 u'#CC79A7': (0.8, 0.4745098039215686, 0.6549019607843137),  
 u'#E5E5E5': (0.8980392156862745, 0.8980392156862745, 0.8980392156862745), 
 u'0.70': (0.7, 0.7, 0.7), 
 u'#009E73': (0.0, 0.6196078431372549, 0.45098039215686275), 
 u'#FBC15E': (0.984313725490196, 0.7568627450980392, 0.3686274509803922), 
 u'#feffb3': (0.996078431372549, 1.0, 0.7019607843137254), 
 u'#56B4E9': (0.33725490196078434, 0.7058823529411765, 0.9137254901960784), 
 u'#e5ae38': (0.8980392156862745, 0.6823529411764706, 0.2196078431372549), 
 u'#348ABD': (0.20392156862745098, 0.5411764705882353, 0.7411764705882353), 
 u'#cbcbcb': (0.796078431372549, 0.796078431372549, 0.796078431372549),  
 u'#D55E00': (0.8352941176470589, 0.3686274509803922, 0.0), 
 u'#81b1d2': (0.5058823529411764, 0.6941176470588235, 0.8235294117647058),
 u'#8EBA42': (0.5568627450980392, 0.7294117647058823, 0.25882352941176473), 
 u'#0072B2': (0.0, 0.4470588235294118, 0.6980392156862745), 
 u'#6d904f': (0.42745098039215684, 0.5647058823529412, 0.30980392156862746),   
 '#00FFCC': (0.0, 1.0, 0.8), 
 u'#fa8174': (0.9803921568627451, 0.5058823529411764, 0.4549019607843137), 
 u'#b3de69': (0.7019607843137254, 0.8705882352941177, 0.4117647058823529)}
colorkeys = list(colorcache.keys())

colorkeys = ('b', 'g', 'r', 'c', 'm', 'y', 'k', 'gray','firebrick','darkgoldenrod','purple','#afeeee','#8EBA42','#7A68A6','#56B4E9','#D55E00')

KmeanIDcolor = [None] * 1600
for i in range(1600):
    KmeanIDcolor[i] = colorkeys[KmeanID[i]]


fig,axes = plt.subplots(2,1)#, sharex=True)
axes[0].plot(msSamples, Data.T)  
axes[0].set_ylabel(' ')
axes[1].vlines(msSamples, [0], GFP, color=KmeanIDcolor)  
axes[1].set_xlabel('ms')
axes[1].set_ylabel('GFP and 5 clusters')
plt.savefig('VGT106kmeans.png')
plt.show()



from pandas import Series, DataFrame 
import pandas as pd 

gsn257 = pd.read_table('/Users/jesong1126/Python3/gs3/nscolor_hgsn/GSN257ToPy.sfp') 
frame = DataFrame(gsn257, columns=['ChLabel', 'X3', 'Y3', 'Z3','X2','Y2'])
x2 = frame.values[:,4] 
y2 = frame.values[:,5] 


sss = pd.read_table('/Users/jesong1126/Python3/gs3/nscolor_hgsn/Seismic.clr', header=None) 
sssframe = DataFrame(sss)
sssframe
aaa = list(open('/Users/jesong1126/Python3/gs3/nscolor_hgsn/Seismic.clr'))

sesmic = sssframe.values
#sesmicr = sesmicframe.values[:,0] 
#sesmicg = sesmicframe.values[:,1] 
#sesmicb = sesmicframe.values[:,2] 
 
Tmaps0 = Tmaps 
Tmaps = Tmaps0 - np.mean(Tmaps0)
 
fig,axes = plt.subplots(1, kOpt, figsize=(15,3))#, sharex=True)
axes[0].scatter(x2, y2, c=Tmaps[:,0], s=20, cmap=plt.get_cmap('seismic'))
axes[0].set_alpha(0.75)
axes[0].set_title('Cluster 1')
axes[0].set_xticks([])
axes[0].set_yticks([])

axes[1].scatter(x2, y2, c=Tmaps[:,1],s=20, cmap=plt.get_cmap('seismic'))
axes[1].set_title('Cluster 2')
axes[1].set_xticks([])
axes[1].set_yticks([])

axes[2].scatter(x2, y2, c=Tmaps[:,2],s=20, cmap=plt.get_cmap('seismic'))
axes[2].set_title('Cluster 3')
axes[2].set_xticks([])
axes[2].set_yticks([])

axes[3].scatter(x2, y2, c=Tmaps[:,3],s=20, cmap=plt.get_cmap('seismic'))
axes[3].set_title('Cluster 4')
axes[3].set_xticks([])
axes[3].set_yticks([])

axes[4].scatter(x2, y2, c=Tmaps[:,4],s=20, cmap=plt.get_cmap('seismic'))
axes[4].set_title('Cluster 5') 
axes[4].set_xticks([])
axes[4].set_yticks([])
plt.savefig('VGT106kTMaps.png')
plt.show()


fig = plt.figure()
ax1 = plt.subplot2grid((3,5), (0,0), colspan=5)
ax1.plot(msSamples, Data.T)  
ax1.set_ylabel(' ')
ax1.set_xlim([msSamples[0], msSamples[-1]])
ax2 = plt.subplot2grid((3,5), (1,0), colspan=5)
ax2.vlines(msSamples, [0], GFP, color=KmeanIDcolor)  
ax2.set_xlim([msSamples[0], msSamples[-1]])
ax2.set_xlabel('ms')
ax2.set_ylabel('GFP ')
ax2.text(0, 4, '1', color='blue',fontsize=15)
ax2.text(500, 3.5, '2', color='green',fontsize=15)
ax2.text(750, 3.5, '3', color='red',fontsize=15)
ax2.text(1000, 3, '4', color='c',fontsize=15)
ax2.text(1300, 6, '5', color='m',fontsize=15)
ax3 = plt.subplot2grid((3,5), (2,0)) #, rowspan=2)
ax3.scatter(x2, y2, c=Tmaps[:,0], s=20, cmap=plt.get_cmap('seismic'))
ax3.set_alpha(0.75)
ax3.set_xlabel('Cluster 1')
ax3.spines['top'].set_color('blue')
ax3.spines['bottom'].set_color('blue')
ax3.spines['left'].set_color('blue')
ax3.spines['right'].set_color('blue')
ax3.set_xticks([])
ax3.set_yticks([])
ax4 = plt.subplot2grid((3,5), (2,1))
ax4.scatter(x2, y2, c=Tmaps[:,1], s=20, cmap=plt.get_cmap('seismic'))
ax4.set_alpha(0.75)
ax4.set_xlabel('Cluster 2')
ax4.spines['top'].set_color('green')
ax4.spines['bottom'].set_color('green')
ax4.spines['left'].set_color('green')
ax4.spines['right'].set_color('green')
ax4.set_xticks([])
ax4.set_yticks([])
ax5 = plt.subplot2grid((3,5), (2,2))
ax5.scatter(x2, y2, c=Tmaps[:,2], s=20, cmap=plt.get_cmap('seismic'))
ax5.set_alpha(0.75)
ax5.set_xlabel('Cluster 3')
ax5.spines['top'].set_color('red')
ax5.spines['bottom'].set_color('red')
ax5.spines['left'].set_color('red')
ax5.spines['right'].set_color('red')
ax5.set_xticks([])
ax5.set_yticks([])
ax6 = plt.subplot2grid((3,5), (2,3))
ax6.scatter(x2, y2, c=Tmaps[:,3], s=20, cmap=plt.get_cmap('seismic'))
ax6.set_alpha(0.75)
ax6.set_xlabel('Cluster 4')
ax6.spines['top'].set_color('c')
ax6.spines['bottom'].set_color('c')
ax6.spines['left'].set_color('c')
ax6.spines['right'].set_color('c')
ax6.set_xticks([])
ax6.set_yticks([])
ax7 = plt.subplot2grid((3,5), (2,4))
ax7.scatter(x2, y2, c=Tmaps[:,4], s=20, cmap=plt.get_cmap('seismic'))
ax7.set_alpha(0.75)
ax7.set_xlabel('Cluster 5')
ax7.spines['top'].set_color('m')
ax7.spines['bottom'].set_color('m')
ax7.spines['left'].set_color('m')
ax7.spines['right'].set_color('m')
ax7.set_xticks([])
ax7.set_yticks([])
plt.savefig('VGT106kButterGfpTMaps.png')
plt.show()

#fig,axes = plt.subplots(3,5)#, sharex=True)
#axes[0].plot(msSamples, Data.T)  
#axes[0].set_ylabel(' ')
#axes[1].vlines(msSamples, [0], GFP, color=KmeanIDcolor)  
#axes[1].set_xlabel('ms')
#axes[1].set_ylabel('GFP and 5 clusters')
#
#plt.savefig('VGT106kmeansAll.png')
#plt.show()


##  
#ax1 = fig.add_subplot(2,1,1)
#ax1.plot(msSamples, Data.T)  
#ax2 = fig.add_subplot(2,1,2)
#ax2.vlines(msSamples, [0], GFP, color=KmeanIDcolor)  
#plt.show()

## 
#plt.plot(msSamples, Data.T)  
#plt.show()


#fig = plt.figure()
#ax1 = fig.add_subplot(2,1,1)
#ax1.vlines(msSamples, [0], GFP, color=KmeanIDcolor)  
#ax2 = fig.add_subplot(2,1,2)
#ax2.scatter(msSamples, KmeanID)
#ax2.ylim([-1, k+1])
#plt.show()
#

#plt.vlines(msSamples, [0], GFP, color=KmeanIDcolor) #,cmap=plt.get_cmap('hsv')) 
#plt.xlim([msSamples[0], msSamples[-1]])
#plt.show()
#
#to_rgba(x, alpha=None, bytes=False)¶

N = 150
r = 2 * np.random.rand(N)
theta = 2 * np.pi * np.random.rand(N)
area = 200 * r**2 * np.random.rand(N)
colors = theta


#KmeanIDcolor = [None] * 1600
#for i in range(1600):
#    KmeanIDcolor[i] = colors[KmeanID[i]]


ax = plt.subplot(111, polar=True)
c = plt.scatter(theta, r, c=colors, s=area, cmap=plt.get_cmap('hsv'))
c.set_alpha(0.75)
plt.show()


    

#import matplotlib.colors 
#
#matplotlib.colors.rgb2hex([.2, .3, .1])
#matplotlib.colors.rgb_to_hsv([.2, .3, .1])
#matplotlib.colors.hsv_to_rgb([0.25, 0.66666667, 0.3]) 
#matplotlib.colors.rgb_to_hsv(colors)
#
#mycmap = plt.get_cmap('hsv') 
#
#imshow(mycmap)
#
#
#
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
#ax.imshow(gradient, aspect='auto', cmap=plt.get_cmap('hsv'))
#


    #    fprintf(sprintf('There are %d conditions:', nCondition));
    #for i=1:nCondition
    #    fprintf(sprintf(' %s ', eventsName{i}));
    #end
    #
    #fprintf(sprintf('\nSampling Rate is %d samples/second. \n', srate));
    #fprintf(sprintf('Total number of samples is %d and total Time period is %d ms. \n', nSamples, nSamples/srate * 1000));
    #fprintf(sprintf('Number of samples of pre-stimulus is %d and pre-stimulus time period is %d ms.\n', nSamplesPre, nSamplesPre/srate * 1000));
    #    
    
#print('nC:',nC, 'nSamples', nSamples, 'nTrials', nTrials) 
#        
#config.s = s 
#print("EEG data is loaded.")

#filePath = '/Users/jesong1126/Python3/gs3/data/SEP_108_0691_blc_ave_aref.mff' 
#filePath = '/Users/jesong1126/Python3/gs3/data/SEP_107_0046_fil_seg_bcr.mff'

#
#
#srate = 1000
#sss = getEpochInfos.getEpochInfos(filePath, srate) 
#hdr = read_mff_header.read_mff_header(filePath)
#nC = hdr['nChans']
#nSamples = hdr['nSamples']
#nSamplesPre = hdr['nSamplesPre']
#nTrials = hdr['nTrials']
#srate = hdr['Fs']
#summaryInfo = hdr['orig'] 
#epochLabels = summaryInfo['epochLabels']  
#
#s = read_mff_data.read_mff_data(filePath, 'epoch', 1, hdr['nTrials'], hdr) 
#



#hdr.keys() 
#['chantype', 'nSamples', 'chanunit', 'nSamplesPre', 'nTrials', 'orig', 'nChans', 'Fs', 'label'] 
 
 
 
 
 
 
 