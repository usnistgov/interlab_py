import scipy as sp
import scipy.stats
import math
import numpy as np

import matplotlib.pyplot as plt
import matplotlib.gridspec as gs

from sklearn import decomposition as skdecomp


def idfunc(*args,**kwargs):
    '''Returns exactly the function's positional arguments, ignoring any keyword arguments'''
    if len(args) == 1:
        return args[0] #We only passed one argument to the function. We want to return arg, not (arg)
    return args

def g_formatter(x):
    a,b = '{:.2e}'.format(x).split('e')
    b = int(b)
    a = float(a)
    if not(b > 2) or (b < -1): 
        mult = 10**b
        a = a * mult
        prec = min(max(2-b,0),2)
        return (r'${:.'+str(prec)+'f}$').format(a)
    return r'${}\times10^{{{}}}$'.format(a,b)

def calculate_minmax(data):
    """Determines the maximum and minimum (excluding zero) of a data set. Rounds to one significant figure.
    """
    data_max = data.max()
    data_min = data[data > 0].min()
    
    #Find the power of 10 
    log_max = math.log(data_max,10)
    log_min = math.log(data_min,10)
    
    decimals_max = int(-1*math.floor(log_max))+1
    decimals_min = int(-1*math.floor(log_min))
    
    significance_max = 10 ** decimals_max
    significance_min = 10 ** decimals_min
    
    distance_max = math.ceil(data_max * significance_max)/significance_max
    distance_min = math.floor(data_min * significance_min)/significance_min
    
    return distance_min,distance_max

class ExperimentGroup(object):
    """Container for spectral data and methods for analysis
    
    This object contains a group of spectral measurements that are related, usually by being of the same physical object and having the same structure.
    
    :key name: The name of the group of measurements.
    :key data: Array of data to be used for the interlab analysis
    :key rawdata: Array of unprocessed data, if different from data
    
    :key data_names: List of data sets (labs) with data for each sample
    
    :key distance_metrics: List of distance metrics. Each metric in this list will be used to create a :py:class:`DistanceMetric` object.
    :key distribution_function: Which distribution will be assumed when assigning Z scores to each measurement of a sample.
    """
    def __init__(self,
                 name=None,
                 xdata=None,
                 rawdata=None,data=None,data_names=None,
                 distance_metrics=None,distribution=None):
        
        self.name = name
        self.xdata = xdata
        self.data = data
        self.data_names = data_names
        self.distance_metrics = []
        self.distribution = distribution
        
        for metric in distance_metrics:
            self.distance_metrics += [
                DistanceMetric(distribution=distribution,**metric)
            ]
        
        return
    
    @property
    def names(self):
        names = {}
        for item in self.distance_metrics:
            names[item.name] = item
        return names
    
    def __str__(self):
        
        str_rep = self.name
        
        return str_rep
    
    def __getitem__(self,x):
        if type(x) is int:
            return self.distance_metrics[x]
        else:
            try:
                return self.names[x]
            except KeyError:
                print('No metric with name ' + x)
        
        return
    
    def distance(self,metric):
        """Pass-through to :py:class:`.DistanceMetric`
        """
        self[metric].distance(self.data)
        return
    
    def fit_zscores(self,metric):
        """Pass-through to :py:class:`.DistanceMetric`
        """
        self[metric].fit_zscores()
        return
    
    def find_outliers(self,metric,**kwargs):
        """Pass-through to :py:class:`.DistanceMetric`
        """
        self[metric].find_outliers(**kwargs)
        return
    
    def histogram(self,metric,*args,**kwargs):
        """Pass-through to :py:class:`.DistanceMetric`
        """
        self[metric].histogram(*args,**kwargs)
        return
    
    def plot_zscores(self,metric,*args,**kwargs):
        """Pass-through to :py:class:`.DistanceMetric`
        """
        ax = self[metric].plot_zscores(*args,**kwargs)
        
        ax.set_xlim(-0.5,len(self.data_names)-0.5)
        ax.set_xticks(np.arange(len(self.data_names)))
        ax.set_xticklabels(self.data_names,rotation='vertical')
        return
    
    def minmax(self,metric):
        return self[metric].population.minmax()
    
    def plot_data(self,ax,linecolor='k',**kwargs):
        """Line plot of the spectral data in this group.
        """
        numsets = self.data.shape[0]
        if self.xdata is not None:
            pl = ax.plot(self.xdata,
                         self.data.T/self.data.max(axis=1)+np.array(range(numsets))-0.5/float(1),
                         linecolor)
            if self.xdata[0] > self.xdata[1]:
                ax.invert_xaxis()
        else:
            pl = ax.plot(self.data.T/self.data.max(axis=1)+np.array(range(numsets))-0.5/float(1),linecolor)
        ax.set_ylim([-0.5,numsets-0.5])
        ax.text(0.1,0.95,self.name,ha='left', va='center',color='k',size=9,
                 bbox=dict(facecolor='w',alpha=1,lw=0,pad=0),transform=ax.transAxes)
        
        
        return
    
    def distance_measure_plot(self,ax_row,plot_data=False,distance_metrics=None,
                              data_kwargs={},distance_kwargs={},show_titles=False):
        """Heatmaps of interspectral distances and (optionally) line plots of data
        """
        
        if distance_metrics is None:
            metrics_to_plot = self.distance_metrics
        else:
            metrics_to_plot = []
            for metric in self.distance_metrics:
                for name in distance_metrics:
                    if metric.name == name: 
                        metrics_to_plot += [metric]
        
        ax_dist = ax_row
        if plot_data:
            ax_dist = ax_row[1:]
            
            self.plot_data(ax_row[0],**data_kwargs)
        for metric,ax in zip(metrics_to_plot,ax_dist):
            
            im = metric.distance_plot(ax,show_titles=show_titles,
                                      **distance_kwargs)
            
            cb = plt.colorbar(im,ax=ax,fraction=0.2)
            cb.locator = plt.MaxNLocator(nbins=4,prune='lower',integer=True,min_n_ticks=3)
            cb.update_ticks()
            
            ax.set_yticks(range(len(self.data_names)))
            ax.set_yticklabels([])
            ax.set_xticks([])
        ax_row[0].set_yticklabels(self.data_names,size=7)
        
    
class DistanceMetric(object):
    """A distance metric and the interspectral distances associated with it
    
    :param metric: The name of the metric function
    :param function: The metric function. Can be a binary function f(x,y) or a string identifying a metric recognized by :py:func:`scipy.spatial.distance.pdist`
    :key distribution_function: Which distribution will be assumed when assigning Z scores to each measurement of a sample.
    :key mahalanobis_dict:
    """
    
    def __init__(self,metric,function,distribution=None,mahalanobis_dict=None):
        self.name = metric
        self.function = function
        self.distribution = distribution
        
        self.distance_matrix = None
        self.distance_max = None
        self.distance_min = None
        self.population = None
        self.mahalanobis_dict = mahalanobis_dict
        
        return
    
    @property
    def zscores(self):
        return self.population.zscores
    
    @property
    def values(self):
        return self.population.values
    
    @property
    def outlier_mask(self):
        return self.population.outlier_mask
    
    def distance(self,data):
        """Calculates the pairwise interspectral distances and the average interspectral distance, then loads the average distances into the :py:class:`Population`.
        
        :param data: The data for which distances will be calculated
        """
        #VI=None It used to work that passing VI=None for distance metrics that didn't require a 
        kwargs = {}
        #significant_components=data.shape[0]

        #Empty containers
        distance_dict = dict()
        distance_square_dict = dict()

        if self.function == 'mahalanobis':
            if not(self.mahalanobis_dict):
                raise ValueError('Must define the inverse covariance matrix if metric is mahalanobis')
            #VI=self.mahalanobis_dict['VI']
            kwargs = self.mahalanobis_dict
            #significant_components=mahalanobis_dict['significant_components']
        
        #distance_condensed = sp.spatial.distance.pdist(data,self.function,VI=VI)
        distance_condensed = sp.spatial.distance.pdist(data,self.function,**kwargs)
        self.distance_matrix = sp.spatial.distance.squareform(distance_condensed)
        
        
        self.population = Population(name=self.name,
                                     distribution=self.distribution,
                                     values=self.distance_matrix.sum(axis=1)/self.distance_matrix.shape[0]
                                    )
        self.distance_min,self.distance_max = self.population.minmax()
        return
    
    def distance_plot(self,ax,vmin=None,vmax=None,cmap='Greys',origin='lower',
                      aspect='equal',show_titles=False,**kwargs):
        """Plots a heatmap of the pairwise interspectral distance.
        
        :param ax: The axis to make the plot
        :key vmin: Passed to imshow
        :key vmax: Passed to imshow
        :key cmap: Passed to imshow
        :key origin: Passed to imshow
        :key aspect: Passed to imshow
        :key show_titles:
        """
        if vmin is None: vmin,vmax = self.minmax()
        
        if show_titles:
            ax.set_title(self.name.replace(' ','\n'),size=10)
        
        im = ax.imshow(self.distance_matrix,cmap=cmap,origin=origin,
                       vmin=vmin,vmax=vmax,aspect=aspect,**kwargs
                      )
        return im
        
    
    def fit_zscores(self,**kwargs):
        """Pass-through to :py:class:`.Population`
        """
        self.population.fit_zscores(**kwargs)
        return
    
    def find_outliers(self,recursive=False,**kwargs):
        """Pass-through to :py:class:`.Population`
        """
        self.population.find_outliers(recursive=recursive,**kwargs)
        if recursive:
            self.distance_min,self.distance_max = calculate_minmax(self.values[self.outlier_mask])
        return
    
    def histogram(self,*args,**kwargs):
        """Pass-through to :py:class:`.Population`
        """
        self.population.histogram(*args,**kwargs)
        return
    
    def plot_zscores(self,*args,**kwargs):
        """Pass-through to :py:class:`.Population`
        """
        ax = self.population.plot_zscores(*args,**kwargs)
        return ax
    
    def minmax(self):
        return calculate_minmax(self.distance_matrix)
    
    

class Population(object):
    """Object containing some values, a distribution function fit to those values, and corresponding scores
    
    :param name: The name assigned to this Population
    :key distribution: Which distribution will be assumed when assigning Z scores to laboratories.
    :key values: The values to which the distribution will be fit
    """
    
    def __init__(self,name,distribution=None,values=None):
        self.name = name
        self.values = values
        self.distribution = distribution
        self.params = None
        self.zscores = None
        self.outlier_mask = None
        return
    
    def fit_zscores(self,data=None,mask=None):
        """Fit the distribution to the data using a provided mask and calculate the scores of the data
        
        :key data: The data that will be fit to the distribution. Normally, this is the values of the Population
        :key mask: If present, masks some elements of data
        """
        
        if data is None: data = self.values
        if mask is None: mask = np.ones_like(data,dtype=np.bool)
        
        std = self.distribution(1)        
        self.params = self.distribution.fit(data[mask],floc=0)
        
        rv = self.distribution(*self.params)
        
        CDF = rv.cdf(data)
        self.zscores = std.ppf(CDF)
        
        #Explicitly calculate Z if the distribution is lognorm, because this will give more reliable results for large Z
        if type(self.distribution) is type(sp.stats.lognorm):
            mean = self.params[2]
            stdev = self.params[0]
            self.zscores = np.exp( (np.log(data) - np.log(mean)) / stdev)
        
        #self.outlier_mask = self.zscores < zlimit
        
        return
    
    def find_outliers(self,recursive=False,support_fraction=0.6,final_screen=False):
        """Finds the outliers of the distribution
        
        :key recursive: If true, finds outliers and refit until all outliers have been removed
        :key support_fraction: Fraction of values that are guaranteed to be retained
        """
        
        max_percentile = 0.9495 #95% confidence limint, respecting support fraction (respecting 3 sig figs rounding)
        max_ignore_support_fraction = 0.9895 #99% confidence limit, ignoring support fraction (3 sig fig rounding)
        
        std = self.distribution(1)
        zlimit = std.ppf(max_percentile)
        zlimit_ignore = std.ppf(max_ignore_support_fraction)
        
        self.fit_zscores(self.values)
        max_num_outliers = int((1-support_fraction)*len(self.values))
        
        self.outlier_mask = np.ones_like(self.values,dtype=np.bool)
        
#         if not recursive:
#             for i in range(max_num_outliers):
#                 self.outlier_mask = self.zscores < zlimit #If not recursive, just get all outliers and leave
#             return
        
        # Note that, because we always use a support fraction, we always want to remove outliers one-by-one
        for i in range(max_num_outliers):
            zmax = self.zscores[self.outlier_mask].max() #Get the biggest z-score still considered in the distribution
            if zmax < zlimit: break #If the biggest z-score is within the zlimit, break from the loop
            self.outlier_mask = self.zscores < zmax
            #Recalculate z-scores only if recursive
            if recursive: self.fit_zscores(self.values,self.outlier_mask) 
        if final_screen:
            self.outlier_mask = self.zscores < zlimit_ignore
            if recursive: self.fit_zscores(self.values,self.outlier_mask) 
        return
    
    def minmax(self):
        return calculate_minmax(self.values)
    
    def histogram(self,ax,num_bins=20,max_bin=None,min_bin=None,xvals=None):
        """Plots a histogram of the populations `values` with the corresponding distribution function
        """
        if max_bin is None: min_bin,max_bin=self.minmax()
        if xvals is None: xvals = np.linspace(min_bin,max_bin,num_bins*10)
        
        binsize = (max_bin - min_bin)/float(num_bins)
        
        hist_norm,bins = np.histogram(self.values,range=(min_bin,max_bin), bins=num_bins,density=True)
        dist_vals = self.distribution.pdf(xvals,*self.params)
        confidence_limit = self.distribution.ppf(0.9495,*self.params)

        ax.plot(xvals,dist_vals,'k')
        ax.bar(bins[:-1]+0.5*binsize,hist_norm,binsize,color='r')
        ax.axvline(x=confidence_limit,color='k',ls=':')
        
        return
    
    def plot_zscores(self,ax,rotation=None):
        """Plots a bar chart of the populations `values` and annotates it with the corresponding zscores
        """
        width = 0.4
        
        #Get the number of data sets and number of spectra
        numsets = len(self.zscores)
        
        bar_colors = ['w'] * numsets
        bar_edge = ['k'] * numsets
        box_trans = [0.7] * numsets
        for int_num,inlier in enumerate(self.outlier_mask):
            #Color and transparency if the point is an outlier
            if not(inlier):
                bar_colors[int_num] = 'r'
                bar_edge[int_num] = 'r'
                box_trans[int_num] = 0
        
        #bar chart of projected distances
        ax.bar(np.arange(numsets),self.values,
               2*width,alpha=1,color=bar_colors,edgecolor=bar_edge)
        
        
        (ymin,ymax) = ax.get_ylim()
        sample_label_pos = 0.05*ymax + 0.95*ymin
        
        for int_num,(zscore,value,inlier) in enumerate(zip(self.zscores,self.values,self.outlier_mask)):
            if(inlier): 
                label_pos = sample_label_pos + value
            else:
                label_pos = sample_label_pos
            ax.text(int_num,label_pos,
                      g_formatter(zscore),
                      ha='center', va='bottom',
                      color='k',size=9,rotation=rotation,
                      bbox=dict(facecolor='w',alpha=box_trans[int_num],lw=0)
                     )
        return ax
        

class InterlabArray(object):
    """Object for executing interlaboratory comparison
    
    :key name: Name for the interlab object
    :key data_array: Z-score array for the measurements in this interlaboratory study
    :key lablist: List of laboratories that have data in this interlab
    :key datasets: List of datasets that are in this interlab
    :key distribution_function: Which distribution will be assumed when assigning Z scores to laboratories. The default is sp.stats.lognorm
    """
    
    def __init__(self,
                 name=None,
                 data_array=None,
                 lablist=None,
                 datasets=None,
                 distribution_function=sp.stats.lognorm,):
        
        #Data array that will be fit
        self.data_array = data_array
        
        #Distribution function to fit the zscores to
        self.distribution_function = distribution_function
        
        #Metric name
        self.name = name
        
        #PCA object
        self.pca = skdecomp.PCA()
        self.pca_scores = None
        self.population = Population('Projected statistical distance',
                                     distribution=self.distribution_function)
        
        return
    
    def __str__(self):
        return self.name
    
    @property
    def zscores(self):
        return self.population.zscores
    
    @property
    def values(self):
        return self.population.values
    
    @property
    def outlier_mask(self):
        return self.population.outlier_mask
    
    def fit_transform(self):
        """Fits the sklearn.pca object and then removes the mean-centering. Calculates projected statistical distance and stores that in :py:class:`Population`
        """
        #Fit PCA
        self.pca.fit(self.data_array.T)
        
        #Find where the zero point of our z-score space lies in PCA space
        effective_zero = self.pca.transform(np.zeros_like(self.pca.mean_.reshape(1,-1)))
        
        #Transform the scores into the translated PCA space
        self.pca_scores = self.pca.transform(self.data_array.T) - effective_zero
        
        #Get the active components (minimum number of components that explains at least 95% of variance)
        #self.active_components = self.pca.explained_variance_ratio_.cumsum() < 0.95
        self.active_components = self.pca.explained_variance_ratio_ > 0.01
        # ^^^ this definition actually has one fewer component than what we want so we have to do some manipulation
        
        #print(self.active_components)
        
        self.active_components[1] = True #Obviously there should be at least one active component
        #self.active_components = np.roll(self.active_components,1) #Move the active component mask one element. Combined with the statement above, this ensures that the second component will be active, and we want two active components.
        #self.active_components[0] = True #Obviously the first component should be active
        
        #print(self.active_components)
        #Store the projected statistical distance
        num_active_components = np.count_nonzero(self.active_components)
        self.population.values = np.linalg.norm(self.pca_scores[:,0:num_active_components],axis=1)

        return
    
    #def plot_zscore_heatmap(self):
        
    
    def fit_zscores(self,**kwargs):
        """Pass-through for :py:class:`.Population.fit_zscores`
        """
        #Fit the projected statistical distance
        self.population.fit_zscores(**kwargs)
        return
    
    def find_outliers(self,**kwargs):
        """Pass-through for :py:class:`.Population.find_outliers`
        """
        self.population.find_outliers(**kwargs)
        return
    
    def plot_components(self,ax):
        """Plots the principal component loadings for the statistical distances
        
        :param ax: A list of the distance metrics for which loadings will be plotted. If None, plot loadings for all metrics in this project
        :key xlabel_buffer: Space allocated for x-axis labels (in inches)
        """
        PCAsign = np.sign(self.pca.components_[0,:].mean())
        
        #Get the number of active components in the PCA
        num_active_components = np.count_nonzero(self.active_components)
        #print(self.name,num_active_components,self.active_components)

        #Calculate the positions where the bars will go
        total_width = 0.4 # total width available for all components
        width = total_width/num_active_components #width allocated to each component
        positions = np.arange(-total_width,total_width,2*width) #list of offsets for all components

        #Create the spectrum of colors that will be used for the plot
        bar_properties = np.arange(1,0,-1.0/num_active_components)

        #Get the number of data sets and number of spectra
        numsets = len(self.lablist)
        numspecs = len(self.datasets)


        #Plot the principal components
        for component_num in range(num_active_components):
            #intensity = bar_properties[component_num]

            #Intensity corresponds to the explained variance ratio
            intensity = 1 - self.pca.explained_variance_ratio_[component_num]

            #Edges should be darker to provide contrast for lighter bars
            edge_intensity = (1 - self.pca.explained_variance_ratio_[component_num])/1.3

            bar_colors = [[intensity] * 3] * len(self.lablist)
            edge_colors = [[edge_intensity] * 3] * len(self.lablist)

            position = positions[component_num]
            im = ax.bar(np.arange(len(self.datasets))+position,
                            self.pca.components_[component_num,:]*PCAsign,
                            2*width,
                            alpha=1,color=bar_colors,edgecolor=edge_colors)
        
        #X tick labels
        ax.set_xlim(-1,numspecs + 2) # We'll want 2 labels for the explained variance ratios
        ax.set_xticks(range(numspecs))
        ax.set_xticklabels(self.datasets,rotation='vertical')
        
        #Draw dividing lines between components
        for spec_num,spectrum in enumerate(self.datasets):
            ax.axvline(x=spec_num+0.5,color='k',ls=':')
        
        #Print the explained variance ratios
        #(ymin,ymax) = ax.get_ylim()
        #center = (ymax+ymin)/2.0
        #interval = (ymin-ymax)/(num_active_components + 1)
        #positions = np.arange(ymax+interval,ymin,interval)
        interval = -1/(num_active_components + 1)
        positions = np.arange(1+interval,0,interval)

        for component_number in range(num_active_components):
            component_label_position = positions[component_number]
            component_horizontal_position = len(self.datasets)/(len(self.datasets) + 1)
            componentscore = "%4.1f" % (self.pca.explained_variance_ratio_[component_number] * 100)
            componentscore_string = componentscore + " %"
            ax.text(component_horizontal_position,component_label_position,componentscore_string,
                        ha='center', va='center',transform=ax.transAxes)
        return
    
    def plot_zscores(self,*args,**kwargs):
        """Plots the projected statistical distances annotated with the corresponding laboratory-level Z scores. Pass-through to :py:class:Population.plot_zscores.
        """ 
        ax = self.population.plot_zscores(*args,**kwargs)
        
        ax.set_xlim(-0.5,len(self.lablist)-0.5)
        ax.set_xticks(np.arange(len(self.lablist)))
        ax.set_xticklabels(self.lablist,rotation='vertical')
        return
    
    
    def plot_outliers(self,ax,y_component=1):
        """Plots the principal component scores for each lab along with the final distribution used to calculate the outliers
        
        :param ax: The axis on which the plot will be made
        :key y_component: Which principal component to use on the Y axis, if not the first
        """
        zscores_projected = self.pca_scores
        z_mask = self.outlier_mask
        
        #Create the scatterplot, blue for outliers and red for inliers
        scatterplot = ax.scatter(self.pca_scores[z_mask,0],
                                     self.pca_scores[z_mask,y_component],
                                         c='r',edgecolors='none',s=60,zorder=3)
        scatterplot = ax.scatter(self.pca_scores[np.invert(z_mask),0],
                                     self.pca_scores[np.invert(z_mask),y_component],
                                         c='b',edgecolors='none',s=60,zorder=4)

        #Get the x and y limits of the plot
        (xmin,xmax) = ax.get_xlim()
        (ymin,ymax) = ax.get_ylim()

        #Make the grid 
        numgrid = 60

        xstep = (xmax-xmin)/numgrid
        ystep = (ymax-ymin)/numgrid

        xdata = np.arange(xmin,xmax,xstep)
        ydata = np.arange(ymin,ymax,ystep)

        xx,yy = np.meshgrid(xdata,ydata)
        dist = np.sqrt(xx ** 2 + yy ** 2)

        rv = self.population.distribution(*self.population.params)
        std = self.population.distribution(1)
        pdf_for_plot = std.ppf(rv.cdf(dist))

        levelsc = [1, 2, 3, 4, 5]
        levelsf = [0.25, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5]

        levelsc = list(np.arange(1,std.ppf(0.95))) + [std.ppf(0.95)]
        levelsf = [0.25] + list(np.arange(0.5,std.ppf(0.95),0.5)) + [std.ppf(0.95)]

        cformats =  ['{:1.0f}'] * (len(levelsc) - 1)
        cformats += ['{:3.1f}']

        cformats_dict = {x:fmt_str.format(x) for x,fmt_str in zip(levelsc,cformats) }
        ccolors = ['w'] * int(len(levelsc)/2)
        ccolors += ['k'] * int( (len(levelsc) + 1)/2)

        ctr1 = ax.contourf(xx,yy,pdf_for_plot,extend='min',cmap='Greys_r',levels=levelsf,zorder=1)
        ctr2 = ax.contour(xx,yy,pdf_for_plot,extend='neither',colors=ccolors,levels=levelsc,zorder=2)
        contourlabels = plt.clabel(ctr2,fmt=cformats_dict,colors=ccolors)
        component1score = "%5.2f" % (self.pca.explained_variance_ratio_[0] * 100)
        component2score = "%5.2f" % (self.pca.explained_variance_ratio_[y_component] * 100)
        xaxislabel = 'Component 1: ' + component1score + " %"
        yaxislabel = 'Component ' + str(y_component + 1) +': ' + component2score + " %"

        ax.set_xlabel(xaxislabel,size=20)
        ax.set_ylabel(yaxislabel,size=20)
        
        for j,labeltext in enumerate(self.lablist):
            if not z_mask[j]:
                ax.text(self.pca_scores[j,0]+1.0,self.pca_scores[j,y_component],
                        labeltext,
                        ha='left',va='center',
                        size='9',
                        bbox = dict(facecolor='w'))
        
