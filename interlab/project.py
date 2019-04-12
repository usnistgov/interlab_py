from utilities import ExperimentGroup,InterlabArray,g_formatter,idfunc

import math
import numpy as np
import scipy as sp
import scipy.spatial 
import scipy.stats
import matplotlib
import matplotlib.pyplot as plt

from sklearn import preprocessing as skpp

import plot_utils as plu


#We use the tqdm notebook function and its widgets in order to run a lot of things, but we don't actually need it
try:
    import tqdm
    tqfunc = tqdm.tqdm_notebook
except ImportError:
    #if tqdm isn't available, just use the idfunc defined below
    tqfunc=idfunc

class Project(object):
    """The top-level project class for the interlaboratory comparison module
    
    :key Sample_names: List of sample names, used as keys for the dictionaries of data set names and processed and raw data. Each key in this list will correspond to a :py:class:`ExperimentGroup` object
    :key Data_set_names: Dictionary of data sets (labs) with data for each sample
    :key data: Dictionary of data to be used for the interlab analysis
    :key rawdata: Dictionary of unprocessed data, if different from data
    
    :key distance_metrics: List of distance metrics. Each metric in this list will be used to create a :py:class:`DistanceMetric` object within each :py:class:`ExperimentGroup` object
    
    :key x_data_list: The list of x data in the data array. For 2D data, this is not used
    :key range_to_use: Used to screen certain parts of the spectral data from consideration in the experimental comparison
    
    :key distribution_function: Which distribution will be assumed when assigning Z scores to each measurement of a sample. The default is sp.stats.lognorm
    :key outlier_dist: Which distribution will be assumed when detecting outliers. The default is the same as distribution_function

    
    """
    
    def __init__(self,
                 data=None,
                 rawdata=None,
                 distance_metrics=None,
                 Sample_names=None,
                 Data_set_names=None,
                 x_data_list=None,
                 range_to_use=None,
                 distribution_function=sp.stats.lognorm,
                 outlier_dist=None):
        
        # Feature names
        self._x_data=None
        if x_data_list is not None: self._x_data = np.array(x_data_list)
        
        # Names for the data sets (or laboratories) and samples ()
        self.Labels = Data_set_names
        self.Spectrum_names = Sample_names
        
        #Distance metrics for measuring the sample distances
        self.distance_metrics = distance_metrics
        
        #Range of data to use for analysis
        self._range_to_use = range_to_use
        
        self.zscores_array = dict
        
        self.pca = dict()
        self.projected_zscores = dict()
        self.outlier_mask = dict()
        self.zscore_pdf = dict()
        self.pdfs = dict()
        self._distribution_function = distribution_function
        self._outlier_distribution = outlier_dist
        if outlier_dist is  None:
            self._outlier_distribution = distribution_function
        #Create the experiment groups
        self.experiment_groups = []
        self._create_groups(xdata=self._x_data,data=data,rawdata=rawdata,data_names=Data_set_names)
        
        
        #Create the list of unique labs
        project_lablist = []
        
        for group in self.experiment_groups:
            project_lablist += list(group.data_names)
        
        self.lablist = list(np.unique(project_lablist))
        
        #Create the interlab arrays that will be used for outlier analysis        
        self.interlab_arrays = {}
        self._create_interlab_arrays()
        
        return
    
    def __str__(self):
        #allitems = self.experiment_groups + [self.interlab_arrays[metric['metric']] for metric in self.distance_metrics]
        exps_rep = ['Experiments:'] + [str(item) for item in self.experiment_groups]
        inls_rep = ['Interlab arrays:'] + [str(self.interlab_arrays[metric['metric']]) for metric in self.distance_metrics]
        
        str_rep = '\n  '.join(exps_rep) + '\n' + '\n  '.join(inls_rep)
        
        return str_rep
    
    @property
    def names(self):
        names = {}
        for item in self.experiment_groups:
            names[item.name] = item
        return names
    
    def __getitem__(self,x):
        if x is None:
            return self.experiment_groups
        elif type(x) is slice:
            return self.experiment_groups[x]
        elif type(x) is tuple:
            exp = x[0]
            metric = x[1]
            group = self._get_experiment(exp)
            if type(metric) is str:
                return group[metric]
            else:
                return [group[metric] for metric in x[1]]
        elif type(x) is list:
            groups = [self._get_experiment(exp) for exp in x]
            return groups
        else:
            return self._get_experiment(x)
        
        return
    
    def _get_experiment(self,exp):
        if type(exp) is int:
            return self.experiment_groups[exp]
        else:
            try:
                return self.names[exp]
            except KeyError:
                print('No experiment with name ' + exp)
    
    def _create_groups(self,xdata=None,data=None,rawdata=None,data_names=None,):
        
        for key in self.Spectrum_names:
            
            group_keyw = dict(name=key,
                              xdata=xdata,
                              data=data[key],
                              rawdata=rawdata[key],
                              data_names=data_names[key],
                              distance_metrics=self.distance_metrics,
                              distribution=self._distribution_function,
                              #mahalanobis_dict=self.mahalanobis_dict
                             )
            if data[key].shape[0] > 2:
                self.experiment_groups += [ExperimentGroup(**group_keyw)]
            else:
                print('Group {} has fewer than 3 measurements, not creating'.format(key))
        
        return
    
    def _create_interlab_arrays(self):
        
        for metric in self.distance_metrics:
            metric_name = metric['metric']
            self.interlab_arrays[metric_name] = InterlabArray(name=metric_name)
        return
    
    def set_distances(self):
        """Calculates the interspectral distances for each experiment group and metric
        """

        for metric in tqfunc(self.distance_metrics,desc='Distances'):
            metric_name = metric['metric']
            for group in tqfunc(self.experiment_groups,desc=metric_name):
                group.distance(metric_name)
    
    def fit_zscores(self):
        """Fits the sample-level zscores for each experiment group and metric.
        """
        for metric in tqfunc(self.distance_metrics,desc='Sample Zscores'):
            metric_name = metric['metric']
            for group in tqfunc(self.experiment_groups,desc=metric_name):
                group.fit_zscores(metric_name)
    
    def find_outliers(self,**kwargs):
        """Finds the sample outliers for each experiment group and metric.
        """
        for metric in tqfunc(self.distance_metrics,desc='Sample Outliers'):
            metric_name = metric['metric']
            for group in tqfunc(self.experiment_groups,desc=metric_name):
                group.find_outliers(metric_name,**kwargs)
    
    def fit_lab_zscores(self):
        """Fits the lab-level zscores for each metric.
        """
        for metric in tqfunc(self.distance_metrics,desc='Lab Zscores'):
            metric_name = metric['metric']
            self.interlab_arrays[metric_name].fit_zscores()
    
    def find_lab_outliers(self,**kwargs):
        """Finds the lab outliers for each metric.
        """
        for metric in tqfunc(self.distance_metrics,desc='Lab outliers'):
            metric_name = metric['metric']
            self.interlab_arrays[metric_name].find_outliers(**kwargs)
    
    def extract_matrices(self,**kwargs):
        """Runs :py:func:`.extract_experimental_matrix` for each experimental matrix
        """
        for metric in tqfunc(self.distance_metrics,desc='Extract matrices'):
            metric_name = metric['metric']
            self.extract_experimental_matrix(metric=metric_name,**kwargs)
    
    def plot_distance_fig(self,plot_range=None,
                          cmap='Greys',linecolor='k',
                          distance_metrics=None,plot_data=True,wspace=0,
                          ylabel_buffer=0.8,rightlabel_buffer=0.6,xlabel_buffer=0.5,**kwargs):
        """For each sample, generates the following plots:
          * A plot of the spectra generated for that sample by each laboratory
          * For each metric, a heat map plot of the interspectral distance matrix
        
        :key plot_range: An iterable of integers specifying which sample labels to plot
        :key cmap: The color map that will be used for the distance heat maps
        :key linecolor: The line color that will be used for the spectral data
        :key distance_metrics: A list of the distance metrics for which heat maps will be plotted. If None, plot heat maps for all metrics in this project
        :key plot_data: Boolean that tells whether the raw spectral data will be plotted
        :key wspace: Horizontal spacing between the heat maps
        :key ylabel_buffer: Space allocated for the y axis label (in inches)
        :key rightlabel_buffer: Space allocated for the colorbar label (in inches)
        :key xlabel_buffer: 
        :returns: distance_fig, the distance measure figure matplotlib object.
        """
        
        #If no distance metrics are specified, plot them all
        if distance_metrics is None: 
            distance_metrics = [metric['metric'] for metric in self.distance_metrics]
        
        sets_to_plot = self[plot_range]#self._scan_plot_range(plot_range)
        
        
        #Number of distance plots, number of sets to plot, and number of columns in the plot
        numdist = len(distance_metrics)
        numspecs = len(sets_to_plot)
        numcols = numdist

        # Height of each row, in inches
        height = 1.5

        #Width of the spectrum column and distance-metric column, in terms of row height
        spectrum_width = 2
        distance_width = 1.5

        widths = []
        total_width = 0

        if plot_data:
            numcols += 1
            widths += [spectrum_width]
            total_width += spectrum_width
        for i in range(numdist):
            widths += [distance_width] # Add the appropriate number of columns to the width list
            total_width += distance_width

        fs = (total_width * height,numspecs * height) #total figure size, inches
        gs_key = dict(width_ratios=widths,hspace=0.075,wspace=wspace)
        
        #Create the figure we're going to plot in
        sfig,saxes = plu.make_grid_plot(numspecs,numcols,figsize=fs,gridspec_kw=gs_key,
                                        sharey='row',#sharex='col',
                                        right_ylabel=r'Pairwise Distance, $D_{ij,k}$',
                                        xlabel='Data set ID',ylabel='Data set ID',
                                        ylabel_buffer=ylabel_buffer,
                                        rightlabel_buffer=rightlabel_buffer,
                                        xlabel_buffer=xlabel_buffer,
                                        add_buffer=True)
        show_titles = True
        
        #Create the distance measure plots
        for group,ax_row in zip(sets_to_plot,saxes):
            group.distance_measure_plot(ax_row,plot_data=plot_data,
                                        data_kwargs=dict(linecolor=linecolor),
                                        distance_kwargs=dict(cmap=cmap,**kwargs),
                                        distance_metrics=distance_metrics,
                                        show_titles=show_titles)          
            show_titles = False
        
        #Change the heat map ranges so that they're all the same
        for metric_num,metric_name in enumerate(distance_metrics):
            
            #metric_name = metric['metric']
            
            #Find the least minimum and greatest maximum for all of the groups for this metric
            plot_min = []
            plot_max = []
            for group in sets_to_plot:
                pmn,pmx = group[metric_name].minmax()
                plot_min += [pmn]
                plot_max += [pmx]
            distance_min = np.array(plot_min).min()
            distance_max = np.array(plot_max).max()
            
            #Expand the limits of the images to be the same as the least min and greatest max
            for row in saxes:
                plot_locator = metric_num + int(plot_data) #need to increment this by 1 if the data is being plotted with the heat maps
                im = row[plot_locator].get_images()[0]
                im.set_clim(vmin=distance_min,vmax=distance_max)
        
        return sfig
    
    
    def process_mahalanobis(self,sets_to_extract=None,threshold=1e-1):
        """Calculates the pooled covariance for use by Mahalanobis distance
        """
        
        #Extract the groups we want to use for the Mahalanobis distance
        groups_to_use = self[sets_to_extract]#self._scan_plot_range(sets_to_extract)
        
        #Default values need defined before running
        covariance = None
        denom = 0

        for group in groups_to_use:
            #Mean-center the data
            data_centered = (group.data.T - group.data.mean(axis=1))
            
            #Covariance matrix for this group
            cov_this = np.cov(data_centered)
            if covariance is None:
                covariance = np.zeros_like(cov_this)
            
            #This appears a lot because statistics
            std_factor = group.data.shape[1] - 1
            
            #Start adding together the weighted covariances
            covariance += std_factor * cov_this
            denom += std_factor
        #Final step of pooling
        covariance = covariance / denom
        
        #Pseudo-inverse
        VI = np.linalg.pinv(covariance,rcond=threshold)
        
        #self.mahalanobis_dict = dict(VI=VI,)
        
        for group in groups_to_use:
            group['Mahalanobis'].mahalanobis_dict = dict(VI=VI,)
        
    
    
    def _extract_matrix(self,sets=None,metric=None,screen_outliers=True):
        '''Turns the dict-of-zscores into a zscores array with experiments on one axis and labs on the other. Replaces missing lab/experiment combinations with np.nan'''
        if metric is None: raise ValueError('Must specify the metric to use')
        
        sets_to_extract = self[sets]#self._scan_plot_range(sets)
        #if sets is None: sets = [group.name for group in sets_to_extract]
            
        #Create the array of zscores for experiments
        shape = (len(sets),len(self.lablist))
        zscores_array = np.empty(shape)
        zscores_array.fill(np.nan)
        
        for group in sets_to_extract:
            #Extract z-scores and outliers            
            zscores = group[metric].zscores
            if screen_outliers:
                outliers = ~group[metric].outlier_mask #outlier_mask MASKS the outliers, we want the opposite
            else:
                #If we're not screening outliers, the outlier mask should just be all zeros, i.e. nothing is an outlier
                outliers = np.zeros_like(group[metric].outlier_mask,dtype=np.bool)
            data_names = group.data_names
            
            #key_location and lab_location tell us where to find the inlier and outlier labs
            key_location = sets.index(group.name)
            for outlier,zscore,data_name in zip(outliers,zscores,data_names):
                lab_location = self.lablist.index(data_name)
                
                #If this element in zscores_array is nan, then we haven't filled it yet
                if np.isnan(zscores_array[key_location,lab_location]):
                    if not outlier:
                        zscores_array[key_location,lab_location] = zscore
                #Otherwise, we've filled it, so take the average
                else:
                    if not outlier:
                        zscores_array[key_location,lab_location] = (zscores_array[key_location,lab_location] + zscore)/2
        return zscores_array
    
    def extract_experimental_matrix(self,
                                    sets=None,
                                    metric=None,
                                    screen_outliers=True,
                                    imputation_axis=0):
        """Extracts the zscore data from the dict-of-vectors format and casts it as a 2D array.
        
        The dictionary of sample-level z-scores is recast as an array, with one dimension corresponding to sample names and the other corresponding to laboratory.
        
        :key sets: Sets to extract for the interlab comparison
        :key metric: Distance metric that will be used
        :key screen_outliers: Whether to remove outlier measurements before imputing missing values
        :key imputation_axis: Axis along which to impute missing values
        """
        
        if sets is None: sets = [group.name for group in self.experiment_groups]
        
        zscores_array = self._extract_matrix(sets=sets,metric=metric,screen_outliers=screen_outliers)
        
        #Find labs that have had all their spectra removed as outlier spectra
        no_exp_for_lab = np.all(np.isnan(zscores_array),axis=0)
        
        #Impute the missing Z scores by replacing all nan values with the median
        imp =  skpp.Imputer(strategy='median',axis=imputation_axis) 
        zscores_imputed = imp.fit_transform(zscores_array.T[~no_exp_for_lab]).T
        
        #Load the imputed score array into the interlab array object
        self.interlab_arrays[metric].data_array = zscores_imputed
        self.interlab_arrays[metric].lablist = list(np.array(self.lablist)[~no_exp_for_lab])
        self.interlab_arrays[metric].datasets = sets
        
        #PCA transformation and projected statistial distance calculation
        self.interlab_arrays[metric].fit_transform()
        return
    
    def plot_zscores_heatmap(self,
                             sets=None,
                             metric=None,
                             screen_outliers=True,
                             cmap='Greys',**subplot_kwargs):
        """Extracts the zscore data from the dict-of-vectors format and casts it as a 2D array, plotting it with missing and outlier data called out.
        
        The dictionary of sample-level z-scores is recast as an array, with one dimension corresponding to sample names and the other corresponding to laboratory. The array is then plotted as a heat map.
        :key sets: Sets to extract for the interlab comparison
        :key metric: Distance metric that will be used
        :key screen_outliers: Whether to remove outlier measurements before imputing missing values
        :key imputation_axis: Axis along which to impute missing values
        """
        if sets is None: sets = [group.name for group in self.experiment_groups]
        
        #Extract the experiment matrix without screening outliers
        zscores_array = self._extract_matrix(sets=sets,metric=metric,screen_outliers=False)
        
        #Find the smallest finite value in the matrix (NaNs are missing) and set the vmin for the heatmap to that value. Then fill the NaNs with zeros
        vmin = zscores_array[np.isfinite(zscores_array)].min()
        zscores_array[np.isnan(zscores_array)] = 0
        
        #By default, we are not flagging outliers
        vmax = None
        extend = 'min'
        
        if screen_outliers:
            extend = 'both'
            #Find the outlier limit, which is the 95th (94.95) percentile, and set vmax to that value
            std = self._distribution_function(1)
            vmax = std.ppf(0.9495)
        
        #Create the figure and axis objects
        fig,ax = plt.subplots(**subplot_kwargs)
        
        #Get the colormap and then set the over- and under-colors
        cmaps = plt.cm.get_cmap(cmap)
        cmaps.set_over('r') #Red for outliers
        cmaps.set_under('0.6') #Grey for missing data
        
        #Plot the heat map and color bar
        im = ax.imshow(zscores_array,cmap=cmaps,origin='lower',aspect='auto',vmax=vmax,vmin=vmin)
        cb = fig.colorbar(im,ax=ax,extend=extend)
        
        #Axis labels
        ax.set_ylabel('Experiment',size=15)
        ax.set_xlabel('Data Set ID',size=15)
        cb.set_label('Z-score',size=15)
        
        #X tick labels for facilities
        a = ax.set_xticks(range(len(self.lablist)))
        a = ax.set_xticklabels(self.lablist,rotation='vertical')
        
        #Y tick labels for experiment groups
        a = ax.set_yticks(range(len(sets)))
        a = ax.set_yticklabels(sets)

    
    def _standard_plot(self,numplots,numcols,gs_kw=dict(wspace=0.1,hspace=0.4),**kwargs):
        '''Creates a standard plot that we use a lot'''
        
        numrows = int(math.floor((numplots + numcols - 1)/ numcols))
        height = 1
        width = 5
        
        fig, ax_array = plu.make_grid_plot(numrows,numcols,
                                           plotsize=(width,height),
                                           gridspec_kw=gs_kw,
                                           **kwargs)
        axes = ax_array.flatten('F')
        return fig,axes
    
    def _scan_plot_range(self,plot_range=None):
        #If no experimental groups are specified, show them all
        if plot_range is None: 
            sets_to_plot = self.experiment_groups
        #Otherwise, find the groups by name
        else:
            sets_to_plot = []
            for name in plot_range:
                for group in self.experiment_groups:
                    if group.name == name: sets_to_plot += [group]
        return sets_to_plot
    
    def plot_histograms(self,metric,plot_range=None,numcols=2,xlabel_buffer=0.8,rotation=None) :
        """Plots a histogram of the average interspectral distance for each sample, along with the corresponding fit
        
        :param metrics: The metric for which the distances will be plotted
        :key plot_range: An iterable of integers specifying which sample labels to plot
        :key numcols: The number of columns in the distance plot
        :key xlabel_buffer: Space allocated for x-axis labels (in inches)
        :key rotation: Specifies the orientation of the z-score labels for individual labs
        :returns pdffig: The distances and scores plot as a matplotlib figure object
        """
        sets_to_plot = self[plot_range]#self._scan_plot_range(plot_range)
        
        numplots = len(sets_to_plot)
        
        xlabel = r"Average Diameter Distance, $\^{D}_{i,k}$"
        ylabel = "Probability Density"
        
        gs_kw = dict(wspace=0.1,hspace=0.4)
        
        pdffig,pdfax = self._standard_plot(numplots,numcols,gs_kw=gs_kw,
                                           xlabel_buffer=xlabel_buffer,
                                           ylabel=ylabel,xlabel=xlabel,
                                           add_buffer=True,
                                           sharex=True,sharey=True,
                                          )
        plot_min = []
        plot_max = []
        for group in sets_to_plot:
            plot_min += [group[metric].distance_min]
            plot_max += [group[metric].distance_max]

        distance_min = np.array(plot_min).min()
        distance_max = np.array(plot_max).max()
        
        #Create the histogram bins and corresponding PDF evaluation points
        num_bins = 20
        max_bin = distance_max
        min_bin = distance_min
        binsize = (max_bin - min_bin)/float(num_bins)
        lognorm_xvals = np.linspace(min_bin,max_bin,num_bins*10)
        
        #Ticks
        pdf_locator = plt.MaxNLocator(nbins=4,prune='both')
        
        #Plot each histogram
        for ax,group in zip(pdfax,sets_to_plot):
            group.histogram(metric,ax,
                            num_bins=num_bins,max_bin=max_bin,min_bin=min_bin)
            
            #Count number of non-outlier samples
            num_samples = (np.count_nonzero(group[metric].outlier_mask),len(group.data_names))
            
            #Annotate the histograms
            plot_text = '\n'.join((group.name,'N = {} / {}'.format(*num_samples)))
            
            ax.text(0.98,0.9,plot_text,ha='right',va='top',transform=ax.transAxes)
            ax.yaxis.set_major_locator(pdf_locator)
        
        return pdffig
    
    def plot_zscore_distances(self,metric,plot_range=None,numcols=2,
                              xlabel_buffer=1,ylabel_buffer=0.6,
                              rotation=None) :
        """Plots a bar chart of the average interspectral distance for each sample, annotated with the generalized Z score for each sample
        
        :param metric: The metric for which the distances will be plotted
        :key plot_range: An iterable of integers specifying which sample labels to plot
        :key numcols: The number of columns in the distance plot
        :key xlabel_buffer: Space allocated for x-axis labels (in inches)
        :key ylabel_buffer: Space allocated for y-axis labels (in inches)
        :key rotation: Specifies the orientation of the z-score labels for individual labs
        :returns zscorefig: The distances and scores plot as a matplotlib figure object
        """
        sets_to_plot = self[plot_range]#self._scan_plot_range(plot_range)
    
        numplots = len(sets_to_plot)
        
        ylabel = r"Average Diameter Distance, $\^{D}_{i,k}$"
        xlabel = 'Data set ID'
        if numplots == 2:
            ylabel = ylabel.replace('r D','r\nD')
        if numplots < 2:
            ylabel = ylabel.replace(' ','\n')
        
        gs_kw = dict(wspace=0.1,hspace=0.4)
        
        #By default, share axes among all plots
        sharex = True
        olm0 = sets_to_plot[0][metric].outlier_mask.shape #Get the shape of the outlier mask of the first group
        #Check to make sure the outlier masks are the same length
        for group in sets_to_plot:
            #Check if the lablists of each experimental group are the same length
            mask_same = group[metric].outlier_mask.shape == olm0
            if not(mask_same):
                sharex = False #Can't share axes if different number of labs
        
        zscorefig,zscoreax = self._standard_plot(numplots,numcols,gs_kw=gs_kw,
                                                 xlabel_buffer=xlabel_buffer,
                                                 ylabel_buffer=ylabel_buffer,
                                                  ylabel=ylabel,xlabel=xlabel,
                                                  add_buffer=True,
                                                  sharex=sharex,sharey=True,)
        
        #Ticks
        pdf_locator = plt.MaxNLocator(nbins=4,prune='both')
        
        for ax,group in zip(zscoreax,sets_to_plot):
            group.plot_zscores(metric,ax,rotation=rotation)
            ax.text(0.98,1.02,group.name,ha='right',va='bottom',transform=ax.transAxes)
            ax.yaxis.set_major_locator(pdf_locator)
            if not sharex: ax.set_xticklabels([])

        return zscorefig
    
    
    def plot_projected_zscores(self,distance_metrics=None,xlabel_buffer=1,rotation=None):
        """Plots the projected statistical distances annotated with the corresponding laboratory-level Z scores.
        
        :key distance_metrics: A list of the distance metrics for which statistical distances will be plotted. If None, plot statistical distances for all metrics in this project
        :key xlabel_buffer: Space allocated for x-axis labels (in inches)
        :key rotation: Specifies the orientation of the z-score labels for individual labs
        :returns: zscorefig, the projected statistical distances plot as a matplotlib figure object
        """       
        #Distance metrics that will be used
        if distance_metrics is None: 
            interlabs = self.interlab_arrays
            #[print(interlab) for interlab in interlabs]
            distance_metrics = [interlab for interlab in interlabs]
        interlabs = []
        for metric in distance_metrics:
            interlabs += [self.interlab_arrays[metric]]
        
        #Number of plots, size of figure
        numplots = len(interlabs)
        numcols = 1
        numrows = int(math.floor((numplots + numcols - 1)/ numcols))
        width = 5
        height = 1
        
        #Xlabel and ylabel
        ylabel = r"Projected Statistical Distance, $||T_L||$"
        xlabel = 'Data set ID'
        
        #Collapese ylabel if there aren't many rows
        if numplots == 2:
            ylabel = ylabel.replace('l D','l\nD')
        if numplots < 2:
            ylabel = ylabel.replace(' ','\n')
        
        gs_kw = dict(wspace=0.1,hspace=0.4)
        
        #By default, share axes among all plots
        sharex = True
        #Check to make sure the outlier masks are the same length
        for interlab in interlabs:
            #Check if the lablists are the same length
            mask_same = interlab.outlier_mask.shape == interlabs[0].outlier_mask.shape
            if not(mask_same):
                sharex = False #Can't share axes if different number of labs
        
        zscorefig,zscoreax_array = plu.make_grid_plot(numrows,numcols,plotsize=(width,height),
                                                  xlabel_buffer=xlabel_buffer,
                                                  ylabel=ylabel,xlabel=xlabel,
                                                  add_buffer=True,
                                                  sharex=sharex,sharey=True,gridspec_kw=gs_kw)
        zscoreax = zscoreax_array.flatten('F')
        
        for ax,interlab,metric in zip(zscoreax,interlabs,distance_metrics):
            interlab.plot_zscores(ax,rotation=rotation)
            
            ax.set_title(metric)
            axis_locator = plt.MaxNLocator(nbins=5,prune='both')
            ax.yaxis.set_major_locator(axis_locator)
            if not sharex: ax.set_xticklabels([])
        
        return zscorefig
    
    def plot_zscore_loadings(self,distance_metrics=None,xlabel_buffer=1):
        """Plots the principal component loadings for the statistical distances
        
        :key distance_metrics: A list of the distance metrics for which loadings will be plotted. If None, plot loadings for all metrics in this project
        :key xlabel_buffer: Space allocated for x-axis labels (in inches)
        :returns: loadfig, the projected statistical loadings plot as a matplotlib figure object
        """
        
        #Distance metrics that will be used
        if distance_metrics is None: 
            interlabs = self.interlab_arrays
            #[print(interlab) for interlab in interlabs]
            distance_metrics = [interlab for interlab in interlabs]
        interlabs = []
        for metric in distance_metrics:
            interlabs += [self.interlab_arrays[metric]]
        
        #Number of plots, size of figure
        numplots = len(interlabs)
        numcols = 1
        numrows = int(math.floor((numplots + numcols - 1)/ numcols))
        width = 5
        height = 1
        
        ylabel = u"Component Value"
        xlabel = 'Cluster ID'
        #if numplots < 2:
        #    ylabel = ylabel.replace(' ','\n')
        
        gs_kw = dict(wspace=0.1,hspace=0.4)
        
        loadfig,loadax_array = plu.make_grid_plot(numrows,numcols,plotsize=(width,height),
                                                  xlabel_buffer=xlabel_buffer,
                                                  ylabel=ylabel,xlabel=xlabel,
                                                  add_buffer=True,
                                                  sharex=True,sharey=True,gridspec_kw=gs_kw)
        loadax = loadax_array.flatten('F')
        
        for ax,interlab,metric in zip(loadax,interlabs,distance_metrics):
            interlab.plot_components(ax)
            
            ax.set_title(metric)
            axis_locator = plt.MaxNLocator(nbins=5,prune='both')
            ax.yaxis.set_major_locator(axis_locator)
        return loadfig
    
    def plot_zscore_outliers(self,metric,y_component=1,text=True):
        """Plots the principal component scores for each lab along with the final distribution used to calculate the outliers
        
        :param metric: The metric used to calculate the interspectral distances
        :key y_component: Which principal component to use on the Y axis, if not the first
        :key text: Whether to label the plot with the name of the 
        :returns: zscore_outliers_fig, the Z score outlier plot as a matplotlib figure object
        """
        zscore_outliers_fig,ax = plt.subplots(1,1,figsize = (10,5))
        self.interlab_arrays[metric].plot_outliers(ax,y_component=y_component)
        
        if text:
            zscore_outliers_fig.text(0.85,0.8,metric,ha='right',va='top',size=15,transform=ax.transAxes)
        return zscore_outliers_fig