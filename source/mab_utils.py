import struct
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.patches as mplp
import mpl_toolkits.axes_grid1.inset_locator as mplax
import os
import pandas as pd
#import tqdm

def idfunc(*args,**kwargs):
    '''Returns exactly the function's positional arguments, ignoring any keyword arguments'''
    if len(args) == 1:
        return args[0] #We only passed one argument to the function. We want to return arg, not (arg)
    return args

#We use the tqdm notebook function and its widgets in order to run a lot of things, but we don't actually need it
try:
    import tqdm
    tqfunc = tqdm.tqdm_notebook
except ImportError:
    #if tqdm isn't available, just use the idfunc defined above
    tqfunc=idfunc
    
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

# File reading and metadata processing utilities

def nmr_read(file):
    '''Reads a binary file from the NMR spectral reader
    
    :parameter file: The file to be read
    :type file: str
    :returns: data, an array of the spectral values stored in file
    '''
    data = [] #Blank list
    with open(file,'rb') as f:
        f.read(2048) # Header information, skip
        while True:
            binary_data = f.read(4) # the NMR data is stored as 4-byte real numbers. Read the 4 bytes
            if not binary_data: break #If reach EOF, break out of the loop
            data += [struct.unpack('f',binary_data)] #interpret
    data = np.array(data)
    return data  

def key_from_name(name):
    '''Takes a tuple of strings and concatenates them together. Used for turning Pandas groupby keys into strings for interlab
    '''
    #print(name)
    key = str(name)#'_'.join(name)
    if type(name)==type(()):
        key =  '-'.join(name)
    return key

def exptype(df):  
    '''Reads the Pandas DataFrame containing the experimental metadata to find the NMR field frequency. Spectra are split into low-field (<700 MHz), high-field (800 MHz and greater), and mid-field (everything else).
    '''
    types = []
    for s in df.TITLE:
        base_type = s[s.find('-')+5:s.find('-')+8] #This is the portion of the title containing the NMR frequency
        gross_type = '600'
        if base_type.startswith('8') or base_type.startswith('9'): gross_type = '800'
        if base_type.startswith('7'): gross_type = '700'
        types += [gross_type]
    return types

# def exptype(df):
#     '''Reads the Pandas DataFrame containing the experimental metadata to find the experimental type, which is the [UN][freq] portion of the title.
    
#     :parameter df: The Pandas DataFrame containing the experimental metadata
#     :returns: types, a list of the experimental types extracted from the df.TITLE column
#     '''
#     types = []
#     for s in df.TITLE:
#         types += [s[s.find('-')+3:s.find('-')+8]]
#     return types
def expcode(df):
    '''Reads the Pandas DataFrame containing the experimental metadata to find the experimental code, which is the [DE][123][ABCD] portion of the title
    
    :parameter df: The Pandas DataFrame containing the experimental metadata
    :returns: types, a list of the experimental codes extracted from the df.TITLE column
    '''
    types = []
    for s in df.TITLE:
        types += [s[:s.find('-')+2]]
    return types

def getfile(df,prefix='./study/interp/'):
    '''Extracts the filename from the Pandas DataFrame containing the experimental metadata'''
    #prefix = './study/interp/'
    filelist = os.listdir(prefix)
    files = []
    for i in df.INDEX:
        #print(i)
        begin_str = '{:0>3}'.format(i)
        for filename in filelist:
            if filename.startswith(begin_str): 
                file = prefix+filename
                break
        #file = prefix+filelist[int(i) - 1]
        files += [file]
    return files

def read_all_metadata(
    metadata_file=None,
    raw_data_path=None,
    keep_nmr_corrections_data=False,
    excode=expcode,
    extype=exptype,
):
    '''Reads the NMR spectral metadata from metadata_file and builds the Pandas dataframe containing all spectral metadata for later analysis. Raw data is found in raw_data_path and is loaded into metadata_table as a pointer.
    '''
    #Read the metadata table into a Pandas dataframe
    df = pd.read_table(metadata_file,
                       comment='#',
                      )
    #Empty list
    parsed_data = []

    #Split the rows of the metadata table by whitespace and build into a Numpy array
    for string_array in df.values[1:]:
        parsed_string = np.array(string_array[0].split())
        parsed_data += [parsed_string]
    parsed_data = np.array(parsed_data).T
    
    #Determine whether we are keeping just the lab ID metadata or if we need to keep the NMR adjustment data
    last_col = 5
    if keep_nmr_corrections_data:
        last_col = None
    column_names = df.columns.values[0].split()[1:last_col]
    column_names
    
    #Create the metadata table from the parsed version of the text metadata table
    metadata_dict = {}
    for column_name,column_data in zip(column_names,parsed_data): #Zip ensures that there will be no mismatched lengths
        metadata_dict[column_name] = column_data #Make a dict
    metadata_table = pd.DataFrame.from_dict(metadata_dict) #Build dataframe from dict

    #Extract the experiment type and experiment codes from the metadata table
    metadata_table = metadata_table.assign(
        #ExpType=mab_utils.exptype(metadata_table),
        ExpType=extype(metadata_table),
        ExpCode=excode(metadata_table),
        File=getfile(metadata_table,prefix=raw_data_path))
    return metadata_table
    
def _read_all_metadata(metadata_file=None,raw_data_path=None):
    '''Reads the NMR spectral metadata from metadata_file and builds the Pandas dataframe containing all spectral metadata for later analysis. Raw data is found in raw_data_path and is loaded into metadata_table as a pointer.
    '''
    #Read the metadata table into a Pandas dataframe
    df = pd.read_table(metadata_file,
                       comment='#',
                      )
    #Empty list
    parsed_data = []

    #Split the rows of the metadata table by whitespace and take the first 16 elements
    for string_array in df.values[1:]:
        parsed_string = np.array(string_array[0].split())
        parsed_data += [parsed_string[:16]]

    #Create the metadata table from the parsed version of the text metadata table
    metadata_table = pd.DataFrame(columns=df.columns.values[0].split()[1:16],data=parsed_data)

    #Extract the experiment type and experiment codes from the metadata table
    metadata_table = metadata_table.assign(
        #ExpType=mab_utils.exptype(metadata_table),
        ExpType=exptype(metadata_table),
        ExpCode=expcode(metadata_table),
        File=getfile(metadata_table,prefix=raw_data_path))
    return metadata_table

def read_all_data(metadata_table,indices_to_split,split_point,shape=None):
    '''Splits the spectra in the Pandas dataframe based on indices_to_split. Spectra are grouped based on indices_to_split[:split_point] and further subdivided by indices_to_split[split_point:]. 
    '''
    
    data_list = []
    data_dict = dict()
    titles_dict = dict()
    
    #Experiments will be sorted by group_indices; measurements with the same group_indices will be comparable
    #lower_indices are to ensure that data are sorted by a unique identifier
    group_indices = indices_to_split[:split_point]#['ExpCode','ExpType']
    lower_indices = indices_to_split[split_point:]
    for name,split_table in tqfunc(metadata_table.groupby(group_indices)):
        description = key_from_name(name)

        #Empty container for this set of measurements
        nmrdata = []
        titles = []
        for codename, explist in tqfunc(split_table.groupby(lower_indices),desc=description,):
            titles += [explist.DIR_NAME.values]
            #Have to iterate over the measurements with this experiment because some are duplicated
            for file in explist.File.values:
                data = nmr_read(file)
                data = data.reshape(shape)
                data = data.reshape(-1,1)
                nmrdata += [data]
        #Make list into array
        nmrdata = np.transpose(np.concatenate(nmrdata,axis=1))
        titles = np.concatenate(titles)
        #
        data_list += [nmrdata]
        data_dict[description] = nmrdata
        titles_dict[description] = titles
    return data_list,data_dict,titles_dict
    
#Data processing utilities

def process_all_data(metadata_table,full_data_dict,indices_to_split,split_point,thresh=0):
    
    data_dict = dict()
    
    for description,nmrdata in full_data_dict.items():
        
        nmrprocessed = []
        for rawdata in nmrdata:
            #Crush values below threshold to 0
            data=np.zeros_like(rawdata)
            data[np.abs(rawdata)>thresh] = rawdata[np.abs(rawdata)>thresh] 

            #Convert to colormap
            data,data_colors = mab_utils.colorize(data,interp_map,cm)

            #Check for zero and negative values (not actually necessary after mapping)
            data = inl.fix_spectrum(data) 

            data = data/data.sum() #Sum-normalize

            #np.dot spits out a 1-D vector, need to add another axis
            nmrprocessed += [data[:,np.newaxis]] 
        nmrprocessed = np.concatenate(nmrprocessed,axis=1).T
        data_dict[description] = nmrprocessed
    
    return data_dict

def save_outlier_data(metadata_table,outlier_project,indices_to_split,split_point,metric):
    #Experiments will be sorted by group_indices; measurements with the same group_indices will be comparable
    #lower_indices are to ensure that data are sorted by a unique identifier
    group_indices = indices_to_split[:split_point]#['ExpCode','ExpType']
    lower_indices = indices_to_split[split_point:]
    #Iterate over 
    for name,data in metadata_table.groupby(group_indices):

        #Correct the key
        description = key_from_name(name)

        #Find the corresponding Experimental Group
        for this_group in outlier_project.experiment_groups:
            if this_group.name == description:break
            #this_group = mab_project[description]
        #print(this_group.name)
        data_metric = this_group[metric]

        #Load the outliers into the metadata table and print the distances and corresponding z scores
        nmr_start = 0
        for codename, explist in data.groupby(lower_indices):
            nmr_end = nmr_start + len(explist)
            this_dist = data_metric.values[nmr_start:nmr_end]
            this_z = data_metric.zscores[nmr_start:nmr_end]
            this_out = ~data_metric.outlier_mask[nmr_start:nmr_end]
            metadata_table.loc[explist.index.values,'Zscore'] = this_z
            metadata_table.loc[explist.index.values,'Outlier'] = this_out
            nmr_start = nmr_end
    

# def make_colormap(color_positions,map_name):
#     #Keys for the color_dict that will be used by matplotlib
#     color_keys = ['red','green','blue','alpha']
#     color_dict = dict(red=None,green=None,blue=None,alpha=None)
    
#     #Iterate over position,color pairs
#     for position,color in color_positions:
#         #Match RGB(A) values in color with the appropriate key
#         #This will make sure that if the color is only three values, the alpha tuple does not get created
#         for key,color_val in zip(color_keys,color):
#             if color_dict[key]:
#                 color_dict[key] += ((position,color_val,color_val),)
#             else:
#                 color_dict[key] = ((position,color_val,color_val),)
#     #Make and return the color map
#     cmap = mcolors.LinearSegmentedColormap(map_name,color_dict)
#     return cmap
    
def colorize(data,interp_map,colormap,greyscale_conversion=[0.3,0.59,0.11,0]):
    '''Converts data to a colormap based on the interpolation map interp_map, and also converts the colormap to a greyscale map using greyscale_conversion.
    '''
    data = np.interp(data, *interp_map) #Use the linear map to map the data between 0 and 1
    data_colors = colormap(data) #Convert to colormap
    scaled_data = np.dot(data_colors,greyscale_conversion) #Convert to greyscale
    return scaled_data,data_colors

#Midpoint Normalize for non-symmetric color maps
class MidpointNormalize(mcolors.Normalize):
    '''This is the custom normalization example taken from https://matplotlib.org/tutorials/colors/colormapnorms.html 
    '''
    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
        self.midpoint = midpoint
        mcolors.Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        # I'm ignoring masked values and all kinds of edge cases to make a
        # simple example...
        x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
        return np.ma.masked_array(np.interp(value, x, y))

def make_diverging_colormap(map_name,
                            midpoint=0.5,crush_point=0.25,
                            min_color=(0,0,0.5,),
                            crush_lo_color=(0,0,1),
                            mid_color=(1,1,1),
                            crush_hi_color=(1,0,0),
                            max_color=(0.3,0,0)
                 ):
    '''Reproduces seismic colormap by default'''
    crush_lo = midpoint - crush_point
    crush_hi = midpoint + crush_point
    zero = 0
    one = 1
    
    color_positions = (
        (zero,     min_color),
        (crush_lo, crush_lo_color),
        (midpoint, mid_color),
        (crush_hi, crush_hi_color),
        (one,      max_color)
    )
    
    cmap = mcolors.LinearSegmentedColormap.from_list(map_name,list(color_positions))
    return cmap

seismic_with_black = make_diverging_colormap('seismic_with_black',mid_color=(0,0,0))
seismic_with_alpha = make_diverging_colormap('seismic_with_black',crush_point=0.05,
                                             mid_color=(1,1,1,1))
RdBu_with_black = make_diverging_colormap('RdBu_with_black',crush_point=0.05,
                                          min_color=(0,0,1),
                                          crush_lo_color=(0,0,0),
                                          crush_hi_color=(0,0,0),
                                          max_color=(1,0,0),
                                          mid_color=(0,0,0))
RdBu_with_white = make_diverging_colormap('RdBu_with_white',crush_point=0.05,
                                          min_color=(0,0,1),
                                          crush_lo_color=(1,1,1),
                                          crush_hi_color=(1,1,1),
                                          max_color=(1,0,0),
                                          mid_color=(1,1,1))
    
    
#Plot utilities

def nmr_plot(data,cmap=None,norm=None,color='purple',
             plot_corners=None,levels=None,ax=None,figsize=None):
    '''Generates a contour plot of the NMR spectral intensities
    
    :parameter data: The array of NMR spectral intensities
    :type data: 2d array
    :key cmap: The colormap that will be used for the contour map
    :key norm: The norm that will be used for the contour map if cmap is used
    :key color: The colors used for the contour plot. Only used if cmap is None
    :key plot_corners: The corners of the plot in data space, used to generate the axes of the contour plot
    :key levels: The levels of the contour plot
    :key ax: a matplotlib axis into which the contour plot will be entered. If None, creates a new figure
    :key figsize: If ax is None, defines the size of the new figure  
    
    '''
    if levels is None:
        raise Exception('Contour levels must be defined')
    
    if ax is None:
        fig,ax=plt.subplots(figsize=figsize)
    
    contour_dict = dict(extent=plot_corners,
                        aspect='auto',
                        origin='lower',
                        levels=levels,
                        )
    if cmap is not None:
        contour_dict['cmap'] = cmap
        contour_dict['norm'] = norm
    else:
        contour_dict['colors'] = color
    im = ax.contour(data,**contour_dict)

    return ax,im

def make_grid_plot(numrows,numcols,figsize=None,plotsize=None,
                   column_width=6,row_height=4,
                   label_buffers=None,
                   ylabel_buffer=0.75,xlabel_buffer=0.5,
                   xlabel=None,ylabel=None,
                   add_buffer=False,
                   **subplots_args):

    if plotsize is not None:
        column_width,row_height = plotsize
    
    if label_buffers is not None:
        xlabel_buffer,ylabel_buffer = label_buffers
    
    full_width = numcols*column_width
    full_height = numrows*row_height
    if add_buffer:
        full_width = numcols*column_width + ylabel_buffer
        full_height = numrows*row_height + xlabel_buffer
    
    bottom_buffer = xlabel_buffer/full_height
    left_buffer = ylabel_buffer/full_width

    ylabel_pos = 0.5*(1+bottom_buffer)
    xlabel_pos = 0.5*(1+left_buffer)
  
    fs = (full_width,full_height)
    if figsize is not None:
        fs = figsize
    fig,axes = plt.subplots(numrows,numcols,figsize=fs,squeeze=False,**subplots_args)
    fig.subplots_adjust(left=left_buffer,right=1,top=1,bottom=bottom_buffer)
    
    if ylabel:
        fig.text(0,ylabel_pos,ylabel,size=15,rotation='vertical',va='center',ha='center')
    if xlabel:
        fig.text(xlabel_pos,0.0,xlabel,ha="center",va="center",size=15)
        
    return fig,axes

def nmrbig(explist,ncols=2,
           cmap=None,color='purple',levels=[1,5,10,15,20],
           plot_corners=None,shape=None,window=True,**grid_args):
    
    if plot_corners is None: plot_corners = [2.002,-0.799,27.513,7.005]
    if shape is None: shape = (1024,1315)

    nrows = len(explist) // ncols
    if nrows < len(explist) / ncols:
        nrows += 1
    fig,axes = make_grid_plot(nrows,ncols,**grid_args)
    bbox = dict(facecolor='w')

    for (index,row),ax in zip(explist.iterrows(),axes.flatten()):
        
        file = row.File
        index = '{:0>3d}'.format(int(row.INDEX))
        etype = row.ExpType
        code = row.ExpCode
        zscore = row.Zscore
        outlier = row.Outlier
        laboutlier = row.LabOutlier
        
        data = nmr_read(file)
        data = data.reshape(shape)
        if window: data = data[256:957,453:1272,]
        #plot_corners = [2.002,-0.799,27.513,7.005]
        #print(data.max())
        nmr_plot(data,plot_corners=plot_corners,levels=levels,ax=ax,color=color,cmap=cmap)
        
        plot_label = '-'.join((index,code,etype))
        zscore_string = 'Z = {:4.3f}'.format(zscore)
        plot_label = ': '.join((plot_label,zscore_string))
        
        ax.text(0.95,0.05,plot_label,ha='right',va='bottom',transform=ax.transAxes,bbox=bbox)
        if outlier:
            ax.text(0.05,0.95,'This is an outlier',ha='left',va='top',transform=ax.transAxes,bbox=bbox)
        elif laboutlier:
            tx = 'Lab outlier, z = {:5.2f}'.format(row.LabZscore)
            ax.text(0.05,0.95,tx,ha='left',va='top',transform=ax.transAxes,bbox=bbox)
    return fig,axes

def plot_the_nmr(explist,ncols=2,
                 color=None,cmap=seismic_with_black,norm=None,
                 colorbar=False,levels=np.array([0,2,5,10,15,20]),
                 plot_corners=None,shape=None,window=True,**grid_kw):

    nrows = len(explist) // ncols
    if nrows < len(explist) / ncols:
        nrows += 1
    fig,axes = make_grid_plot(nrows,ncols,**grid_kw)
    bbox = dict(facecolor='w')
    
    if plot_corners is None: plot_corners = [2.002,-0.799,27.513,7.005]
    if shape is None: shape = (1024,1315)
        
    try:
        all_laboutliers = explist.LabOutlier
        laboutliers = True
    except AttributeError:
        laboutliers = False

    for (index,row),ax in zip(explist.iterrows(),axes.flatten()):
        
        file = row.File
        index = '{:0>3d}'.format(int(row.INDEX))
        etype = row.ExpType
        code = row.ExpCode
        outlier = row.Outlier
        if laboutliers: 
            laboutlier = row.LabOutlier
        else:
            laboutlier = False
        


        
        data = nmr_read(file)
        data = data.reshape(shape)
        if window: data = data[256:957,453:1272,]
        #print(data.max())
        
        plot_contours = levels[levels > data.min()]
        
        nmr_plot_dict = dict(plot_corners=plot_corners,
                             levels=plot_contours,
                             ax=ax,
                            )
        if norm is not None: nmr_plot_dict['norm'] = norm

        if color is None:
            plot_color = cmap
#             if outlier: 
#                 plot_color='Reds'
#             elif laboutlier:
#                 plot_color='Oranges'
            nmr_plot_dict['cmap'] = plot_color
        else:
            if len(color) > 1: 
                plot_color = color[levels > data.min()]
            else:
                plot_color = color
            nmr_plot_dict['color'] = plot_color
            
        _,im = nmr_plot(data,**nmr_plot_dict)
        
        labtext = '-'.join((index,code,etype))
        
        zscore_for_print = g_formatter(row.Zscore)
        zscorestring = 'z = {}'#'z = {:5.2f}'
        
        tx = zscorestring.format(zscore_for_print)#row.Zscore)
        
        if colorbar: 
            cb = fig.colorbar(im,ax=ax)
            cb.set_label('NMR field intensity (% of max)',size=15)
        
        ax.text(0.95,0.05,labtext,ha='right',va='bottom',transform=ax.transAxes,bbox=bbox)
        if outlier:
            tx = 'This is an outlier, '+tx
        #    ax.text(0.05,0.95,tx,ha='left',va='top',transform=ax.transAxes,bbox=bbox)
        elif laboutlier:
            tx = 'Lab outlier, '+tx
        ax.text(0.05,0.95,tx,ha='left',va='top',transform=ax.transAxes,bbox=bbox)
            
    return fig,axes

def better_mark_inset(big_ax,inset_ax,corner_1,corner_2,inset_spine_color
                     ):
        #mplax.mark_inset(ax,ax_sub,loc1=4,loc2=3,ls=':',color='w')
    
    #Get edges of inset axis in data space
    (inset_start_x,inset_end_x) = inset_ax.get_xlim()
    (inset_start_y,inset_end_y) = inset_ax.get_ylim()
    inset_width = inset_end_x - inset_start_x
    inset_height = inset_end_y - inset_start_y
    
    rect = mplp.Rectangle((inset_start_x,inset_start_y),
                          inset_width,inset_height,
                          edgecolor=inset_spine_color,facecolor='none',ls=':')
    big_ax.add_artist(rect)
    big_ax.figure.canvas.draw()
    corner_of_axis = big_ax.transData.inverted().transform(inset_ax.transAxes.transform(corner_1))
    corner_of_data = inset_ax.transData.inverted().transform(inset_ax.transAxes.transform(corner_1))
    big_ax.plot((corner_of_axis[0],corner_of_data[0]),(corner_of_axis[1],corner_of_data[1]),
            ls=':',color=inset_spine_color)
    
    corner_of_axis = big_ax.transData.inverted().transform(inset_ax.transAxes.transform(corner_2))
    corner_of_data = inset_ax.transData.inverted().transform(inset_ax.transAxes.transform(corner_2))
    big_ax.plot((corner_of_axis[0],corner_of_data[0]),(corner_of_axis[1],corner_of_data[1]),
            ls=':',color=inset_spine_color)

def nmr_color_scaled(file,fig,axes,colormap,true_colormap,extent=None,shape=None,interp_map=None,
                     inset_start=None,inset_width=None,inset_height=None,
                     colormap_spine_color='k',corner_1=[0,0],corner_2=[1,0]):
    data = nmr_read(file)
    data = data.reshape(shape)
    
    
    #Create the scaled and colorized data for plotting
    scaled,colorized = colorize(data,interp_map,true_colormap)
    junk,colorized_inset = colorize(data,interp_map,colormap)
    
    
    #Working in the right axis
    ax = axes.flatten()[1]
    
    
    #Plot the scaled data (this is what goes into the outlier analysis)
    im=ax.imshow(scaled,cmap='Greys_r',
                 extent=extent,aspect='auto',vmin=0#,vmax=1
                )
    nmr_annotate(im,ax,label='Greyscale brightness',cb_ticks=[0,scaled.max()])
    
    #Create the subplot
    if inset_start is not None:
        inset_start_x,inset_start_y = inset_start
        ax_sub = mplax.zoomed_inset_axes(ax,zoom=3,loc=1)
        im=ax_sub.imshow(scaled,#cmap='RdBu',
                         cmap='Greys_r',
                         extent=extent,aspect='auto',vmin=0,#vmax=1,
                        )
        #Set limits of subplot to the size of the desired viewing window
        ax_sub.set_xlim(inset_start_x,inset_start_x+inset_width)
        ax_sub.set_ylim(inset_start_y,inset_start_y+inset_height)

        inset_spine_color='w'
        for spine in ax_sub.spines.items():
            #print(spine)
            spine[1].set_color(inset_spine_color)

        ax_sub.set_xticks([])
        ax_sub.set_yticks([])
        fig.canvas.draw()
        better_mark_inset(ax,ax_sub,corner_1,corner_2,inset_spine_color)

    #Working in the left axis
    ax = axes.flatten()[0]
    
    
    #Plot the scaled data (this is what goes into the outlier analysis)
    im=ax.imshow(colorized_inset,cmap=colormap,
                 extent=extent,aspect='auto',vmin=0#,vmax=1
                )
    cb = nmr_annotate(im,ax,cb_ticks=np.array(interp_map[1])[[0,1,3,4]],
                      label='NMR field intensity (% of max)')
    cb.set_ticklabels(['{:.0f}'.format(interp_map[0][0]),
                       '{:.1f}'.format(interp_map[0][1]),
                       '{:.1f}'.format(interp_map[0][3]),
                       '{:.0f}'.format(interp_map[0][4])])
    
    #Create the subplot
    if inset_start is not None:
        ax_sub = mplax.zoomed_inset_axes(ax,zoom=3,loc=1)
        im=ax_sub.imshow(colorized_inset,#cmap='RdBu',
                         extent=extent,aspect='auto',vmin=0,#vmax=1,
                        )
        #Set limits of subplot to the size of the desired viewing window
        ax_sub.set_xlim(inset_start_x,inset_start_x+inset_width)
        ax_sub.set_ylim(inset_start_y,inset_start_y+inset_height)

        inset_spine_color='k'
        for spine in ax_sub.spines.items():
            #print(spine)
            spine[1].set_color(colormap_spine_color)

        ax_sub.set_xticks([])
        ax_sub.set_yticks([])
        fig.canvas.draw()
        better_mark_inset(ax,ax_sub,corner_1,corner_2,colormap_spine_color)

def nmr_annotate(im,ax,cb_ticks=[0,1],label=None):
    ax.invert_yaxis()
    ax.set_xlabel('$^1$H shift (ppm)',size=15)
    ax.set_ylabel('$^{13}$C shift (ppm)',size=15)
    ax.set_xticks([1,0])
    ax.set_yticks([15,20,25])
    cb=plt.colorbar(im,ax=ax)
    cb.set_ticks(cb_ticks)
    cb.ax.tick_params(labelsize=15)
    if label: cb.set_label(label,size=15)
    return cb