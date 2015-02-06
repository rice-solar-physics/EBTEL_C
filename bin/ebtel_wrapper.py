#Name: ebtel_wrapper.py
#Author: Will Barnes
#Date: 1 February 2015

#Description: Wrapper for the EBTEL-C program

#Import modules
import os,os.path
import subprocess
import numpy as np
import matplotlib.pyplot as plt

def plot_ebtel_dem(data_directory,data_file,**kwargs):
    """Plot differential emission measure parameters generated by the EBTEL-C model.
    
    Arguments:
    data_directory -- directory that contains EBTEL-C output files
    data_file -- file that contains EBTEL-C output to be plotted
    
    Optional keyword arguments
    print_fig_filename -- specify filename to print figure to (default: show on screen)
    
    """
    
    #Load DEM data
    data = np.loadtxt(data_directory+data_file)
    #Slice array to get necessary vectors
    temp = data[:,0]
    dem_tr = data[:,1]
    dem_cor = data[:,2]
    dem_tot = data[:,3]
    
    #Set up the figure
    fig = plt.figure(figsize=(10,10))
    ax = fig.gca()
    fs = 18
    ax.plot(temp,dem_tr,label=r'TR')
    ax.plot(temp,dem_cor,label=r'Corona')
    ax.plot(temp,dem_tot,label=r'Total')
    ax.legend()
    ax.set_title(r'EBTEL DEM',fontsize=fs)
    ax.set_xlabel(r'$\log(T_{DEM})$ (K)',fontsize=fs)
    ax.set_ylabel(r'$\log($DEM$)$ (cm$^{-5}$ K$^{-1}$)',fontsize=fs)
    ax.set_xlim([5.5,7.5])
    
    #Check if output filename is specified
    if 'print_fig_filename' in kwargs:
        plt.savefig(kwargs['print_fig_filename'],format='eps',dpi=1000)
    else:
        plt.show()
    

def plot_ebtel(data_directory,data_file,**kwargs):
    """Plot plasma parameters generated by the EBTEL-C model.
    
    Arguments:
    data_directory -- directory that contains EBTEL-C output files
    data_file -- file that contains EBTEL-C output to be plotted
    
    Optional keyword arguments
    print_fig_filename -- specify filename to print figure to (default: show on screen)
    
    """
    
    #Load the data
    data = np.loadtxt(data_directory+data_file)
    #Slice the array to get the vectors we want
    time = data[:,0]
    temp = data[:,1]
    dens = data[:,2]
    temp_apex = data[:,5]
    dens_apex = data[:,6]
    heat = data[:,12]
    
    #Set up the figure
    fig,axes = plt.subplots(3,1,figsize=(15,10))
    #Set the fontsize
    fs = 18
    #Plot the heating
    axes[0].plot(time,heat)
    axes[0].set_ylabel(r'$q$ (erg cm$^{-3}$ s$^{-1}$)',fontsize=fs)
    axes[0].set_title(r'EBTEL Plasma Parameters',fontsize=fs)
    axes[0].set_xlim([time[0],time[-1]])
    #Plot the temperatures
    axes[1].plot(time,temp/10**6,label=r'$T$')
    axes[1].plot(time,temp_apex/10**6,'r--',label=r'$T_a$')
    axes[1].set_ylabel(r'$T$ (MK)',fontsize=fs)
    axes[1].legend(loc=1)
    axes[1].set_xlim([time[0],time[-1]])
    #Plot the densities
    axes[2].plot(time,dens/10**8,label=r'$n$')
    axes[2].plot(time,dens_apex/10**8,'r--',label=r'$n_a$')
    axes[2].set_xlabel(r'$t$ (s)',fontsize=fs)
    axes[2].set_ylabel(r'$n$ (10$^8$ cm$^{-3}$)',fontsize=fs)
    axes[2].set_xlim([time[0],time[-1]])
    axes[2].legend(loc=1)
    
    #Check if output filename is specified
    if 'print_fig_filename' in kwargs:
        plt.savefig(kwargs['print_fig_filename'],format='eps',dpi=1000)
    else:
        plt.show()
         

def run_ebtel(exec_directory,config_directory,**kwargs):
    """Run ebtel executable for a single configuration or a whole directory of configurations
    
    Arguments:
    exec_directory -- path to directory that contains executable
    config_directory -- path to config file directory
    
    Optional keyword arguments:
    config_file -- specific config file (run only one instance of EBTEL)
    
    """
    #Change to executable directory
    subprocess.call(['cd',exec_directory],shell=True)
    
    #Check if we want to run a single file or a whole directory
    if 'config_file' in kwargs:
        #Get full output when running a single config file
        output = subprocess.check_output([exec_directory+'ebtel',config_directory+kwargs['config_file']],shell=True)
    else:
        for name in os.listdir(config_directory):
            if os.path.isfile(config_directory+name):
                #Only get exit code when running many configurations
                output = subprocess.call(['./ebtel',name],shell=True)
    
    #Print the output of the subprocess call
    print output

def print_xml_config(config_dictionary,**kwargs):
    """Print XML configuration file for the EBTEL-C model
    
    Arguments:
    config_dictionary -- dictionary that holds config file inputs
    
    Optional keyword arguments:
    config_file -- specify configuration filename (default: 'ebtel_config.xml')
    
    """
    
    #Check if we have passed a filename
    #If not, pass a default filename
    if 'config_file' in kwargs:
        config_file = kwargs['config_file']
    else:
        config_file = 'ebtel_config.xml'
        
    #Open the file
    f = open(config_file,'w')
    
    #Print necessary header info
    f.write('<?xml version="1.0" ?>\n')
    f.write('<input>\n')

    #Loop through dictionary and print to xml file
    for key in config_dictionary:
        #Print tab delimiter, brackets and keyword
        f.write('\t<')
        f.write(key)
        f.write('>')
        #Check if entry is a list
        #If so print it as a list
        if isinstance(config_dictionary[key],list) or type(config_dictionary[key]).__name__ == 'ndarray':
            #Make temporary list
            temp = config_dictionary[key]
            #Skip to new line
            f.write('\n')
            #Begin loop over list
            for i in range(len(config_dictionary[key])):
                f.write('\t\t<')
                f.write(key+str(i))
                f.write('>')
                f.write(str(temp[i]))
                f.write('</')
                f.write(key+str(i))
                f.write('>\n')
            #Print additional tab to preserve alignment
            f.write('\t')
        else:
            #Print value
            f.write(str(config_dictionary[key]))
        #Close the brackets and print newline
        f.write('</')
        f.write(key)
        f.write('>\n')
    
    #Close the main node of the file
    f.write('</input>')
    
    #Close the file
    f.close()