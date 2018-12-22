# Module for photometric data reduction
# Written by Alan Braeley and Christian Ruiz

import numpy as np
import matplotlib.pyplot as plt
import astropy
from astropy.io import fits 
import os 
import scipy.ndimage.interpolation as interp


def filesorter(filename, foldername, fitskeyword_to_check, keyword):
    '''
    Gets a fits file. Reads the header. Reads the type. Checks if the folder we want to put the file in already
    exists. If not, creates a new one. Checks if file type matches the desired output file type. If it does,
    it moves the file into the folder.
    '''
    if os.path.exists(filename):
        pass
    else:
        print(filename + " does not exist or has already been moved.")
        return
    
    header = fits.getheader(filename)
    fits_type = header[keyword]
    
    if os.path.exists(foldername):
        pass
    else:
        print("Making new directory: " + foldername)
        os.mkdir(foldername)
        
    if fits_type == fitskeyword_to_check:
        destination = foldername + '/'
        print("Moving " + filename + " to: ./" + destination + filename)
        os.rename(filename, destination + filename)  
    return

def mediancombine(filelist):
    '''
    Takes a list of files and makes a master median file out of these files.
    '''
    n = len(filelist)
    first_frame_data = fits.getdata(filelist[0])
    imsize_y, imsize_x = first_frame_data.shape
    fits_stack = np.zeros((imsize_y, imsize_x , n), dtype = np.float32) 
    for ii in range(0, n):
        im = fits.getdata(filelist[ii])
        fits_stack[:,:,ii] = im
    med_frame = np.median(fits_stack, axis=2)
    return med_frame
    

def bias_subtract(filename, master_bias):
    '''
    subtracts the bias from the dark frame
    '''

    data = fits.getdata(filename)
    header = fits.getheader(filename)
    dark_bias = data - master_bias
    fits.writeto('b_' + filename, dark_bias, header, clobber = True) # to save the FITS
    return

def dark_subtract(filename, master_dark):
    '''
    subtracts the dark from the flat frame
    '''
    # Your code goes here.
    data = fits.getdata(filename)
    header = fits.getheader(filename)
    
    flat_dark_subtracted = data - master_dark
    
    fits.writeto('d' + filename, flat_dark_subtracted, header, clobber = True) # Also to save the FITS
    return


def norm_combine_flats(filelist, filter):
    '''
    normalizes the bias and dark subtracted flat frames then median combines them.
    '''
    n = len(filelist)
    first_frame_data = fits.getdata(filelist[0])
    imsize_y, imsize_x = first_frame_data.shape
    fits_stack = np.zeros((imsize_y, imsize_x , n), dtype = np.float32) 
    for ii in range(0, n):
        im = fits.getdata(filelist[ii])
        norm_im = im / np.median(im)
        fits_stack[:,:,ii] = norm_im
    med_frame = np.median(fits_stack, axis=2)
    flat_header = fits.getheader(filelist[0])
    fits.writeto('Master_Flat_'+filter+'band.fit', med_frame, flat_header, clobber = True)
    return med_frame



def flat_field(filename, master_flat):
    '''
    subtracts the dark from the flat frame
    '''
    data = fits.getdata(filename)
    header = fits.getheader(filename)
    
    flat_dark_subtracted = data/master_flat
    
    fits.writeto('f' + filename, flat_dark_subtracted, header, clobber = True) # Also to save the FITS
    return



def centroid_shifter(images, x1, x2, y1, y2, bgndx1, bgndx2, bgndy1, bgndy2):
    """
    Returns a list of shifts based on the centroid of a given region containing a star
    """
    bgnd = fits.getdata(images[0])[bgndy1:bgndy2,bgndx1:bgndx2]
    med_bgnd = np.median(bgnd)
    dev = np.std(bgnd)
    three_sig = 3*dev
    count = 0
    
    star_cent = []
    
    for image in images:
        region = fits.getdata(image)[y1:y2,x1:x2]
        region = region - three_sig - med_bgnd
        plt.imshow(region)
        
        x_coord = 0
        y_coord = 0
        pix_sum = 0
        
        #normalize the region that doesn't include the star to 0
        for i in range (region.shape[0]):
            for j in range (region.shape[1]):
                if region[i][j]<0:
                    region[i][j]=0
                    
                pix_sum += region[i][j]
                x_coord += (j+1) * region[i][j]
                y_coord += (i+1) * region[i][j]
                    
        star_cent.append([y_coord/pix_sum,x_coord/pix_sum])
    
    shifts_arr = []
    
    for i in star_cent:
        shifts_arr.append([np.around(star_cent[0][0]-i[0]),np.around(star_cent[0][1]-i[1])])
        print('Shift image ', count, 'by', shifts_arr[count])
        count+=1
    return shifts_arr



def shift_combine(images, sh_arr, pad_size, filter):
    """
    Take in list of images, array of shifts, and padding size. Default padding size set to 50 because that is what I am
    using to test.
    """
    counter = 0
    n = len(images)
    first_frame_data = fits.getdata(images[0])
    imsize_y, imsize_x = first_frame_data.shape
    fits_stack = np.zeros((imsize_y+(2*pad_size), imsize_x+(2*pad_size) , n), dtype = np.float32) 
    
    for image in images:
        
        newIm = np.pad(fits.getdata(image), pad_size, 'constant', constant_values = -0.001)
        newerIm = interp.shift(newIm, sh_arr[counter], cval = -0.001)
        print('Generating shifted image')
        fits.writeto('Shifted_'+str(counter)+'_'+filter+'.fit', newerIm, fits.getheader(images[counter]), clobber = True)
        print('Shifted image successfully made')
        fits_stack[:,:,counter] = newerIm
        counter += 1
        
    med_frame = np.median(fits_stack, axis=2)
    med_header = fits.getheader(images[0])
    plt.imshow(med_frame)
    
    print('Generating median combined image')
    fits.writeto('Med_Combine_'+filter+'.fit', med_frame, med_header, clobber = True)
    print('Generated image Med_Combine_'+filter+'.fit')

