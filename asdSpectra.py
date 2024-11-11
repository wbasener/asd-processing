import specdal
from spectral import envi
from spectral import resampling
import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import tkinter as tk
from tkinter import filedialog
from numpy import where
from numpy import isnan
from numpy import zeros
from PIL import Image
import os


def plot_asd_spectra():
    files = filedialog.askopenfilenames(parent=parent, 
                                        title='Choose .asd files to plot')
    if files=='':
        return

    fig, ax = plt.subplots()
    for fname in files:
        s = specdal.Spectrum(filepath=fname)
        wl = s.measurement.index
        
        # Basener's fix around detector plane edges:
        # note: This is similar to the jump_correct method is SpecDal, but shifts to match 1st derivative, 
        # and best performace comes from using this then jump_correc
        # Fix 1: shift 0<wl<1000 range up/down to smooth jump at 1000
        i1 = where(wl==1000)[0][0]
        if not isnan(s.measurement.iloc[i1]):
            dp = ( ((s.measurement.iloc[i1+1]-s.measurement.iloc[i1+2]) + (s.measurement.iloc[i1-1]-s.measurement.iloc[i1]))/2 )
            d1 = (s.measurement.iloc[i1]-s.measurement.iloc[i1+1])
            s.measurement.iloc[:(i1+1)] = s.measurement.iloc[:(i1+1)] + dp - d1
        # Fix 2: shift 1800<wl<2500 range up/down to smooth jump at 1800
        i2 = where(wl==1800)[0][0]
        if not isnan(s.measurement.iloc[i2]):
            dp = ( ((s.measurement.iloc[i2+1]-s.measurement.iloc[i2+2]) + (s.measurement.iloc[i2-1]-s.measurement.iloc[i2]))/2 )
            d1 = (s.measurement.iloc[i2]-s.measurement.iloc[i2+1])
            s.measurement.iloc[(i2+1):] = s.measurement.iloc[(i2+1):] + d1 - dp
        
        plt.plot(wl,
                s.measurement.iloc[:],
                label=os.path.basename(fname))

    # Show the major grid lines with dark grey lines
    plt.grid(True, which='major', color='#666666', linestyle='-')
    # Show the minor grid lines with very faint and almost transparent grey lines
    plt.minorticks_on()
    plt.grid(True, which='minor', color='#999999', linestyle='-', alpha=0.2)
    plt.legend()
    plt.title('ASD Spectra')
    plt.xlabel('Wavelength')
    plt.ylabel('Reflectance')
    plt.show()


def create_jpg_quicklooks():
    dir_name = filedialog.askdirectory(parent=parent, 
                                  title='Choose a directory with envi hyperspectral images to create RGB jpgs')
    if dir_name=='':
        return
    for im_fname in os.listdir(dir_name):
        if im_fname.endswith('.hdr'):
            fname_full_hdr = dir_name+'\\'+im_fname
            print('Creating JPG from '+fname_full_hdr+'.')
            try:
                fname_full_jpg = fname_full_hdr[:-4]+'_RGB.jpg'
                im = envi.open(fname_full_hdr)
                RGB_im = np.zeros((im.nrows,im.ncols,3)).astype(np.uint8)
                RGB_wl = np.array((630,532,465))
                wl = np.asarray(im.bands.centers)
                if min(wl)<10:
                    RGB_wl = RGB_wl/1000.0
                for i in range(3):
                    band_idx = np.argmin(np.abs(wl-RGB_wl[i]))
                    band_im = im.read_band(band_idx).astype(float)
                    low = np.percentile(band_im, 2)
                    high = np.percentile(band_im, 98)
                    band_im = np.clip(band_im, low, high)
                    band_im = band_im-low
                    band_im = 255*band_im/(high-low)
                    RGB_im[:,:,i] = band_im.astype(np.uint8)
                im_data = Image.fromarray(RGB_im, 'RGB') 
                im_data.save(fname_full_jpg) 
            except:
                print('Unable to create JPG. This is likely because the hyperspectra data for the located header could npot be found.')
    

def asd_2_sli():
    files = filedialog.askopenfilenames(parent=parent, 
                                        title='Choose .asd files to include in spectral library')
    if files=='':
        return
    fname_save = filedialog.asksaveasfile(parent=parent, 
                                     title='Save as Spectral Library', 
                                     filetypes = [('Spectral Library', '*.sli')], 
                                     defaultextension = ('Spectral Library', '*.sli'))
    if fname_save=='':
        return

    spec_names = []
    nSpec = len(files)
    s = specdal.Spectrum(filepath=files[0])
    nBands = len(s.measurement[:])
    spectra_arr = zeros((nSpec,nBands))
    wl = s.measurement.index
    fwhm = []
    for val in wl:
        width=8
        if val < 1000:
            width=3
        fwhm.append(width)
    
    for idx, fname in enumerate(files):
        s = specdal.Spectrum(filepath=fname)
        spec_names.append(s.name)
        
        # Basener's fix around detector plane edges:
        # note: This is similar to the jump_correct method is SpecDal, but shifts to match 1st derivative, 
        # and best performace comes from using this then jump_correc
        # Fix 1: shift 0<wl<1000 range up/down to smooth jump at 1000
        i1 = where(wl==1000)[0][0]
        if not isnan(s.measurement.iloc[i1]):
            dp = ( ((s.measurement.iloc[i1+1]-s.measurement.iloc[i1+2]) + (s.measurement.iloc[i1-1]-s.measurement.iloc[i1]))/2 )
            d1 = (s.measurement.iloc[i1]-s.measurement.iloc[i1+1])
            s.measurement.iloc[:(i1+1)] = s.measurement.iloc[:(i1+1)] + dp - d1
        # Fix 2: shift 1800<wl<2500 range up/down to smooth jump at 1800
        i2 = where(wl==1800)[0][0]
        if not isnan(s.measurement.iloc[i2]):
            dp = ( ((s.measurement.iloc[i2+1]-s.measurement.iloc[i2+2]) + (s.measurement.iloc[i2-1]-s.measurement.iloc[i2]))/2 )
            d1 = (s.measurement.iloc[i2]-s.measurement.iloc[i2+1])
            s.measurement.iloc[(i2+1):] = s.measurement.iloc[(i2+1):] + d1 - dp
        
        spectra_arr[idx,:] = s.measurement
        
    # save the collection as an ENVI spectral library
    header = {}
    header['spectra names'] = spec_names
    header['wavelength'] = list(wl)
    header['fwhm'] = fwhm
    lib = envi.SpectralLibrary(spectra_arr, header, [])
    lib.save(fname_save.name[:-4])
    

def compute_confusers():    
    nSpecOut = 1000
    fname_sli_tgt = filedialog.askopenfilename(parent=parent, 
                                        title='Choose spectral library with target spectra', 
                                        filetypes = [('Spectral Library', '*.sli')])    
    if fname_sli_tgt=='':
        return
    fname_sli_lib = filedialog.askopenfilename(parent=parent, 
                                        title='Choose spectral library with general spectra', 
                                        filetypes = [('Spectral Library', '*.sli')])     
    if fname_sli_lib=='':
        return

    # open the files
    lib_tgt = envi.open(fname_sli_tgt.rpartition('.')[0]+'.hdr') 
    lib = envi.open(fname_sli_lib.rpartition('.')[0]+'.hdr') 
    # get the number of spectra in each library
    nSpecTgt = lib_tgt.spectra.shape[0]    
    nSpec = lib.spectra.shape[0]
    nBands = lib.spectra.shape[1]
    nSpecOut = np.min([nSpecOut,nSpec])
    
    # comput the pairwise correlations
    means_tgt = np.mean(lib_tgt.spectra, axis=1)
    std_tgt = np.std(lib_tgt.spectra, axis=1)
    spec_tgt_dnormalize = np.zeros(lib_tgt.spectra.shape)
    for i in range(nSpecTgt):
        spec_tgt_dnormalize[i,:] = (lib_tgt.spectra[i,:] - means_tgt[i] )/std_tgt[i]   
    means_lib = np.mean(lib.spectra, axis=1)   
    std_lib  = np.std(lib.spectra, axis=1)
    spec_normalize = np.zeros(lib.spectra.shape)
    for i in range(nSpec):
        spec_normalize[i,:] = (lib.spectra[i,:] - means_lib[i])/std_lib[i]
    cors = np.matmul(spec_tgt_dnormalize,spec_normalize.T)/nBands

    max_correlation = np.max(cors, axis=0)
    sorted_indices = np.flip(np.argsort(max_correlation))
    spectra_out = np.zeros((nSpecOut,nBands))
    spectra_out_pairs = np.zeros((nSpecOut*2,nBands))
    spec_names_out = []
    spec_names_out_pairs = []
    for count, idx in enumerate(sorted_indices[:nSpecOut]):
        # spectra for confuser library
        spectra_out[count,:] = lib.spectra[idx,:]
        spec_names_out.append(lib.names[idx])
        # spectra for confuser-target pairs library
        spectra_out_pairs[count*2,:] = lib.spectra[idx,:]
        spec_names_out_pairs.append(lib.names[idx])
        idx_tgt = np.argmax(cors[:,idx])
        spectra_out_pairs[count*2+1,:] = lib_tgt.spectra[idx_tgt,:]
        spec_names_out_pairs.append(lib_tgt.names[idx_tgt])
    
    # save the collection as an ENVI spectral library
    fname_save = filedialog.asksaveasfile(parent=parent, 
                                    title='Save as Spectral Library', 
                                    filetypes = [('Spectral Library', '*.sli')], 
                                    defaultextension = ('Spectral Library', '*.sli'))
    # save the spectral library of confusers
    header = {}
    header['spectra names'] = spec_names_out
    header['wavelength'] = list(lib.bands.centers)
    lib = envi.SpectralLibrary(spectra_out, header, [])
    lib.save(fname_save.name[:-4])
    # save the spectral library of confuser-target pairs
    header = {}
    header['spectra names'] = spec_names_out_pairs
    header['wavelength'] = list(lib.bands.centers)
    lib = envi.SpectralLibrary(spectra_out_pairs, header, [])
    lib.save(fname_save.name[:-4]+'_pairs')
    
    
    

def resample_library():  
    fname_sli_data = filedialog.askopenfilename(parent=parent, 
                                        title='Choose spectral library with spectral data', 
                                        filetypes = [('Spectral Library', '*.sli')])  
    if fname_sli_data=='':
        return  
    fname_data_bands = filedialog.askopenfilename(parent=parent, 
                                        title='Choose spectral library or image with resmpling wavelengths', 
                                        filetypes = [('Spectral Library', '*.sli'),('ENVI Spectrall Image', '*.*')])  
    if fname_data_bands=='':
        return  

    # open the files
    sli_data = envi.open(fname_sli_data.rpartition('.')[0]+'.hdr') 
    data_bands = envi.open(fname_data_bands.rpartition('.')[0]+'.hdr') 
    # convert library to image wavelength units
    if np.mean(data_bands.bands.centers) < 100:
        # image units are micrometers
        if np.mean(sli_data.bands.centers) < 100:
            wl_scale = 1.
        else:
            wl_scale = 1/1000.
    else:
        # image units are nanometers
        if np.mean(sli_data.bands.centers) < 100:
            wl_scale = 1000.
        else:
            wl_scale = 1.
    sli_data.bands.centers = list(wl_scale*np.asarray(lib.bands.centers))
    # resample library to image wavelengths
    # compute resampling matrix
    # compute the library bandwidths if they are not present
    #  - bandwiths should be distance to nearest band, and at least 3nm
    if sli_data.bands.bandwidths is None:
        sli_data.bands.bandwidths = []
        for i in range(len(sli_data.bands.centers)):
            if i == 0:
                w = np.abs(sli_data.bands.centers[1] - sli_data.bands.centers[0])
            elif i == (len(sli_data.bands.centers)-1):                    
                w = np.abs(sli_data.bands.centers[i] - sli_data.bands.centers[i-1])
            else:
                w = np.min([np.abs(sli_data.bands.centers[i] - sli_data.bands.centers[i-1]),
                            np.abs(sli_data.bands.centers[i+1] - sli_data.bands.centers[i]) ])
            w = np.max([w,3])
            sli_data.bands.bandwidths.append(w)
    if data_bands.bands.bandwidths is None:
        data_bands.bands.bandwidths = []
        for i in range(len(data_bands.bands.centers)):
            if i == 0:
                w = np.abs(data_bands.bands.centers[1] - data_bands.bands.centers[0])
            elif i == (len(data_bands.bands.centers)-1):                    
                w = np.abs(data_bands.bands.centers[i] - data_bands.bands.centers[i-1])
            else:
                w = np.min([np.abs(data_bands.bands.centers[i] - data_bands.bands.centers[i-1]),
                            np.abs(data_bands.bands.centers[i+1] - data_bands.bands.centers[i]) ])
            w = np.max([w,3])
            data_bands.bands.bandwidths.append(w)
    # compute the resampling parameters
    resample = resampling.BandResampler(sli_data.bands.centers, 
                                        data_bands.bands.centers,
                                        fwhm1=sli_data.bands.bandwidths,
                                        fwhm2=data_bands.bands.bandwidths)
    # resample
    sli_data.spectra = np.matmul(sli_data.spectra,resample.matrix.T)
    sli_data.bands.centers = data_bands.bands.centers
    
    fname_save = filedialog.asksaveasfile(parent=parent, 
                                    title='Save as Spectral Library', 
                                    filetypes = [('Spectral Library', '*.sli')], 
                                    defaultextension = ('Spectral Library', '*.sli'))
    sli_data.save(fname_save)    
    
    
    
    
    
    
    
    
    
    
parent = tk.Tk()
frame = tk.Frame(parent)
frame.pack()

plot_spec_btn = tk.Button(frame, 
                   text="Plot ASD Spectra", 
                   command=plot_asd_spectra
                   )

plot_spec_btn.pack(side=tk.LEFT)

create_jpg_quicklooks_btn = tk.Button(frame, 
                   text="Create JPGs from HSI", 
                   command=create_jpg_quicklooks
                   )
create_jpg_quicklooks_btn.pack(side=tk.LEFT)

asd_2_sli_btn= tk.Button(frame, 
                   text="Create Spectral Library", 
                   command=asd_2_sli
                   )

asd_2_sli_btn.pack(side=tk.LEFT)

compute_confusers_btn= tk.Button(frame, 
                   text="Compute Confusers", 
                   command=compute_confusers
                   )

compute_confusers_btn.pack(side=tk.LEFT)

resample_library_btn= tk.Button(frame, 
                   text="Resample a library", 
                   command=resample_library
                   )

resample_library_btn.pack(side=tk.LEFT)

exit_button = tk.Button(frame,
                   text="Exit",
                   fg="green",
                   command=parent.destroy)
exit_button.pack(side=tk.RIGHT)

parent.mainloop()