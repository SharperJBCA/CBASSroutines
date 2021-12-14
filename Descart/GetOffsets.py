from Descart import DescartLoadFuncs
import numpy as np
from matplotlib import pyplot

from astropy.io import fits
import sys
import glob


def GetChunk(i, nSamples, nScans):
    """
    Get the start and end indices
    """
    step = nSamples//nScans
    
    return step*i, step*(i+1)

def ReadOffsets(offsetfiles):
    """
    Calls wrapper for the offset reading code in descart
    """

    allvals = []
    for offsetfile in offsetfiles:
        count, length = DescartLoadFuncs.countvals(offsetfile)
        values = np.zeros(count) # Not ideal!
        DescartLoadFuncs.wrapper(offsetfile,values)
        allvals += [values]
    allvals = np.concatenate(allvals)
    return allvals

def GetOffsets(filename, offsetdir = 'offsets', otheritems= [],
               nodes=None, headerScanName='NAXIS2', offset = 500):
    """
    For a given cbass FITS file return the destriped offsets

    kwargs
    node - Each output offset filename is formated as: _fitsfilename_scanNo_node.off
    where if you destripe in I only node 1 == I1, node 2 == I2, etc...
    """

    # define nodes for most common case
    if isinstance(nodes, type(None)):
        nodes = {1:'I1',
                 2:'I2',
                 3:'Q1',
                 4:'U1',
                 5:'Q2',
                 6:'U2'}
        
    if not 'reload' in filename:
        filename_parts = filename.split('.fits')
        filename_offset = filename_parts[0]+'_reload.fits'
    else:
        filename_offset = filename

    code = filename_offset.split('/')[-1]
    hdu = fits.open(filename)

    mask = hdu[1].data['FLAG'].astype(bool) | hdu[1].data['DAYFLAG'].astype(bool)

    output = {}
    for node, key in nodes.items():
        search_str = '{}/{}*_{}.off'.format(offsetdir,code,node)
        offsetfiles = np.array(glob.glob(search_str))
        if len(offsetfiles) == 0:
            continue
        scannum = np.array([int(filename.split('_')[-2]) for filename in offsetfiles])
        offsetfiles = offsetfiles[np.argsort(scannum)]
        offsets = ReadOffsets(offsetfiles)

       # if not key in hdu[1].data:
        #    continue

        data = hdu[1].data[key]

        Nscans = hdu[2].header[headerScanName] 
        N = data.size//Nscans * Nscans # Crop to integer number of scans

        data_indices = np.zeros((Nscans,data.size//Nscans))
        mask_cube    = mask[:N].reshape((Nscans,mask.size//Nscans))

        offset_indices = np.zeros(data_indices.size,dtype=int) - 1
        data_test = np.zeros(data_indices.size)
        weights   = np.zeros(data_indices.size)
        data_cube    = data[:N].reshape((Nscans,mask.size//Nscans))


        last = 0
        lastOffset = 0
        for i in range(Nscans):
            start = hdu[2].data['START'][i]
            end   = hdu[2].data['END'][i]

            ntod = end-start #+ 2
            na= ntod//offset 
            N = na * offset
            d = np.arange(N,dtype=int)
            m = mask[start:start+N]

            d = d[~m]
            
            Noffsets = d.size//offset 
            N2 = d.size//offset * offset
            d  = d[:N2]

            #print(start, N2, d.size, lastOffset,offset_indices.size)
            if (start + d.size) <= offset_indices.size:
                offset_indices[start:start+N2] = d//offset + lastOffset
            else:
                N3 = offset_indices.size-start
                offset_indices[start:] = d[:N3]//offset + lastOffset

            if np.max(offset_indices) >= offsets.size:
                print('WARNING: Index exceeds size, why?')
                offset_indices[offset_indices > offsets.size-1] = offsets.size-1

            
            z = data[start:start+N2] - np.nanmean(data[start:start+N2])
            zN2 = z.size//2 * 2
            z[z == 0] = np.nan
            rms = np.nanstd(z[1::2] - z[:-1:2])
            z[np.isnan(z)] = 0

            if (start + d.size) <= offset_indices.size:
                data_test[start:start+N2] = z
                weights[start:start+N2] = 1/rms**2
            else:
                data_test[start:start+N2] = z[:N3]
                weights[start:start+N2] = 1/rms**2

            if False:
                ps1 = np.abs(np.fft.fft(data_test[start:start+N2]))**2
                ps2 = np.abs(np.fft.fft(data_test[start:start+N2]-offsets[offset_indices[start:start+N2]]))**2
                nu  = np.fft.fftfreq(N2,d=1/100.)

                pyplot.subplot(311)
                pyplot.plot(data_test[start:start+N2])
                pyplot.plot(offsets[offset_indices[start:start+N2]])
                pyplot.subplot(312)
                pyplot.plot(nu[1:nu.size//2], ps1[1:nu.size//2])
                pyplot.xscale('log')
                pyplot.yscale('log')
                pyplot.grid()
                pyplot.subplot(313)
                pyplot.plot(nu[1:nu.size//2], ps1[1:nu.size//2])
                pyplot.plot(nu[1:nu.size//2], ps2[1:nu.size//2])
                pyplot.xscale('log')
                pyplot.yscale('log')
                pyplot.grid()
                pyplot.show()
            lastOffset = np.max(offset_indices)+1
            last += N2
        

        if not 'offsets' in output:
            output['offsets'] = {}
        output['offsets'][key] = np.zeros(data_test.size)
        output['offsets'][key][:] = offsets[offset_indices]
        output['offsets'][key][offset_indices == -1] = 0

        if not 'data' in output:
            output['data'] = {}
        output['data'][key] = data_test
        if not 'weights' in output:
            output['weights'] = {}
        output['weights'][key] = weights

        # pyplot.subplot(211)
        # pyplot.plot(offsets)
        # pyplot.title(np.max(offset_indices))
        # pyplot.subplot(212)
        # pyplot.plot( data_test)
        # pyplot.plot(output['offsets'][key])
        # pyplot.show()

    for key in otheritems:
        data = hdu[1].data[key]


        Nscans = hdu[2].header[headerScanName] 
        N = data.size//Nscans * Nscans # Crop to integer number of scans

        data_indices = np.zeros((Nscans,data.size//Nscans))
        mask_cube    = mask[:N].reshape((Nscans,mask.size//Nscans))

        offset_indices = np.zeros(data_indices.size,dtype=int) - 1
        data_test = np.zeros(data_indices.size)
        data_cube    = data[:N].reshape((Nscans,mask.size//Nscans))


        last = 0
        for i in range(Nscans):
            start = hdu[2].data['START'][i]
            end = hdu[2].data['END'][i]

            ntod = end-start + 1
            na= ntod//offset 
            N = na * offset
            d = np.arange(N,dtype=int)
            m = mask[start:start+N]
            d = d[~m]            
            N2 = d.size//offset * offset
            if (start + data.size) <= data_test.size:
                data_test[start:start+N2] = data[start:start+N2]
            else:
                N3 = data_test.size-start
                data_test[start:] = data[start:start+N3]

            last += N2
        
        if not 'data' in output:
            output['data'] = {}
        output['data'][key] = data_test

    hdu.close()

    return output


if __name__ == "__main__":
    
    filename = sys.argv[1]
    

    output = GetOffsets(filename, offsetdir = 'offsets/', nodes=None, headerScanName='NAXIS2', offset = 1000)
    print(output)
    pyplot.plot(output['data']['I1'], label='Data')
    pyplot.plot(output['offsets']['I1'],label='Offsets')
    pyplot.title('{}'.format(filename.split('/')[-1]))
    pyplot.xlabel('Sample')
    pyplot.ylabel('K')
    pyplot.savefig('{}_Full.png'.format(filename.split('/')[-1]))
    pyplot.clf()
    pyplot.plot(output['data']['I1'], label='Data')
    pyplot.plot(output['offsets']['I1'],label='Offsets')
    pyplot.title('{}'.format(filename.split('/')[-1]))
    pyplot.xlabel('Sample')
    pyplot.ylabel('K')
    pyplot.xlim(0,100000)
    pyplot.savefig('{}_zoom.png'.format(filename.split('/')[-1]))
    pyplot.show()
