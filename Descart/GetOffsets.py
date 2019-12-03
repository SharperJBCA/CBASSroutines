import DescartLoadFuncs
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

def GetOffsets(filename, offsetdir = 'offsets/', nodes=None, headerScanName='NAXIS2', offset = 1000):
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
        
    code = filename.split('/')[-1]
    hdu = fits.open(filename)

    mask = hdu[1].data['FLAG'].astype(bool) | hdu[1].data['DAYFLAG'].astype(bool)

    output = {}
    for node, key in nodes.items():
        offsetfiles = glob.glob('offsets/{}*_{}.off'.format(code,node))
        
        if len(offsetfiles) == 0:
            continue

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
        data_cube    = data[:N].reshape((Nscans,mask.size//Nscans))


        last = 0
        lastOffset = 0
        for i in range(Nscans):
            start = hdu[2].data['START'][i]
            end = hdu[2].data['END'][i]

            ntod = end-start + 1
            na= ntod//offset 
            N = na * offset
            d = np.arange(N,dtype=int)
            m = mask[start:start+N]
            d = d[~m]
            
            Noffsets = d.size//offset 
            N2 = d.size//offset * offset
            d = d[:N2]
            offset_indices[start:start+N2] = d//offset + lastOffset

            
            data_test[start:start+N2] = data[start:start+N2] - np.nanmean(data[start:start+N2])
            lastOffset = np.max(offset_indices)+1
            last += N2


        
        #print(key)
        if not 'offsets' in output:
            output['offsets'] = {}
        output['offsets'][key] = offsets[offset_indices]
        output['offsets'][key][offset_indices==-1]=np.nan
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
