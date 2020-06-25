# This tries out the distribution of processes across multiple cores, and
# writing the results without collisions to a single file

import multiprocessing as mp
import numpy as np
import pyfits
import glob
import os
import time

fn = './junk.txt'

def worker(file_name, q):
    # reads in a single FITS file, calculates some quantities, and returns them

    f = pyfits.open(file_name)
    img_before_padding_before_FT = np.copy(f[0].data)

    size_y = np.shape(img_before_padding_before_FT)[0]
    size_x = np.shape(img_before_padding_before_FT)[1]

    # take FFT; no padding for now
    '''
    AmpPE, ArgPE = fft_img(img_before_padding_before_FT).fft(padding=0)
    '''

    # save image to check
    '''
    hdu = pyfits.PrimaryHDU(img_before_padding_before_FT)
    hdulist = pyfits.HDUList([hdu])
    hdu.writeto("log_images/img_seen_prepradding_" + file_name_base + ".fits", clobber=True)
    '''

    string_basename = os.path.basename(file_name)
    string_size_y = str(size_y)
    string_size_x = str(size_x)

    with open(fn, 'rb') as f:
        size = len(f.read())
    list1 = string_basename, string_size_y, string_size_x
    res = ', '.join(list1) # strip extra parentheses etc.
    q.put(res)

    return res

def listener(q):
    '''listens for messages on the q, writes to file. '''

    with open(fn, 'w') as f:
        while 1:
            m = q.get()
            if m == 'kill':
                f.write('killed')
                break
            f.write(str(m) + '\n')
            f.flush()

def main():

    path_stem = "./synthetic_fizeau/trial1_opd_tip_tilt/"

    # get list of file names together
    files_list = glob.glob(path_stem + "*.fits")
    print("Files to analyze: ")
    print(files_list)

    #must use Manager queue here, or will not work
    manager = mp.Manager()
    q = manager.Queue()
    pool = mp.Pool(mp.cpu_count() + 2)

    #put listener to work first
    watcher = pool.apply_async(listener, (q,))

    #fire off workers
    jobs = []
    for i in range(len(files_list)):
        file_name_fits = files_list[i]
        job = pool.apply_async(worker1, (file_name_fits, q))
        jobs.append(job)

    # collect results from the workers through the pool result queue
    for job in jobs:
        job.get()

    #now we are done, kill the listener
    q.put('kill')
    pool.close()
    pool.join()

if __name__ == "__main__":
   main()
