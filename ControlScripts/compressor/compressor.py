#! python 2.7
# --- Restart the C-BASS South compressor
#
# version 1: Stuart Harper, October 2nd 2019
#
import sys
import getopt
from communication import com
import time

def restart():
    """
    Runs through set commands needed to restart the compressor
    """
    print('Restarting compressor...')

    # Initilise the compressor communicator
    h = com()
    h.open()

    h.compressor.disable()
    # add a loop here?
    # notRestarted = True
    # while notRestarted:
    #    notRestarted = h.compressor.enable()
    #    time.sleep(30) # wait 30 seconds before checking again
    # ??
    h.compressor.enable()

    print('Successfully restarted')

def usage():
    usage = """
    :> python compressor.py <options>

    Options:
    --help : Print this output 
    --restart : Restarts compressor
    """
    print(usage)

if __name__ == "__main__":
    try:
        # add options to this option list for extra functions
        optlist = ['help','restart']
        opts, args = getopt.getopt(sys.argv[1:], '', optlist)
    except getopt.GetoptError, err:
        print str(err)
        usage()
        sys.exit(2)

    # add function + option name to this dictionary for it be called
    funcs = {'help':usage, 'restart': restart}
    options = {}
    # Parse the options:
    for opt in opts:
        key = opt[0].split('--')[-1]
        options[key] = key in optlist

    for k,v in options.iteritems():
        if v:
            funcs[k]()
