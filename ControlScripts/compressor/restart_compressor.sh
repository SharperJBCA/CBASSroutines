# Script to be run to restart the compressor

pypath=/home/cbassuser/monitor/source/pyserial-2.7

# want to check the PYTHONPATH is looking to pyserial-2.7
if [ -z "$PYTHONPATH" ]
then
    # PYTHONPATH not set so set it
    export PYTHONPATH="$pypath"
else
    # PYTHONPATH set, but does it contain pyserial already?
    if [[ $PYTHONPATH == *"$pypath"* ]];
    then
	echo "PYTHONPATH set"
    else
	export PYTHONPATH="$pypath:$PYTHONPATH"
    fi
fi
 
# now call the python script
python compressor.py --restart
