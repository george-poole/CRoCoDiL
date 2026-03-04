PID=$1
NICE=${2:-19}
renice -n $NICE -p $PID