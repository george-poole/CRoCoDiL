# Kills zombie Python processes left running by IPython notebooks, 
# without killing any Python processes currently running from a Python script.
#
# `pkill -f /Users/George/miniconda3/envs/lucifex/bin/python -m ipykernel_launcher`

CONFIRM=false

while [[ "$1" == --* ]]; do
    case "$1" in
        --confirm)
        CONFIRM=true
        shift
        ;;
    esac
done

PKILL_ARG=${1:-"-m ipykernel_launcher"}
USER_ENV=${2:-"/Users/George/miniconda3/envs/lucifex/bin/python"}

echo "Python processes:"
ps -eo pid,ni,args | grep python 

if ! $CONFIRM; then
    echo "Exiting unconfirmed run"
    exit
fi

echo ""
echo "Killing zombie python processes..." 
pkill -f "$USER_ENV $PKILL_ARG"

echo ""
echo "Python processes not killed:"
ps -eo pid,ni,args | grep python  