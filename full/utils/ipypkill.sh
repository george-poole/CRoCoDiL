# Kills zombie Python processes left running by IPython notebooks, 
# without killing any Python processes currently running from a Python script.

DRY=false

while [[ "$1" == --* ]]; do
    case "$1" in
        --dry)
        DRY=true
        shift
        ;;
    esac
done

PKILL_ARG=${1:-"-m ipykernel_launcher"}
USER_ENV=${2:-"/Users/George/miniconda3/envs/lucifex/bin/python"}

echo "Python processes before:"
ps -eo pid,ni,args | grep python 

if $DRY; then
    echo "Exiting dry run"
    exit
fi

echo "Killing python processes..." 
pkill -f "$USER_ENV $PKILL_ARG"
echo "Python processes after:"
ps -eo pid,ni,args | grep python  