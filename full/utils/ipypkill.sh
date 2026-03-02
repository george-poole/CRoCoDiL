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

USER_ENV="/Users/George/miniconda3/envs/lucifex/bin/python"
echo "Python processes before:"
ps -ef | grep python 

if $DRY; then
    echo "Exiting dry run"
    exit

fi
echo "Killing python processes..." 
pkill -f "$USER_ENV -m ipykernel_launcher"
echo "Python processes after:"
ps -ef | grep python  