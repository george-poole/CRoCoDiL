# kills zombie Python processes left running due a bug in VSCode / Jupyter
echo "Python processes before:"
ps -ef | grep python 
echo "Killing python processes..." 
pkill -f "/Users/George/miniconda3/envs/lucifex/bin/python -m ipykernel_launcher"
echo "Python processes after:"
ps -ef | grep python  