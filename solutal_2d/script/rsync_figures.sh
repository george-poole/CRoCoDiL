DIR_NAME='./data/'
DEST="${DIR_NAME}/backup"

rsync -a $DIR_NAME $DEST \
 --exclude '*.h5' --exclude '*.xdmf' --exclude '*.npz' --exclude '*.npy' \
 --exclude '*.csv' --exclude '*.ipynb' --exclude '*.pickle' \
 --exclude 'backup/'
