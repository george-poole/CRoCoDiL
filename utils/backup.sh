DIR_NAME=$1
BACKUP_NAME=${2:-"backup"}
DEST="${DIR_NAME}/${BACKUP_NAME}"

rsync -a $DIR_NAME $DEST \
 --exclude "${BACKUP_NAME}/" \
 --exclude '*.h5' \
 --exclude '*.xdmf' \
 --exclude '*.npz' \
 --exclude '*.npy' \
 --exclude '*.csv' \
 --exclude '*.ipynb' \
 --exclude '*.pickle' \
