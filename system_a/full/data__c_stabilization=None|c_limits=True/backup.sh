DIR_NAME=${1:-"./figures/"}
BACKUP_NAME=${2:-"backup"}
DEST="${DIR_NAME}/${BACKUP_NAME}"

rsync -a $DIR_NAME $DEST \
 --exclude '*.pickle' \
 --exclude '*.pdf' \
 --exclude "${BACKUP_NAME}/"
