PNG_ONLY=false

while [[ "$1" == --* ]]; do
    case "$1" in
        --png_only)
        PNG_ONLY=true
        shift
        ;;
    esac
done

DIR_NAME=${1:-"./figures/"}
BACKUP_NAME=${2:-"backup"}
DEST="${DIR_NAME}/${BACKUP_NAME}"

if $PNG_ONLY; then
    rsync -a $DIR_NAME $DEST --exclude "${BACKUP_NAME}/" --include '*.png'
else
    rsync -a $DIR_NAME $DEST --exclude "${BACKUP_NAME}/"
fi
