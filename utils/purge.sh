CONFIRM=false

while [[ "$1" == --* ]]; do
    case "$1" in
        --confirm)
        CONFIRM=true
        shift
        ;;
    esac
done

GLOB_FILE=$1
GLOB_DIR=${2:-'./system_a/*'}

EXTS=("xdmf" "h5" "pickle")
FILES=()

for ext in "${EXTS[@]}"
    do 
        FOUND=($(find . -name "$GLOB_FILE.$ext" -path "./$GLOB_DIR"))
        FILES=("${FILES[@]}" "${FOUND[@]}")
    done

for file in "${FILES[@]}"
    do 
        echo Found file to cleanup $file 
    done

if ! $CONFIRM; then
    echo "Exiting unconfirmed run"
    exit
fi

for i in "${FILES[@]}"
    do 
        echo Purging $i 
        rm $i
    done