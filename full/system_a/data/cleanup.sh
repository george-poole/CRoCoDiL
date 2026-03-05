CONFIRM=false

while [[ "$1" == --* ]]; do
    case "$1" in
        --confirm)
        CONFIRM=true
        shift
        ;;
    esac
done

GLOB=$1
DIR_NAME=${2:-'*'}

EXTS=("xdmf" "h5")
FILES=()

for ext in "${EXTS[@]}"
    do 
        FOUND=($(find . -name "$GLOB.$ext" -path "./$DIR_NAME"))
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
        echo Cleaning up $i 
        rm $i
    done