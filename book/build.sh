REMOTE=false
DRY=false

while [[ "$1" == --* ]]; do
    case "$1" in
        --remote)
        REMOTE=true
        shift
        ;;
    esac
    case "$1" in
        --dry)
        DRY=true
        shift
        ;;
    esac
done

GLOB=$1
N_STOP=${2:-200}
TARGET_DIR=${3:-""}
NBCONVERT_ARGS=${4:-"--allow-errors"}
BUILD_ARGS=${5:-""}

if [ -z "$TARGET_DIR" ]; then
    DIRS=("system_a" "system_aMu" "system_aTheta" "system_b" "system_c" "system_d" "system_x" "tutorial")
else
    DIRS=($TARGET_DIR)
fi

for dir in "${DIRS[@]}"; do
    unlink $dir
    ln -s "../$dir" $dir 
    ipynb_paths=($(find .. -name "$GLOB.ipynb" -path "../$dir/*"))
    for i in "${ipynb_paths[@]}"; do
        echo Found notebook to execute $i
    done
    if ! $DRY; then
        for ipynb in "${ipynb_paths[@]}"; do
            echo Executing notebook $ipynb 
            export IPYNB_FILE_PATH="$ipynb"
            export N_STOP="$N_STOP"
            echo Beginning execution "$(date)"
            jupyter nbconvert --execute --to notebook --inplace $ipynb $NBCONVERT_ARGS  
            echo Finished execution "$(date)"
        done
    fi    
done

if $DRY; then
    echo "Exiting dry run"
    exit
fi

echo Making gallery ...
python make_gallery.py
echo Gallery made

jupyter-book build . $BUILD_ARGS
ln -sf "./_build/html/index.html" alias.html

echo Making extra html paths ...
python make_html_paths.py
echo Extra html paths made

if $REMOTE; then
    ghp-import -n -p -f ./_build/html
fi

ERROR_KEYWORD="Traceback"
for dir in "${DIRS[@]}"; do
    ipynb_paths=($(grep -rl $ERROR_KEYWORD ./$dir/* --include="*.ipynb"))
    for ipynb in "${ipynb_paths[@]}"; do
        echo ""
        echo "WARNING! Error found in $ipynb"
    echo ""
    done
done