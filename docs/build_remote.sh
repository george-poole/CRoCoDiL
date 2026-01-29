GLOB=$1
BUILD=$2

bash ./build_local.sh $GLOB $BUILD
ghp-import -n -p -f ./_build/html