N_PROC=$1 
N_STOP=${2:-'15000'}
CELL=(${3:-'quadrilateral' 'left' 'right' 'left_right'})
NX=${4:-160}
NY=${5:-160}

echo "n_proc = $N_PROC"
echo "n_stop = $N_STOP"
echo "cell = ${CELL[@]}"
echo "Nx = $NX"
echo "Ny = $NY"

function python_simulate {
    python simulate.py \
    --cell $1 \
    --n_stop $2 \
    --write_step 0.01 \
    --write_file '("FunctionSeries", "ConstantSeries")' \
    --dir_base '"./appendix/mesh_resolution"' \
    --dir_labels '("cell", "Nx", "Ny")' \
    --dir_timestamp True \
    --dt_init 0.000001 \
    --n_init 10 \
    --t_stop 100.0 \
    --texec True \
    --Ra 1000 \
    --Da 1000 \
    --sr 0.1 \
    --cr 0.0 \
    --h0 0.9 \
    --epsilon 0.01 \
    --Nx $NX \
    --Ny $NY \
    --c_stabilization None \
    --c_limits Ellipsis
}
export -f python_simulate

if [ $N_PROC -eq 0 ]
then
# single job
python_simulate ${CELL[0]} $N_STOP
else
# parallel jobs
parallel -j $N_PROC python_simulate {1} $N_STOP ::: ${CELL[@]}
fi