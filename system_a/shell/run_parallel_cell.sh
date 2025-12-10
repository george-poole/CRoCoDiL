N_PROC=$1 
T_STOP=${2:-60.0}
CELL=(${3:-'quadrilateral' 'left' 'right' 'left_right'})
NX=${4:-160}
NY=${5:-160}
N_STOP=${6:-20000}

echo "n_proc = $N_PROC"
echo "t_stop = $T_STOP"
echo "n_stop = $N_STOP"
echo "cell = ${CELL[@]}"
echo "Nx = $NX"
echo "Ny = $NY"

function python_simulate {
    python simulate.py \
    --cell $1 \
    --t_stop $T_STOP \
    --n_stop $N_STOP \
    --write_step 0.01 \
    --write_file '("FunctionSeries", "ConstantSeries")' \
    --dir_base '"./appendix/mesh_cell"' \
    --dir_params '("cell", "Nx", "Ny")' \
    --dir_timestamp True \
    --dt_init 0.000001 \
    --n_init 10 \
    --timing True \
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
export T_STOP
export N_STOP
export -f python_simulate

if [ $N_PROC -eq 0 ]
then
# single job
python_simulate ${CELL[0]}
else
# parallel jobs
parallel -j $N_PROC python_simulate {1} ::: ${CELL[@]}
fi