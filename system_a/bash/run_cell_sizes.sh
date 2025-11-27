N_PROC=$1 
T_STOP=${2:-60.0}
NX=(${3:-120 140 160 160})
NY=(${4:-120 140 160 200})
N_STOP=${5:-20000}

echo "n_proc = $N_PROC"
echo "t_stop = $T_STOP"
echo "n_stop = $N_STOP"
echo "Nx = ${NX[@]}"
echo "Ny = ${NY[@]}"

function python_simulate {
    python simulate.py \
    --write_file '("FunctionSeries", "ConstantSeries")' \
    --write_step 0.01 \
    --dir_base '"./appendix/mesh_resolution"' \
    --dir_params '("Nx", "Ny", "Ra", "Da", "sr")' \
    --dir_label "'t_stop=${T_STOP}'" \
    --dir_timestamp True \
    --Nx $1 \
    --Ny $2 \
    --Ra 1000 \
    --Da 1000 \
    --epsilon 0.01 \
    --h0 0.9 \
    --sr 0.1 \
    --cr 0.0 \
    --c_stabilization None \
    --c_limits Ellipsis \
    --t_stop $T_STOP \
    --n_stop $N_STOP \
    --dt_init 0.000001 \
    --n_init 10 \
    --timing True \
}
export T_STOP
export N_STOP
export -f python_simulate


if [ $N_PROC -eq 0 ]
then
# single job
python_simulate ${NX[0]} ${NY[0]}
else
# parallel jobs
parallel --link -j $N_PROC python_simulate {1} {2} ::: ${NX[@]} ::: ${NY[@]}
fi