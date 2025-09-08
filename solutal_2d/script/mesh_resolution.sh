N_PROC=$1 
N_STOP=${2:-15000}
NX=(${3:-120 140 160})
NY=(${4:-120 140 160})

echo "n_proc = $N_PROC"
echo "n_stop = $N_STOP"
echo "Nx = ${NX[@]}"
echo "Ny = ${NY[@]}"

function python_simulate {
    python simulate.py \
    --Nx $1 \
    --Ny $2 \
    --n_stop $3 \
    --write_step 0.01 \
    --write_file '("FunctionSeries", "ConstantSeries")' \
    --dir_base '"./appendix/mesh_resolution"' \
    --dir_labels '("Nx", "Ny", "Ra", "Da", "sr")' \
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
    --c_stabilization None \
    --c_limits Ellipsis
}
export -f python_simulate

if [ $N_PROC -eq 0 ]
then
# single job
python_simulate ${NX[0]} ${NY[0]} $N_STOP
else
# parallel jobs
parallel --link -j $N_PROC python_simulate {1} {2} $N_STOP ::: ${NX[@]} ::: ${NY[@]}
fi