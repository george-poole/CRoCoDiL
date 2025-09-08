N_PROC=$1 
N_STOP=${2:-15000}
RA=(${3:-100 500 1000})
DA=(${4:-1 100 500 1000})
SR=(${5:-0.1 0.05 0.2})

echo "n_proc = $N_PROC"
echo "n_stop = $N_STOP"
echo "Ra = ${RA[@]}"
echo "Da = ${DA[@]}"
echo "sr = ${SR[@]}"

function python_simulate {
    python simulate.py \
    --Ra $1 \
    --Da $2 \
    --sr $3 \
    --n_stop $4 \
    --write_step 0.01 \
    --write_file '("FunctionSeries", "ConstantSeries")' \
    --dir_base '"./data' \
    --dir_labels '("Ra", "Da", "sr")' \
    --dir_timestamp True \
    --dt_init 0.000001 \
    --n_init 10 \
    --t_stop 120.0 \
    --texec True \
    --Nx 160 \
    --Ny 200 \
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
python_simulate ${RA[0]} ${DA[0]} ${SR[0]} $N_STOP
else
# parallel jobs
parallel -j $N_PROC python_simulate {1} {2} {3} $N_STOP ::: ${RA[@]} ::: ${DA[@]} ::: ${SR[@]}
fi