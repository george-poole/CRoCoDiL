N_PROC=$1 
T_STOP=${2:-120.0}
RA=(${3:-100 500 1000})
DA=(${4:-1 100 500 1000})
SR=(${5:-0.1 0.05 0.2})
N_STOP=${6:-20000}

echo "n_proc = $N_PROC"
echo "t_stop = $T_STOP"
echo "n_stop = $N_STOP"
echo "Ra = ${RA[@]}"
echo "Da = ${DA[@]}"
echo "sr = ${SR[@]}"

function python_simulate {
    python simulate.py \
    --Ra $1 \
    --Da $2 \
    --sr $3 \
    --t_stop $T_STOP \
    --n_stop $N_STOP \
    --write_step 0.01 \
    --write_file '("FunctionSeries", "ConstantSeries")' \
    --dir_base '"./data"' \
    --dir_params '("Ra", "Da", "sr")' \
    --dir_timestamp True \
    --dt_init 0.000001 \
    --n_init 10 \
    --texec True \
    --Nx 160 \
    --Ny 200 \
    --cr 0.0 \
    --h0 0.9 \
    --epsilon 0.01 \
    --c_stabilization None \
    --c_limits Ellipsis
}
export T_STOP
export N_STOP
export -f python_simulate

if [ $N_PROC -eq 0 ]
then
# single job
python_simulate ${RA[0]} ${DA[0]} ${SR[0]}
else
# parallel jobs
parallel -j $N_PROC python_simulate {1} {2} {3} ::: ${RA[@]} ::: ${DA[@]} ::: ${SR[@]}
fi