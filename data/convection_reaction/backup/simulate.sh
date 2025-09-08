NPROC=2 

CR=0.0
H0=0.9
EPSILON=0.01

NX=150
NY=150

T_STOP=120
N_STOP=12000
DT_INIT=0.000001
N_INIT=10

if [ $NPROC -eq 0 ]
then
# serial job
RA=500
DA=500
SR=0.1
python simulate.py \
 --write_step 0.01 --write_file '("FunctionSeries", "ConstantSeries")' \
 --dir_base '"./data/convection_reaction"' --dir_timestamp True \
 --dir_labels '("Ra", "Da", "Sr")' \
 --dt_init $DT_INIT --n_init $N_INIT --t_stop $T_STOP --n_stop $N_STOP --texec True \
 --Nx $NX --Ny $NY --c_stabilization None --c_limits '(True, True)' \
 --Ra $RA --Da $DA --Sr $SR --Cr $CR --h0 $H0 --epsilon $EPSILON
elif [ $NPROC -lt 0 ]
then
# test parallel job
let NPROC=-$NPROC
let N_STOP=2*$N_INIT
parallel -j $NPROC "python simulate.py \\
--write_step 0.01 --write_file '(\"FunctionSeries\", \"ConstantSeries\")' \\
--dir_base '\"./data/convection_reaction\"' --dir_timestamp True --dir_tag '\"test\"' \\
--dir_labels '(\"Ra\", \"Da\", \"Sr\")' \\
--dt_init $DT_INIT --n_init $N_INIT --t_stop $T_STOP --n_stop $N_STOP --texec True \\
--c_stabilization None --c_limits '(True, True)' \\
--Nx $NX --Ny $NY --h0 $H0 --epsilon $EPSILON \\
--Ra {1} --Da {2} --Sr {3}" ::: 100 200 ::: 300 400 ::: 0.1 0.2
else
# full parallel job
parallel -j $NPROC "python simulate.py \\
--write_step 0.01 --write_file '(\"FunctionSeries\", \"ConstantSeries\")' \\
--dir_base '\"./data/convection_reaction\"' --dir_timestamp True \\
--dir_labels '(\"Ra\", \"Da\", \"Sr\")' \\
--dt_init $DT_INIT --n_init $N_INIT --t_stop $T_STOP --n_stop $N_STOP --texec True \\
--c_stabilization None --c_limits '(True, True)' \\
--Nx $NX --Ny $NY --h0 $H0 --epsilon $EPSILON \\
--Ra {1} --Da {2} --Sr {3}" ::: 10 100 500 1000 ::: 1 10 100 1000 ::: 0.1 0.2 0.05
fi