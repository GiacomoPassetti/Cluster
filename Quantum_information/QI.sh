for L in 20 30
do
for eps in 0.01 
do
for iter in 5
do
for rr in 1 2 3 4
do 
for tmax in 10
do
 

echo "#!/usr/local_rwth/bin/zsh




#SBATCH --job-name=erx$L _$rr   # The job name.
#SBATCH -c 4                    # The number of cpu cores to use.
#SBATCH --time=47:59:00              # The time the job will take to run.
#SBATCH --mem-per-cpu=3200M        # The memory the job will use per cpu core.
#SBATCH --account=rwth0722
#SBATCH --mail-type=NONE
#SBATCH --mail-user=passetti@physik.rwth-aachen.de
#SBATCH --output=$HOME/output/qi.txt


### Change to the work directory
cd $HOME/quantum_information
### Execute your application
MKL_NUM_THREADS=4
export MKL_NUM_THREADS


export PYTHONPATH=$HOME/TeNPy
python3 erx.py $L $eps $iter $rr $tmax


echo \$?

date
" >$HOME/quantum_information/Bash/erx$L$rr

sbatch <$HOME/quantum_information/Bash/erx$L$rr

done
done
done
done
done

