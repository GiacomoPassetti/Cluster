for Delta in 0.5 1 2
do 

echo "#!/usr/local_rwth/bin/zsh




#SBATCH --job-name=XXZ_prop$Delta  # The job name.
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
python3 XXZ_SPIN_PROP.py $Delta


echo \$?

date
" >$HOME/quantum_information/Bash/erx$Delta

sbatch <$HOME/quantum_information/Bash/erx$Delta

done


