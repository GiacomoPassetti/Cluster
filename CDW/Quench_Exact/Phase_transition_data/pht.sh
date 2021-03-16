for Omega in 5 10 20
do
echo "#!/usr/local_rwth/bin/zsh




#SBATCH --job-name=pht  # The job name.
#SBATCH -c 4                    # The number of cpu cores to use.
#SBATCH --time=47:59:00              # The time the job will take to run.
#SBATCH --mem-per-cpu=3200M        # The memory the job will use per cpu core.
#SBATCH --account=rwth0722
#SBATCH --mail-type=NONE
#SBATCH --mail-user=passetti@physik.rwth-aachen.de


#SBATCH --output=$HOME/output/cdw_ex
#

### Change to the work directory
cd $HOME/Spinless_Boson/CDW/Phase_transition
### Execute your application
MKL_NUM_THREADS=4
export MKL_NUM_THREADS


python3 main.py $Omega


echo \$?

date
" >$HOME/Spinless_Boson/CDW/Bash/addBash_ex

sbatch <$HOME/Spinless_Boson/CDW/Bash/addBash_ex

done
