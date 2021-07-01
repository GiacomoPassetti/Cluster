for L in 30
do
for h in 0 1 2 5 10 20 50     
do 
for tmax in 20
do
for dt in 0.005
do
for RI in 1 2 3 4 5 6 7 9 10
do
echo "#!/usr/local_rwth/bin/zsh




#SBATCH --job-name=trs$h$RI  # The job name.
#SBATCH -c 4                    # The number of cpu cores to use.
#SBATCH --time=47:59:00              # The time the job will take to run.
#SBATCH --mem-per-cpu=3200M        # The memory the job will use per cpu core.
#SBATCH --account=rwth0722
#SBATCH --mail-type=NONE
#SBATCH --mail-user=passetti@physik.rwth-aachen.de


#SBATCH --output=$HOME/output/qi
#

### Change to the work directory
cd $HOME/MBL/transitions
### Execute your application
MKL_NUM_THREADS=4
export MKL_NUM_THREADS


export PYTHONPATH=$HOME/TeNPy
python3 trs.py $L $h $tmax $dt $RI


echo \$?

date
" >$HOME/Spinless_Boson/Quench_stark/Bash/addBash_STARK

sbatch <$HOME/Spinless_Boson/Quench_stark/Bash/addBash_STARK

done
done
done
done
done