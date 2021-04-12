for epsilon in 0 0.1 0.001 0.0001 
do
for J in 0 1    
do 
for hmax in 1 5 10 20
do
 

echo "#!/usr/local_rwth/bin/zsh




#SBATCH --job-name=Quench_stark_$L.$g0  # The job name.
#SBATCH -c 4                    # The number of cpu cores to use.
#SBATCH --time=47:59:00              # The time the job will take to run.
#SBATCH --mem-per-cpu=3200M        # The memory the job will use per cpu core.
#SBATCH --account=rwth0722
#SBATCH --mail-type=NONE
#SBATCH --mail-user=passetti@physik.rwth-aachen.de


#SBATCH --output=$HOME/output/ws.$L.$g0
#

### Change to the work directory
cd $HOME/Spinless_Boson/QI
### Execute your application
MKL_NUM_THREADS=4
export MKL_NUM_THREADS


export PYTHONPATH=$HOME/TeNPy
python3 main.py $epsilon $J $hmax


echo \$?

date
" >$HOME/Spinless_Boson/Quench_stark/Bash/addBash_STARK$L.$g0

sbatch <$HOME/Spinless_Boson/Quench_stark/Bash/addBash_STARK$L.$g0

done
done
done

