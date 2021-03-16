
for g0 in 0 0.2 0.4 0.6 0.8 1
do 
for chi in 100 
do
for dt in 0.005 
do
for V in 4
do
echo "#!/usr/local_rwth/bin/zsh




#SBATCH --job-name=rampV_CDW.$g0  # The job name.
#SBATCH -c 4                    # The number of cpu cores to use.
#SBATCH --time=47:59:00              # The time the job will take to run.
#SBATCH --mem-per-cpu=3200M        # The memory the job will use per cpu core.
#SBATCH --account=rwth0722
#SBATCH --mail-type=NONE
#SBATCH --mail-user=passetti@physik.rwth-aachen.de


#SBATCH --output=$HOME/output/ws.$g0.$chi.$dt.$V
#

### Change to the work directory
cd $HOME/Spinless_Boson/CDW
### Execute your application
MKL_NUM_THREADS=4
export MKL_NUM_THREADS


export PYTHONPATH=$HOME/TeNPy
python3 tev_CDW.py $g0 $chi $dt $V


echo \$?

date
" >$HOME/Spinless_Boson/CDW/Bash/addBash_Rv$g0

sbatch <$HOME/Spinless_Boson/CDW/Bash/addBash_Rv$g0

done
done 
done
done
