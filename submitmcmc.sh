ml julia

sbatch \
    --ntasks=353 \
    --cpus-per-task=1 \
    --mem-per-cpu=2G \
    --time=00:15:00 \
    --output=/scratch/users/jgottf/CME/output/%j.out \
    --error=/scratch/users/jgottf/CME/output/%j.out \
    --partition=normal,hns \
    --mail-type=ALL \
    --mail-user=juliansgottfried@gmail.com \
    /scratch/users/jgottf/CME/run_mcmc.jl
