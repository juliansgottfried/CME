ml julia

sbatch \
    --nodes=1 \
    --ntasks-per-node=10 \
    --mem=8G \
    --time=00:15:00 \
    --output=/scratch/users/jgottf/CME/output/%j.out \
    --error=/scratch/users/jgottf/CME/output/%j.out \
    --partition=normal,hns \
    --mail-type=ALL \
    --mail-user=juliansgottfried@gmail.com \
    /scratch/users/jgottf/CME/run.jl
