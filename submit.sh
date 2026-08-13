ml julia

sbatch \
    --nodes=1 \
    --ntasks-per-node=1 \
    --mem=20G \
    --time=12:00:00 \
    --output=/scratch/users/jgottf/CME/output/%j.out \
    --error=/scratch/users/jgottf/CME/output/%j.out \
    --partition=normal,hns \
    --mail-type=ALL \
    --mail-user=juliansgottfried@gmail.com \
    /scratch/users/jgottf/CME/run.jl
