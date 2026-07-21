ml julia

sbatch \
    --nodes=5 \
    --ntasks-per-node=20 \
    --mem=8G \
    --time=24:00:00 \
    --output=/scratch/users/jgottf/CME/output/%j.out \
    --error=/scratch/users/jgottf/CME/output/%j.out \
    --partition=normal,hns \
    --mail-type=ALL \
    --mail-user=juliansgottfried@gmail.com \
    /scratch/users/jgottf/CME/run.jl
