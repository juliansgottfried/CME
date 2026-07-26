ml julia

sbatch \
    --nodes=8 \
    --ntasks-per-node=25 \
    --mem=128G \
    --time=03:00:00 \
    --output=/scratch/users/jgottf/CME/output/%j.out \
    --error=/scratch/users/jgottf/CME/output/%j.out \
    --partition=normal,hns \
    --mail-type=ALL \
    --mail-user=juliansgottfried@gmail.com \
    /scratch/users/jgottf/CME/retrodict.jl
