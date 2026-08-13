ml julia

sbatch \
    --nodes=2 \
    --ntasks-per-node=2 \
    --mem-per-cpu=1G \
    --time=00:05:00 \
    --output=/scratch/users/jgottf/CME/output/%j.out \
    --error=/scratch/users/jgottf/CME/output/%j.out \
    --partition=normal,hns \
    --mail-type=ALL \
    --mail-user=juliansgottfried@gmail.com \
    /scratch/users/jgottf/CME/run2.jl
