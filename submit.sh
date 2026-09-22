ml julia

sbatch \
    --ntasks=114 \
    --cpus-per-task=1 \
    --mem-per-cpu=16G \
    --time=00:10:00 \
    --output=/scratch/users/jgottf/CME/output/%j.out \
    --error=/scratch/users/jgottf/CME/output/%j.out \
    --partition=normal,hns \
    --mail-type=ALL \
    --mail-user=juliansgottfried@gmail.com \
    /scratch/users/jgottf/CME/convert.jl
