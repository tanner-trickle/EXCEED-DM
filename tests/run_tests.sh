#!/usr/bin/env bash

echo "Beginning tests..."

FoBiS.py clean
FoBiS.py build -mode local-gnu

files=(\
    "Si/inputs/vc_scatter_minimal.txt"\
    "Si/inputs/vc_scatter_example.txt"\
    "Si/inputs/vc_scatter_spin.txt"\
    "Si/inputs/vc_scatter_spin_TFF.txt"\
    "Si/inputs/cc_scatter_example.txt"\
    "Si/inputs/cf_scatter_example.txt"\
    "Si/inputs/vf_scatter_example.txt"\
    "Si/inputs/ps_absorp_example.txt"\
    "Si/inputs/vector_absorp_example.txt"\
    "Si/inputs/scalar_absorp_example.txt"\
    "Si/inputs/scalar_absorp_spin_example.txt"\
)

# files=(\
#     "Si/inputs/ps_absorp_example.txt"\
#     "Si/inputs/vector_absorp_example.txt"\
#     "Si/inputs/scalar_absorp_example.txt"\
# )

for file in "${files[@]}"; do
    echo ""
    echo "Testing ./tests/${file}..."
    echo ""
    mpirun -np 12 --use-hwthread-cpus ./build/exdm "./tests/${file}"
done

echo "Testing complete!"
echo ""
echo "Time to complete : "
