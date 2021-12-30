# Examples

## Info

- ID's correspond to folders.

- All input data files needed for example calculations can be found in `examples/data`.

- Each example folder has a readme file with its description, `output`, and `analysis` folder. 

- To perform an example calculation, go to the main folder and run

        mpirun -np <n_proc> ./build/exdm ./examples/<example_ID>/input.txt

where `<n_proc>` is the number of processors to run on and `<example_ID>` is the example ID you want to compute for.

{!../examples/1/readme.md!}
