# Examples

## Info

- ID's correspond to folders.

- All input data files needed for example calculations can be found in `examples/data`.

- Each example folder has a readme file with its description, `output`, and `analysis` folder. 

- To perform an example calculation, go to the main folder and run

        mpirun -np <n_proc> ./build/exdm ./examples/<example_ID>/input.txt

where `<n_proc>` is the number of processors to run on and `<example_ID>` is the example ID you want to compute for.

- To run all the examples, from the main folder run

        ./examples/run_all_examples.sh <n_proc>

where `<n_proc>` is the number of processors to run on. WARNING: might take a while.

{!../examples/1/readme.md!}
{!../examples/2/readme.md!}
{!../examples/3/readme.md!}
{!../examples/4/readme.md!}
{!../examples/5/readme.md!}
{!../examples/6/readme.md!}
{!../examples/7/readme.md!}
{!../examples/8/readme.md!}
{!../examples/9/readme.md!}
{!../examples/10/readme.md!}
{!../examples/11/readme.md!}
{!../examples/12/readme.md!}
{!../examples/13/readme.md!}
{!../examples/14/readme.md!}
{!../examples/15/readme.md!}
