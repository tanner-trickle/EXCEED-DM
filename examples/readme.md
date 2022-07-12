# Examples

The `examples` folder contains a pleathora of self-contained examples to get started with running (and processing the results of) `EXCEED-DM`. Each folder contains its own: 

- `readme.md`: explaining in words what calculation is doing. Repeated below for convenience.

- Input file: `input.txt`. This is where the run parameters (dark matter masses, velocity distribution, etc.) are defined. 

- Electronic configuration file: `elec_config.hdf5`. This contains all the information about the electronic states used in the calculation. For example, in example 1 only the valence and conduction bands in Si are included. Other examples will include more, or other, states. 

- Analysis Jupyter notebook: `analysis.ipynb`, to study the output in more detail.

- `output` folder, where the `EXCEED-DM` calculation output is stored.

## Running

- To perform an example calculation, go to the main folder and run

        mpirun -np <n_proc> ./build/exdm ./examples/<example_ID>/input.txt

where `<n_proc>` is the number of processors to run on and `<example_ID>` is the example ID you want to compute for.

- To run all the examples, from the main folder run

        ./examples/run_all_examples.sh <n_proc>

where `<n_proc>` is the number of processors to run on. 

```{warning}
Running all of the examples can take some time.
```

---

```{include} ../../examples/1/readme.md
```
```{include} ../../examples/2/readme.md
```
```{include} ../../examples/3/readme.md
```
```{include} ../../examples/4/readme.md
```
```{include} ../../examples/5/readme.md
```
```{include} ../../examples/6/readme.md
```
```{include} ../../examples/7/readme.md
```
```{include} ../../examples/8/readme.md
```
```{include} ../../examples/9/readme.md
```
```{include} ../../examples/10/readme.md
```
```{include} ../../examples/11/readme.md
```
```{include} ../../examples/12/readme.md
```
```{include} ../../examples/13/readme.md
```
