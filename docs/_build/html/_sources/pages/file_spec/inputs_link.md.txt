```{toctree}
   :maxdepth: 2
   :collapse: False
```

# Input File

```{contents} Table of Contents
   :depth: 3
   :local:
   :backlinks: none
```

## Introduction

The **input file** is simply a text file containing runtime variables. Variables are collected in to *groups* of similar variables. The group name is inside brackets, [], and the corresponding variables are below it. For example, the **control** group might contain,

    [control]
        calculation     = 'binned_scatter_rate'
        run_description = 'example_1'
        out_folder      = './examples/1/output/'

The description of all the groups, along with their their variables and default values, is given below.

While most variables take a singular value, some accept an array. The dimensions of the input array must match that specified by **Dim** in the variable description. For example, **a_vecs_Ang** in the **material** group has **Dim: [3, 3]**, and therefore the input must be specified as a $3 \times 3$ matrix. This can be done in two ways: placing the rows side by side, e.g.,

    a_vecs_Ang = 1, 0, 0, 0, 1, 0, 0, 0, 1

for the identity matrix, or with the **+=** operator (the preferred method for matrix inputs for readability),

    a_vecs_Ang  = 1, 0, 0
    a_vecs_Ang += 0, 1, 0
    a_vecs_Ang += 0, 0, 1

A colon in **Dim** signifies that the variable accepts an array of any size. For example, **mX** in the **dm_model** group has **Dim: [:]**, meaning that a list of variables is accepted,

    mX = 1e6, 1e7, 1e8, ...

A mixed dimension, e.g., **Dim: [:, 3]** means that any number of length 3 vectors may be specified. For example, **v_e_km_per_sec** in the **astroph_model** group may be specified as

    v_e_km_per_sec  = 0, 0, 230
    v_e_km_per_sec += 0, 0, 240
    v_e_km_per_sec += 0, 0, 250
    $\vdots$

See <https://github.com/jannisteunissen/config_fortran> for more details concerning the input format.

## Groups

```{include} inputs.md
```
