# Companion code to the MSc thesis

Clone the repository and install Julia. Then `cd` to the directory and do
```
julia --proj
] activate .
instantiate
```
Opening a Jupyter session should then suffice to be able to run the notebook.

The notebook contains a pipeline to simulate demographic histories and analyse quantities
related to the statistics of Identity by Descent. The `Segments` object defined in the `Simulations.jl` module can be serialized to disk in an efficient binary format: a few simulations results are provided in this form.
