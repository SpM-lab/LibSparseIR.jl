# LibSparseIR.jl

[![Ask DeepWiki](https://deepwiki.com/badge.svg)](https://deepwiki.com/SpM-lab/LibSparseIR.jl)

Julia interface for libsparseir

> [!WARNING]
> This Julia project is still under construction. Please use other repositories:
> - https://github.com/SpM-lab/sparse-ir
> - https://github.com/SpM-lab/SparseIR.jl
> - https://github.com/SpM-lab/sparse-ir-fortran


## Set up

1. Install Julia
1. Clone libsparseir and build C-API

   ```sh
   git clone https://github.com/SpM-lab/libsparseir.git
   bash build_capi.sh
   ```
1. Clone this project(LibSparseIR.jl)
   ```sh
   git clone https://github.com/SpM-lab/LibSparseIR.git
   ```
1. Build this project:

   ```sh
   julia -e 'using Pkg; Pkg.build()'
   ```
1. It will create `src/C_API.jl` that wraps C-API of libsparseir for Julia interface.
1. To test our Julia project, run:
   ```sh
   julia -e 'using Pkg; Pkg.test()'
   ```
We use [ReTestItems.jl](https://github.com/JuliaTesting/ReTestItems.jl) as a test framework.
