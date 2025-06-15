using BinaryBuilder

name = "libsparseir"
version = v"0.3.1"

# Collection of sources required to complete build
sources = [
    GitSource("https://github.com/SpM-lab/libsparseir.git", "c04c1fe5de9f02012ffe654a0a157c36f3b921a3"),
]

# Bash recipe for building across all platforms
script = raw"""
cd ${WORKSPACE}/srcdir/libsparseir/bundle
./build.sh
cd dist/libsparseir-0.3.1
make
make install PREFIX=${prefix}
"""

# These are the platforms we will build for by default, unless further
# platforms are passed in on the command line
platforms = supported_platforms()
platforms = expand_cxxstring_abis(platforms)
filter!(p -> !(Sys.iswindows(p) && arch(p) == "i686"), platforms)
filter!(p -> !(Sys.islinux(p) && arch(p) == "powerpc64le"), platforms)

# The products that we will ensure are always built
products = [
    LibraryProduct("libsparseir", :libsparseir),
]

# Dependencies that must be installed before this package can be built
dependencies = [
    Dependency("CompilerSupportLibraries_jll"),
]

# Build the tarballs, and possibly a `build.jl` as well.
build_tarballs(ARGS, name, version, sources, script, platforms, products, dependencies;
    julia_compat="1.10", compilers=[:c, :cxx])
