using BinaryBuilder

name = "libsparseir"
version = v"0.3.2"

# Collection of sources required to complete build
sources = [
    GitSource(
        "https://github.com/SpM-lab/libsparseir.git", 
        "8ceaf3b3f4843d6ff3451fc5106928c6977c2f68",
    ),
]

# Bash recipe for building across all platforms
script = raw"""
cd ${WORKSPACE}/srcdir/libsparseir/bundle
./build.sh
cd dist/libsparseir-0.3.2
make
make install PREFIX=${prefix}
"""

platforms = supported_platforms()
platforms = expand_cxxstring_abis(platforms)

products = [
    LibraryProduct("libsparseir", :libsparseir),
]

dependencies = [
    Dependency("CompilerSupportLibraries_jll"),
]

build_tarballs(ARGS, name, version, sources, script, platforms, products, dependencies;
    julia_compat="1.10", compilers=[:c, :cxx])
