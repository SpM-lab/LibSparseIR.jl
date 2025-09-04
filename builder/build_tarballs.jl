using BinaryBuilder, Pkg

name = "libsparseir"
version = v"0.4.2"

# Collection of sources required to complete build
sources = [
    GitSource(
        "https://github.com/SpM-lab/libsparseir.git", 
        "bb5147da806c82e82695da7701b9421182105765",
    ),
    GitSource(
        "https://github.com/tuwien-cms/libxprec.git",
        "d35f3fa9a962d3f96a1eef63132030fd869c183a"
    )
]

# Bash recipe for building across all platforms
script = raw"""
cd ${WORKSPACE}/srcdir/libsparseir/
install_license LICENSE
if [[ "${target}" == *mingw* ]]; then
  LBT="${libdir}/libblastrampoline-5.dll"
else
  LBT="-lblastrampoline"
fi

c++ -O2 -Wall -fPIC -shared -std=c++11 -I${includedir}/eigen3/ -Iinclude -I../libxprec/include ${LBT} src/*.cpp -o ${prefix}/lib/libsparseir.${dlext}
mkdir -p ${prefix}/include/sparseir/
cp include/sparseir/sparseir.h include/sparseir/spir_status.h include/sparseir/version.h ${prefix}/include/sparseir/
"""

platforms = [
        # glibc Linuces
        Platform("i686", "linux"),
        Platform("x86_64", "linux"),
        Platform("aarch64", "linux"),
        Platform("armv6l", "linux"),
        Platform("armv7l", "linux"),
        # Platform("powerpc64le", "linux"), # does not work
        # Platform("riscv64", "linux"), # does not work

        # musl Linuces
        # Platform("i686", "linux"; libc="musl"), # does not work
        Platform("x86_64", "linux"; libc="musl"),
        Platform("aarch64", "linux"; libc="musl"),
        Platform("armv6l", "linux"; libc="musl"),
        Platform("armv7l", "linux"; libc="musl"),

        # BSDs
        Platform("x86_64", "macos"),
        Platform("aarch64", "macos"),
        Platform("x86_64", "freebsd"),
        Platform("aarch64", "freebsd"),

        # Windows
        Platform("i686", "windows"),
        Platform("x86_64", "windows"),
    ]

platforms = [Platform("aarch64", "macos"),]

platforms = expand_cxxstring_abis(platforms)

products = [
    LibraryProduct("libsparseir", :libsparseir),
]

dependencies = [
    Dependency("CompilerSupportLibraries_jll"),
    BuildDependency("Eigen_jll"),
    Dependency(PackageSpec(name="libblastrampoline_jll", uuid="8e850b90-86db-534c-a0d3-1478176c7d93"), compat="5.4.0"),
]

build_tarballs(ARGS, name, version, sources, script, platforms, products, dependencies;
    julia_compat="1.10", compilers=[:c, :cxx])
