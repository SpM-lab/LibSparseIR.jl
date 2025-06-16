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

platforms = expand_cxxstring_abis(platforms)

products = [
    LibraryProduct("libsparseir", :libsparseir),
]

dependencies = [
    Dependency("CompilerSupportLibraries_jll"),
]

build_tarballs(ARGS, name, version, sources, script, platforms, products, dependencies;
    julia_compat="1.10", compilers=[:c, :cxx])
