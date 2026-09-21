# libcccrt Max/MSP externals

Author: Chris Kiefer

* `cccrpc~` — Random Projection Complexity, a port of the `CccRPC` SuperCollider UGen

The externals are built on the Eigen-free `core/` headers, so the only
dependency is the Max SDK.

### Requirements

- CMake >= 3.19
- [max-sdk-base](https://github.com/Cycling74/max-sdk-base) (or a full
  [max-sdk](https://github.com/Cycling74/max-sdk) checkout, which contains it at
  `source/max-sdk-base`)
- macOS: Xcode command line tools. Windows: Visual Studio.

### Cross-compiling for Windows from Linux/macOS

No Windows machine needed; MinGW-w64 links against the SDK's MSVC import
libraries fine, and the Max API is plain C so the different C runtime is
not a problem for these objects.

    sudo apt install g++-mingw-w64-x86-64     # or: brew install mingw-w64
    plugins/max/build-mingw.sh

This clones max-sdk-base beside the repository if needed and writes
`cccrpc~.mxe64` to `plugins/max/build-mingw/externals/` (or the directory
given as the first argument). The result still needs to be tested in Max on
Windows; the MSVC build via `build-windows.ps1` below is the reference.

### Building on Windows from scratch

`build-windows.ps1` installs Git, CMake and the Visual Studio 2022 Build Tools
(C++ workload) via winget, clones max-sdk-base, and builds the external. In an
elevated PowerShell:

    powershell -ExecutionPolicy Bypass -File plugins\max\build-windows.ps1 -Install

`-Install` copies the result into `Documents\Max 9\Library` (or Max 8). The
script can also be downloaded on its own and will clone this repository. Run
with `-SkipToolInstall` if the tools are already present.

### Building

From the external's directory:

    cd plugins/max/cccrpc~
    mkdir build && cd build
    cmake .. -DMAX_SDK_BASE_PATH=/path/to/max-sdk-base
    cmake --build . --config Release

`MAX_SDK_BASE_PATH` defaults to `../max-sdk-base` relative to this repository
(i.e. a sibling checkout), so it can be omitted if the SDK is cloned there.
On macOS pass `-DCMAKE_OSX_ARCHITECTURES="arm64;x86_64"` for a universal
binary; the SDK defaults to x86_64 only.

The external (`cccrpc~.mxo` / `cccrpc~.mxe64`) is written to
`plugins/max/externals/`. Copy it and `cccrpc~/cccrpc~.maxhelp` into a folder
in your Max search path (e.g. `~/Documents/Max 8/Library`).

### cccrpc~

    cccrpc~ [highDim] [lowDim] [maxWinSize] [maxLowDim]

| argument     | default | description                                                  |
|--------------|---------|--------------------------------------------------------------|
| `highDim`    | 10      | projection window length in samples (`h` in the paper)       |
| `lowDim`     | 2       | initial number of projection dimensions (`l`)                |
| `maxWinSize` | 500     | maximum analysis window in ms (sets buffer size)             |
| `maxLowDim`  | 8       | upper limit for `@lowdim` (raised to `lowDim` if larger)     |

`highDim`, `maxWinSize` and `maxLowDim` are fixed at creation.

| attribute  | default | description                                            |
|------------|---------|--------------------------------------------------------|
| `@lowdim`  | `lowDim` arg | number of projection dimensions, `l` (1 .. `maxLowDim`) |
| `@winsize` | 25      | analysis window in ms (clamped to `maxWinSize`)        |
| `@hopsize` | 0.5     | analysis hop as a fraction of `winsize`                |
| `@res`     | 5       | histogram resolution per dimension                     |
| `@rpchop`  | 0.5     | hop between projection windows, as a fraction of `highDim` |
| `@downsample` | 1    | average this many input samples into one before analysis   |

`@downsample N` reduces CPU use by a factor of N: the window still spans
`winsize` milliseconds but holds N times fewer points, so the analysis is N
times cheaper. Because the samples are averaged (a crude low-pass), RPC then
measures the complexity of the smoothed signal rather than the raw waveform.
Raising `@hopsize` is the other way to save CPU: it analyses less often
without changing what is measured.

In the notation of the paper (Kiefer 2023): `h` = `highDim`, `l` =
`@lowdim`, `alpha` = `@rpchop` × `highDim` samples, `beta` = `@res`.

`@lowdim` can change while running: the projection matrix is generated once
for `maxLowDim` rows and a projection into `l` dimensions uses its first `l`
rows, so no reallocation happens on the audio thread.

Signal in. Two outlets:

* left, signal — holds the latest complexity value (the number of occupied
  histogram cells), updated once per analysis hop
* right, float — the same value sent as a message after each hop (at most once
  per signal vector, so with hops shorter than the vector size only the last
  value in each vector is sent)

The projection matrix is generated with a fixed seed, so an object with the
same arguments always uses the same matrix. It is not the same matrix as the
SuperCollider/Python builds prior to the Eigen-free core (see
`core/README.md`).
