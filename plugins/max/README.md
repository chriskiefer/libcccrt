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

    cccrpc~ [highDim] [lowDim] [maxWinSize]

| argument     | default | description                                      |
|--------------|---------|--------------------------------------------------|
| `highDim`    | 10      | projection window length in samples              |
| `lowDim`     | 2       | number of projection dimensions                  |
| `maxWinSize` | 500     | maximum analysis window in ms (sets buffer size) |

These are fixed at creation, as in the SuperCollider UGen.

| attribute  | default | description                                            |
|------------|---------|--------------------------------------------------------|
| `@winsize` | 25      | analysis window in ms (clamped to `maxWinSize`)        |
| `@hopsize` | 0.5     | analysis hop as a fraction of `winsize`                |
| `@res`     | 5       | histogram resolution per dimension                     |
| `@rpchop`  | 0.5     | hop between projection windows, as a fraction of `highDim` |

Signal in, signal out. The output holds the latest complexity value (the
number of occupied histogram cells) and updates once per analysis hop; use
`snapshot~` to read it as a float.

The projection matrix is generated with a fixed seed, so an object with the
same arguments always uses the same matrix. It is not the same matrix as the
SuperCollider/Python builds prior to the Eigen-free core (see
`core/README.md`).
