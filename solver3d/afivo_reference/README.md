# Afivo-Streamer 3D Backend Reference

Stage D1 keeps Afivo-streamer as an external backend. The source is not vendored
into this repository.

- External path: `/home/helianthusczc/projects/afivo-streamer`
- Upstream: `https://github.com/MD-CWI/afivo-streamer.git`
- Pinned commit: `a50b5508775086e90dfe423455fb58d812578410`
- Internal Afivo core: source directory `afivo/` is part of the pinned tree; no
  independent Git submodule is present at this commit.
- License: GPL-3.0 (`LICENSE` and `afivo/LICENSE`)

Build command:

```sh
cd /home/helianthusczc/projects/afivo-streamer
make
```

On GCC 15, Hypre 2.31 may fail under the compiler default C23 mode. The minimal
observed compatibility repair is to rebuild the official local Hypre dependency:

```sh
cd /home/helianthusczc/projects/afivo-streamer/afivo/external_libraries
CFLAGS=-std=gnu17 ./build_hypre.sh
cd /home/helianthusczc/projects/afivo-streamer
make
```

Official test command:

```sh
cd /home/helianthusczc/projects/afivo-streamer
PATH=/home/helianthusczc/projects/streamer-rf-replica/.venv/bin:$PATH bash run_test.sh
```

True 3D smoke command:

```sh
cd /home/helianthusczc/projects/afivo-streamer/programs/standard_3d/tests
OMP_NUM_THREADS=1 ../streamer test_3d.cfg
```

Stage D2 common benchmark metadata and configs live in `common_benchmark/`.
