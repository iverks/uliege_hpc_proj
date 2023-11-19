import typer
from pathlib import Path
from typing_extensions import Annotated
from typing import Optional
import numpy as np
from contextlib import ExitStack

import struct
app = typer.Typer()

@app.command()
def combine_datfiles(folder: Annotated[Optional[Path], typer.Option()] = Path("example_inputs/testing3d")):
    x_coord_min = float("inf")
    y_coord_min = float("inf")
    z_coord_min = float("inf")
    x_coord_max = -float("inf")
    y_coord_max = -float("inf")
    z_coord_max = -float("inf")
    filenames = {}
    for pathy in folder.iterdir():
        if not ("out_p" in pathy.name and pathy.name[-4:] == ".dat" ):
            continue
        # print(pathy)
        x, y, z, _out, _p = pathy.name.split("/")[-1].split("_")
        x, y ,z = int(x), int(y), int(z)
        x_coord_min = min(x, x_coord_min)
        y_coord_min = min(y, y_coord_min)
        z_coord_min = min(z, z_coord_min)
        x_coord_max = max(x, x_coord_max)
        y_coord_max = max(y, y_coord_max)
        z_coord_max = max(z, z_coord_max)
        filenames[(x, y, z)] = pathy

    with ExitStack() as stack, open("p_all.dat", "wb") as outfile:
        files = {index: stack.enter_context(open(fname, "rb")) for index, fname in filenames.items()}
        totnx, totny, totnz = 0,0,0
        totminx, totminy, totminz = float("inf"), float("inf"), float("inf")
        totmaxx, totmaxy, totmaxz = float("inf"), float("inf"), float("inf")

        x_ranges_set = set()
        y_ranges_set = set()
        z_ranges_set = set()
        # Read and combine headers 
        for idx, inp in files.items():
            nnx, = struct.unpack("i", inp.read(4))
            nny, = struct.unpack("i", inp.read(4))
            nnz, = struct.unpack("i", inp.read(4))
            
            minx, = struct.unpack("d", inp.read(8))
            maxx, = struct.unpack("d", inp.read(8))
            miny, = struct.unpack("d", inp.read(8))
            maxy, = struct.unpack("d", inp.read(8))
            minz, = struct.unpack("d", inp.read(8))
            maxz, = struct.unpack("d", inp.read(8))

            if minx not in x_ranges_set:
                totnx += nnx
                x_ranges_set.add(minx)
            if miny not in y_ranges_set:
                totny += nny
                y_ranges_set.add(miny)
            if minz not in z_ranges_set:
                totnz += nnz
                z_ranges_set.add(minz)

            totminx = min(totminx, minx)
            totminy = min(totminy, miny)
            totminz = min(totminz, minz)
            totmaxx = max(totmaxx, maxx)
            totmaxy = max(totmaxy, maxy)
            totmaxz = max(totmaxz, maxz)
        
        # Write header
        outfile.write(struct.pack("i", totnx))
        outfile.write(struct.pack("i", totny))
        outfile.write(struct.pack("i", totnz))
        outfile.write(struct.pack("d", totminx))
        outfile.write(struct.pack("d", totmaxx))
        outfile.write(struct.pack("d", totminy))
        outfile.write(struct.pack("d", totmaxy))
        outfile.write(struct.pack("d", totminz))
        outfile.write(struct.pack("d", totmaxz))

        # Read and combine arrays
        while True:
            arrays = {}
            for idx, inp in files.items():
                # Header per timestep
                ts, = struct.unpack("i", inp.read(4))
                readtime = inp.read(8)
                time, = struct.unpack("d", readtime)
                print(readtime)
                print(struct.pack("d", time))
                print(idx)
                mat: np.ndarray = np.fromfile(inp, dtype="double", count=nnx*nny*nnz)
                mat = mat.reshape((nnx, nny, nnz))
                arrays[idx] = mat
            
            array_stacks = {}
            for x in range(x_coord_min, x_coord_max + 1):
                for y in range(y_coord_min, y_coord_max + 1):
                    stacky = np.concatenate([arrays[(x, y, z)] for z in range(z_coord_min, z_coord_max + 1)], axis=2)
                    array_stacks[(x, y)] = stacky
            array_layers = {}
            for x in range(x_coord_min, x_coord_max + 1):
                layery = np.concatenate([array_stacks[(x,y)] for y in range(y_coord_min, y_coord_max + 1)], axis=1)
                array_layers[x] = layery
            
            totarr: np.ndarray = np.concatenate([array_layers[x] for x in range(x_coord_min, x_coord_max + 1)], axis=0)

            # Write timestep header
            outfile.write(struct.pack("i", ts))
            outfile.write(struct.pack("d", time))

            flattened: np.ndarray = totarr.flatten()
            for val in flattened:
                outfile.write(struct.pack("d", val))
            # flattened.tofile(outfile)
            print(flattened.shape, totnx, totny, totnz, totnx * totny * totnz) 

if __name__ == "__main__":
    app()