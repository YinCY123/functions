#!/usr/bin/env python3

import argparse
import pandas as pd
import pyarrow.parquet as pq
from shapely import wkb
from shapely.geometry import Polygon, MultiPolygon

def decode_wkb_geometry(cell_id, binary_blob):
    try:
        raw = binary_blob.as_buffer().to_pybytes()
        geom = wkb.loads(raw)
        if geom.is_valid and not geom.is_empty:
            coords = []
            if isinstance(geom, Polygon):
                coords = list(geom.exterior.coords)
                for interior in geom.interiors:
                    coords.extend(interior.coords)
            elif isinstance(geom, MultiPolygon):
                for poly in geom.geoms:
                    coords.extend(poly.exterior.coords)
                    for interior in poly.interiors:
                        coords.extend(interior.coords)
            return [{"cell_id": cell_id, "x": x, "y": y} for x, y in coords]
    except Exception as e:
        print(f"Cell {cell_id}: decode failed → {e}")
    return []

def main():
    parser = argparse.ArgumentParser(description="Decode WKB polygons from a Parquet file")
    parser.add_argument("--input", required=True, help="Path to polygons.parquet")
    parser.add_argument("--output", required=True, help="Path to output CSV file")
    args = parser.parse_args()

    table = pq.read_table(args.input)
    geometry_column = table.column("geometry")

    rows = []
    for cell_id, cell in enumerate(geometry_column):
        rows.extend(decode_wkb_geometry(cell_id, cell))

    if rows:
        df = pd.DataFrame(rows)
        df.to_csv(args.output, index=False)
        print(f"✅ Exported {len(rows)} coordinate points to {args.output}")
    else:
        print("⚠️ No valid polygons found or decoded.")

if __name__ == "__main__":
    main()