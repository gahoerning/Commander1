#!/usr/bin/env python3
"""Convert the supplied C-BASS linear-gain CSV into a Commander Stokes I bandpass."""
import argparse
import csv
import hashlib
import math
from pathlib import Path


def convert(source, destination):
    rows = []
    with source.open(newline='') as stream:
        for row in csv.reader(stream):
            if not row or not row[0].strip() or row[0].lstrip().startswith('#'):
                continue
            if len(row) != 13:
                raise ValueError('Expected frequency, six gains and six uncertainties (13 columns).')
            values = list(map(float, row))
            if not all(math.isfinite(value) for value in values):
                raise ValueError('All input values must be finite.')
            nu = values[0]
            if nu <= 0 or (rows and nu <= rows[-1][0]):
                raise ValueError('Frequencies must be positive and strictly increasing.')
            # The input is already linear. Follow the supplied plotting recipe:
            # take magnitudes, then average raw I1 (column 2) and I2 (column 7).
            rows.append((nu, (abs(values[1]) + abs(values[6])) / 2.0))
    if len(rows) < 2 or max(g for _, g in rows) <= 0:
        raise ValueError('At least two samples and nonzero intensity response are required.')
    peak = max(g for _, g in rows)
    digest = hashlib.sha256(source.read_bytes()).hexdigest()
    if source.resolve() == destination.resolve():
        raise ValueError('Input and output must be different files.')
    with destination.open('w') as stream:
        stream.write('# C-BASS North Stokes I, measured 2013-06-25; input gains are linear, not dB.\n')
        stream.write('# Input SHA256: ' + digest + '\n')
        stream.write('# Response = (abs(I1)+abs(I2))/2, divided by its maximum.\n')
        stream.write('# No spectral or calibration weighting has been applied here.\n')
        stream.write('# Frequency [GHz]  Linear relative response\n')
        for nu, gain in rows:
            stream.write(f'{nu:.12e} {gain/peak:.16e}\n')
    return len(rows)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('source', type=Path)
    parser.add_argument('destination', type=Path)
    args = parser.parse_args()
    print(f'Wrote {convert(args.source, args.destination)} Stokes I samples to {args.destination}')
