#!/usr/bin/python3

import sys
import os

from convert_to_edges import convert_file

if __name__ == "__main__":

	if len(sys.argv) != 3:
		print("usage: {0} <input directory> <output directory>".format(sys.argv[0]))
		exit()

	indir = sys.argv[1]
	outdir = sys.argv[2]


	for root, dirs, files in os.walk(indir, topdown=True):
		for dirname in dirs:
			dirpath = os.path.join(outdir, dirname)
			if not os.path.exists(dirpath):
				print("Making {0}".format(dirpath))
				os.makedirs(dirpath)
			else:
				print("{0} exists".format(dirpath))

		for fname in files:
			inpath = os.path.join(root, fname)
			if not inpath.endswith(".obj"):
				print("Skipping {0}".format(inpath))
				continue
			outpath = os.path.relpath(inpath, indir)
			outpath = os.path.join(outdir, outpath)
			outpath = outpath.replace(".obj", ".loop")
			convert_file(inpath, outpath)



