#!/usr/bin/python3

import sys

def convert_file(fname, outname):
	verts = []
	faces = []

	print("Converting {0} -> {1}".format(fname, outname))
	with open(fname, 'r') as objfile:
		for line in objfile:
			if line[0] == 'v' and line[1] == ' ':
				verts.append(line[1:].strip())
			elif line[0] == 'f' and line[1] == ' ':
				faces.append(line[1:].strip())

	with open(outname, 'w') as outfile:
		for line in verts:
			outfile.write("v {0}\n".format(line))
		for face in faces:
			indices = face.split(' ')
			count = len(indices)
			for i in range(count):
				i1 = indices[i].split("//")[0]
				i2 = indices[(i + 1) % count].split("//")[0]
				outfile.write("l {0} {1}\n".format(i1, i2))

if __name__ == "__main__":
	if len(sys.argv) < 2:
		print("usage: {0} <input OBJ> <optional output name>".format(sys.argv[0]))
		exit()

	fname = sys.argv[1]

	if not fname.endswith("obj"):
		print("Input must be a .obj file")
		exit()

	if len(sys.argv) > 2:
		outname = sys.argv[2]
	else:
		outname = fname.replace(".obj", ".loop")

	convert_file(fname, outname)

