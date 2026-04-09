import pyvista as pv
import pykep as pk
import numpy as np
from pyvista.core.utilities import lines_from_points
import spiceypy
import matplotlib.pyplot as plt

pk.utils.load_spice_kernels("../data/de440.bsp")
spiceypy.furnsh("../data/earth_1962_250826_2125_combined.bpc")

emb = pk.udpla.spice("earth", "ECLIPJ2000", "ssb")
moon = pk.udpla.spice("moon", "ECLIPJ2000", "earth")

trajfile = "../res/output.txt"
depotfile = "../res/output2.txt"

words = []
wordsdepot = []
with open(trajfile, "r") as f:
	for line in f:
		words.append(line.split())

with open(depotfile, "r") as f:
	for line in f:
		wordsdepot.append(line.split())

del words[0]
del wordsdepot[0]

words = list(map(lambda x:list(map(float, x)), words))
wordsdepot = list(map(lambda x:list(map(float, x)), wordsdepot))

times = []
xs = []
ys = []
zs = []

timesdepot = []
xsdepot = []
ysdepot = []
zsdepot = []

moonpts = []

pmags = []

print("Reading...")
for word in words:
	# pykep undoes this in its code lmao
	epos = emb.eph((word[0] + 43135.816087188054)/86400)[0]
	mpos = moon.eph((word[0] + 43135.816087188054)/86400)[0]

	times.append(word[0])
	xs.append(word[1] - 0.001 * epos[0])
	ys.append(word[2] - 0.001 * epos[1])
	zs.append(word[3] - 0.001 * epos[2])

	moonpts.append([0.001 * x for x in mpos])

	pmags.append(np.sqrt(word[7] * word[7] + word[8] * word[8] + word[9] * word[9]))

for word in wordsdepot:
	# pykep undoes this in its code lmao
	epos = emb.eph((word[0] + 43135.816087188054)/86400)[0]
	mpos = moon.eph((word[0] + 43135.816087188054)/86400)[0]

	timesdepot.append(word[0])
	xsdepot.append(word[1] - 0.001 * epos[0])
	ysdepot.append(word[2] - 0.001 * epos[1])
	zsdepot.append(word[3] - 0.001 * epos[2])

for time, pmag in zip(times, pmags):
	print(time, pmag)

plt.plot(times, pmags)
plt.show()

rotate = spiceypy.pxform("ITRF93", "ECLIPJ2000", times[0])
pole = np.array([rotate[0][2], rotate[1][2], rotate[2][2]])

sphere = pv.Sphere(radius=6371)

pl = pv.Plotter()
pl.add_mesh(sphere, color="blue", show_edges=False)

pts = np.dstack((xs, ys, zs))[0]
line = lines_from_points(pts)

pl.add_mesh(line, show_edges=True, scalars=np.array(times), clim=[min(times), max(times)], cmap="turbo")

polepts = np.array([pole * 10000, pole * -10000])
poleline = lines_from_points(polepts)

pl.add_mesh(poleline, show_edges=True, color="black")

moonline = lines_from_points(np.array(moonpts))
pl.add_mesh(moonline, show_edges=True, scalars=np.array(times), clim=[min(times),max(times)], cmap="turbo")

ptsdepot = np.dstack((xsdepot, ysdepot, zsdepot))[0]
linedepot = lines_from_points(ptsdepot)

pl.add_mesh(linedepot, show_edges=True, scalars=np.array(timesdepot), clim=[min(timesdepot), max(timesdepot)], cmap="turbo")

pl.show()