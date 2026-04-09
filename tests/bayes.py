import matplotlib.pyplot as plt
import numpy as np
from scipy.special import erf

X = [
	9,
	8,
	8,
	8,
	7,
	6,
	8,
	4,
	5,
	0,
	5,
	6,
	8.5,
	8.5,
	7,
	4,
	8.288,
	3.105,
	7.5,
	8.944,
	8.986,
	7.2,
	7.255,
	8.4,
	8.978,
	-2.849,
	1.311,
	-1.73,
	7.4,
	6.748,
	3.353,
	2.968,
	-0.678,
	1.398, -3.474
]
Y = [
	0,
	0,
	1,
	4,
	5,
	6,
	-4,
	8,
	7,
	-9,
	-7,
	-6,
	-2.5,
	2.5,
	-5.5,
	-8,
	3.508,
	-8.447,
	-4.97,
	-1,
	-0.5,
	-4.9,
	-5.325,
	-1.296,
	0.629,
	-6.775,
	-7.232,
	-8.791,
	-2.987,
	2.914,
	-6.54,
	6.723,
	-7.318,
	-8.85, -8.259
]
S = [
	63.770,
	4.070,
	5.629,
	27.260,
	11.346,
	9.356,
	50.139,
	20.850,
	5.839,
	29.944,
	21.274,
	24.066,
	30.963,
	27.213,
	55.854,
	33.491,
	32.788,
	30.905,
	63.140,
	63.387,
	65.001,
	28.487,
	67.547,
	15.612,
	60.079,
	0, 0,
	15.1,
	1.25,
	0,
	0, 0, 0, 0, 0
]
avg = sum(S) / len(S)
S = np.array(S)

S -= avg

def l(x):
	return np.sqrt((x[0]) ** 2 + (x[1]) ** 2)

def kernel(x, y):
	out = []
	for i in range(len(x)):
		row = []
		for j in range(len(y)):
			row.append(200 * np.exp(-1 * l(x[i] - y[j]) ** 2))
		out.append(row)
	return np.array(out)

def EI(mean, stddev):
	return stddev * np.exp(-mean * mean / 2 / stddev / stddev) / np.sqrt(2 * np.pi) + mean / 2 - mean * erf(-mean / stddev / np.sqrt(2)) / 2

testxs = []
for r in np.arange(7, 9.1, 0.07):
	if r > 9:
		continue
	for theta in np.arange(0, 6.3, 0.0113234):
		x = r * np.cos(theta)
		y = r * np.sin(theta)
		if x * x + y * y < 53.1467247343:
			continue

		dot = x * 0.917477 - y * 0.397789

		if dot < 0:
			continue

		testxs.append([x, y])
testxs = np.array(testxs)

pts = np.dstack((X, Y))[0]

uncertainty = 0

inverse = np.linalg.inv(kernel(pts, pts) + uncertainty * uncertainty * np.eye(len(pts)))
mean = kernel(testxs, pts) @ inverse @ S
cov = kernel(testxs, testxs) - kernel(testxs, pts) @ inverse @ kernel(pts, testxs)

S += avg
mean += avg

ei = EI(mean - np.max(S), np.sqrt(np.diagonal(cov)))

b = 0
bp = []
for i in range(len(ei)):
	print(testxs[i], ei[i])

	if ei[i] > b:
		b = ei[i]
		bp = testxs[i]
print(bp, b)

ax = plt.figure().add_subplot(projection='3d')

#ax.scatter(X, Y, S)
ax.scatter(testxs[:,0], testxs[:,1], ei)

ax.set_xlabel("X")
ax.set_ylabel("Y")
ax.set_zlabel("Z")

plt.show()