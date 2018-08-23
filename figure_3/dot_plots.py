#!/usr/bin/env python
import sys
import numpy as np
import matplotlib.pyplot as plt
from numpy import random
from shapely.geometry.polygon import LinearRing, Polygon
from descartes import PolygonPatch
from scipy.interpolate import spline

def generate_cells(n_species, n_functions):
	out=[]
	for sp in range(n_species):
		n = random.randint(2, 4)
		functions=[]
		while len(functions)<n:
			f = random.randint(n_functions)
			if f in functions: continue
			functions.append(f)
		out.append(functions)
	return out


def fill_gaps(organisms, n_functions):
	merge=[]
	for o in organisms:
		merge+=o
	missing=[]
	for n in range(n_functions):
		if n not in merge:
			missing.append(n)
	if len(missing)==1:
		missing.append(random.randint(n_functions))
	if missing!=[]:
		organisms.append(missing)
	return organisms


def generate_abundances(n):
	out=[]; a=1
	for i in range(n):
		#out.append(random.randint(1, max_abund))
		out.append(a)
		a+=1
	out.sort()
	return out


def generate_colors(n, new=False):
	if new==False:
		colors = ["blue", "green", "red", "cyan", "magenta", "gold", "seagreen", "salmon", "grey", "teal", "deeppink"]
	else:
		colors = ["crimson","slateblue","moccasin","y","steelblue","lightpink","brown","orange","olive","plum","wheat"]
	return colors[:n]


def modify_spines(ax):
	opacity=0.5
	ax.spines['left'].set_alpha(0)
	ax.spines['right'].set_alpha(0)
	ax.spines['top'].set_alpha(0)
	ax.spines['bottom'].set_alpha(opacity)


def draw_dots(ax, organisms, abundances, colors):
	y_coord=[]
	y=1
	for i,org in enumerate(organisms):
		x_list=[];y_list=[]
		st=y
		for n in range(abundances[i]):
			for x in org:
				x_list.append(x)
				y_list.append(y)
				totals[x]+=1
			y+=1
		fi=y
		ax.scatter(x_list, y_list, c=colors[i], linewidth=0.5, edgecolors='k', alpha=1)
		y_coord.append([st, fi])
	for xmaj in range(len(abundances)+2): ax.axvline(x=xmaj, ls='--', linewidth=1, alpha=0.1, c="k")
	ax.get_xaxis().set_ticks([])
	ax.get_yaxis().set_ticks([])
	ax.set_xlabel("Community functions", fontsize=12)
	ax.set_xlim([-1, n_functions])
	modify_spines(ax)
	return y_coord


def draw_landscape(ax, totals):
	x=range(len(totals))
	xnew = np.linspace(min(x),max(x),300)
	ynew = spline(x, totals, xnew)
	
	ax.plot(xnew, ynew)
	#ax.scatter(x, totals)
	for xmaj in range(len(totals)+1): ax.axvline(x=xmaj, ls='--', linewidth=1, alpha=0.1, c="k")
	ax.get_xaxis().set_visible(False)
	ax.get_yaxis().set_ticks([])
	ax.set_xlim([-1, n_functions])
	ax.set_ylim([-5, max(ynew)+5])
	modify_spines(ax)


def draw_connections(ax, before_coord, after_coord, colors):
	for i,before in enumerate(before_coord):
		after = after_coord[i]
		coordinates=[	(0,before[0]),
				(1,after[0]),
				(1,after[1]),
				(0,before[1])]
		trap = Polygon(coordinates)
		patch = PolygonPatch(trap, facecolor=colors[i], linewidth=0, alpha=1)
		ax.add_patch(patch)
	plt.axis('off')


def compute_distance(list1, list2):
	dist=0
	for i, a in enumerate(list1):
		b=list2[i]
		dist+=abs(a-b)
	return dist


def brute_force_rearrangement(totals, n_species):
	n=0; best_dist=1000
	while True:
		n+=1
		organisms = generate_cells(n_species, len(totals))
		abundances = generate_abundances(len(organisms))
		new_totals=[]
		for i in totals: new_totals.append(0)
		for i,org in enumerate(organisms):
			a=abundances[i]
			for f in org:
				new_totals[f]+=a

		if compute_distance(totals, new_totals)<best_dist:
			best_totals = new_totals[:]
			best_abundances = abundances[:]
			best_organisms = organisms[:]
			best_dist = compute_distance(totals, new_totals)
		if n%10000==0:
			if best_dist<20:
				break
			if n==100000:
				break
	return best_organisms, best_abundances, best_totals
		

def get_rearrangement_coordinates(new_organisms, new_abundances, max_y):
	rear_before=[]
	rear_after=[]
	y=0
	for i,org in enumerate(new_organisms):
		rear_before.append([1,1])
		rear_after.append([y, y+new_abundances[i]])
		y+=new_abundances[i]
	
	scaling = float(rear_after[-1][-1])/max_y
	for i in range(len(rear_after)):
		rear_after[i][0] = int(rear_after[i][0]/scaling+1)
		rear_after[i][1] = int(rear_after[i][1]/scaling+1)
	rear_after[-1][-1]-=1
	return rear_before, rear_after



# generate data
seed = random.randint(2**32-1)
seed = 435687341
random.seed(seed)
print("Seed was:", seed)

n_species=8
n_functions=8

organisms = generate_cells(n_species, n_functions)
organisms = fill_gaps(organisms, n_functions)
abundances = generate_abundances(len(organisms))
colors = generate_colors(len(organisms))
fig = plt.figure(figsize=(10, 10))
font = {'family': 'arial', 'weight': 'normal', 'size': 12}
plt.rc('font', **font)



# plot before
totals=[]
for i in range(n_functions): totals.append(0)
ax = fig.add_axes([0.05, 0.55, 0.35, 0.3])
ax.set_ylabel("Taxa")
before_coord = draw_dots(ax, organisms, abundances, colors)
ax = fig.add_axes([0.05, 0.85, 0.35, 0.1])
ax.set_ylabel("Func. landscape", fontsize=12)
draw_landscape(ax, totals)

ax = fig.add_axes([0.05, 0.05, 0.35, 0.3])
ax.set_ylabel("Taxa")
before_coord = draw_dots(ax, organisms, abundances, colors)
ax = fig.add_axes([0.05, 0.35, 0.35, 0.1])
ax.set_ylabel("Func. landscape", fontsize=12)
draw_landscape(ax, totals)




# plot gradual change
random.shuffle(abundances)
totals=[]
for i in range(n_functions): totals.append(0)
ax = fig.add_axes([0.6, 0.55, 0.35, 0.3])
after_coord = draw_dots(ax, organisms, abundances, colors)
ax = fig.add_axes([0.6, 0.85, 0.35, 0.1])
draw_landscape(ax, totals)


# plot joining lines
ax = fig.add_axes([0.4, 0.55, 0.2, 0.3])
draw_connections(ax, before_coord, after_coord, colors)
ax.set_ylim(0, sum(abundances)+2)



# plot rapid rearrangement
new_organisms, new_abundances, new_totals = brute_force_rearrangement(totals, n_species)
new_colors = generate_colors(len(new_organisms), new=True)
ax = fig.add_axes([0.6, 0.05, 0.35, 0.3])
after_coord = draw_dots(ax, new_organisms, new_abundances, new_colors)
ax = fig.add_axes([0.6, 0.35, 0.35, 0.1])
draw_landscape(ax, new_totals)


# plot rearangement lines
after_coord=[]
for i in before_coord:
	after_coord.append([before_coord[-1][-1], before_coord[-1][-1]])
rear_before, rear_after = get_rearrangement_coordinates(new_organisms, new_abundances, before_coord[-1][-1])
ax = fig.add_axes([0.4, 0.05, 0.2, 0.3])
draw_connections(ax, before_coord, after_coord, colors)
draw_connections(ax, rear_before, rear_after, new_colors)
ax.set_ylim(0, sum(abundances)+2)


plt.savefig("figure_3e.png", dpi=300)
plt.show()





