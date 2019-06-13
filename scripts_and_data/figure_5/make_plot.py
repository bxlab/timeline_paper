#!/usr/bin/env python
import sys
import numpy as np
import matplotlib.pyplot as plt
from numpy import random
from shapely.geometry.polygon import LinearRing, Polygon
from descartes import PolygonPatch
from scipy.interpolate import spline

def generate_cells(n_species, n_functions, silent=False):
	if silent==False:
		print "generating community..."
	out=[]
	for sp in range(n_species):
		n = random.randint(2, 4)
		functions=[]
		while len(functions)<n:
			f = random.randint(n_functions)
			if f in functions: continue
			functions.append(f)
		out.append(functions)
	all_f=[]
	for sp in out:
		all_f += sp
	for f in range(1, n_functions):
		if f not in all_f:
			out = generate_cells(n_species, n_functions, silent=True)
			break
	return out


def fill_gaps(organisms, n_functions, silent=False):
	if silent==False:
		print "adding species to make sure each function is represented..."
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


def generate_abundances(n, silent=False):
	if silent==False:
		print "generating abundances..."
	out=[]; a=1
	for i in range(n):
		#out.append(random.randint(1, max_abund))
		out.append(a)
		a+=1
	out.sort()
	return out


def generate_colors(n, new=False):
	print "generating colors..."
	if new==False:
		colors = ["firebrick", "purple", "red", "teal", "magenta", "blue", "seagreen", "salmon", "grey", "teal", "deeppink"]
	else:
		colors = ["crimson","slateblue","orchid","y","steelblue","olive","brown","orange","olive","plum","wheat"]
	return colors[:n]


def modify_spines(ax):
	opacity=0.5
	ax.spines['left'].set_alpha(0)
	ax.spines['right'].set_alpha(0)
	ax.spines['top'].set_alpha(0)
	ax.spines['bottom'].set_alpha(opacity)


def draw_dots(ax, orig_organisms, orig_abundances, orig_colors, n, labelside, seedbank):
	if seedbank=="bot":
		organisms = [range(n)] + orig_organisms
		abundances = [1] + orig_abundances
		colors = ["k"] + orig_colors
		ax.set_ylim(0, sum(abundances)+2)
	if seedbank=="top":
		organisms = orig_organisms + [range(n)]
		abundances = orig_abundances + [1]
		colors = orig_colors + ["k"]
		ax.set_ylim(0, sum(abundances)+2)

	totals=[]
	for i in range(n): totals.append(0)

	draw_boxes(ax, abundances, colors, n, labelside, seedbank)
	print "drawing 'dot' plot of community..."
	y_coord=[]
	y=1
	for i,org in enumerate(organisms):
		st=y
		x_list=[];y_list=[]
		for n in range(abundances[i]):
			for x in org:
				x_list.append(x)
				y_list.append(y)
				totals[x]+=1
			y+=1
		fi=y
		#ax.scatter(x_list, y_list, c=colors[i], linewidth=0.5, edgecolors='k', alpha=1, zorder=10)
		width=0.5
		for x in org:
			coordinates=[ (x-0.5,st-0.5), (x+0.5,st-0.5), (x+0.5,fi-0.5), (x-0.5,fi-0.5) ]
			trap = Polygon(coordinates)
			patch = PolygonPatch(trap, facecolor=colors[i], linewidth=0.5, alpha=1, zorder=1, edgecolor="white")
			ax.add_patch(patch)
		if i!=0:
			y_coord.append([st, fi])
	#for xmaj in range(len(abundances)+1): ax.axvline(x=xmaj-0.5, ls='--', linewidth=1, alpha=0.1, c="k")
	ax.get_xaxis().set_ticks([])
	ax.get_yaxis().set_ticks([])
	ax.set_xlabel("Functional categories", fontsize=12)
	ax.set_xlim([-1, n_functions])
	modify_spines(ax)
	return y_coord, totals


def draw_boxes(ax, abundances, colors, n, labelside, seedbank):
	print "drawing background boxes and labels..."
	bot=0.5

	if seedbank=="top":
		labels=["Seedbank","Taxon A","Taxon B","Taxon C","Taxon D","Taxon E","Taxon F","Taxon G","Taxon H"]
	if seedbank=="bot":
		labels=["Taxon I","Taxon J","Taxon K","Taxon L","Taxon M","Taxon N","Taxon A","Taxon B", "Seedbank"]

	for i, abund in enumerate(abundances):
		top=bot+abund
		coordinates=[ (-1,top), (n,top), (n,bot), (-1,bot) ]
		trap = Polygon(coordinates)
		patch = PolygonPatch(trap, facecolor=colors[i], linewidth=1, edgecolor="white", alpha=0.5, zorder=1)
		ax.add_patch(patch)


		y=(top+bot)/2-0.5
		taxa=labels[len(labels)-i-1]
		if taxa=="Seedbank" and seedbank=="bot":
			y-=0.5

		if labelside=="left":
			x=-1.1
			ha="right"
		if labelside=="right":
			x=n+0.1
			ha="left"
		ax.text(x, y, taxa, ha=ha, fontsize=10)
		bot=top


def draw_landscape(ax, totals, label):
	print "drawing functional landscape of community..."
	totals = [totals[0]] + totals + [totals[-1]]

	x=range(len(totals))
	xnew = np.linspace(min(x),max(x),300)
	ynew = spline(x, totals, xnew)
	
	for i,x in enumerate(xnew):
		y=ynew[i]
		if y<0:
			return "bad"
		ax.plot([x,x], [0,y], '-', c="lightgray", alpha=1)
	ax.plot(xnew, ynew, c="k", linewidth=0.5)
	ax.scatter(range(1,len(totals[1:-1])+1), totals[1:-1], facecolors='white', edgecolors='k', zorder=10)

	#for xmaj in range(len(totals)+1): ax.axvline(x=xmaj-0.5, ls='--', linewidth=1, alpha=0.1, c="k")
	ax.get_xaxis().set_visible(False)
	ax.get_yaxis().set_ticks([])
	plt.axis("off")
	ax.set_xlim([0, n_functions+1])
	ax.set_ylim([-5, max(ynew)+5])
	modify_spines(ax)
	ax.set_title(label, fontsize=14)
	ttl = ax.title
	ttl.set_position([.5, 1.0])


def draw_connections(ax, before_coord, after_coord, colors):
	print "drawing lines connecting showing abundance changes..."
	for i,before in enumerate(before_coord):
		after = after_coord[i]
		coordinates=[	(0,before[0]-0.5),
				(1,after[0]-0.5),
				(1,after[1]-0.5),
				(0,before[1]-0.5)]
		trap = Polygon(coordinates)
		patch = PolygonPatch(trap, facecolor=colors[i], linewidth=1, edgecolor="white", alpha=0.5)
		ax.add_patch(patch)
	plt.axis('off')


def compute_distance(list1, list2):
	dist=0
	for i, a in enumerate(list1):
		b=list2[i]
		dist+=abs(a-b)
	return dist


def brute_force_rearrangement(totals, n_species, old_organisms, distance_max):
	print "generating random community of equivalent functional potential..."
	n=0; best_dist=1000
	while True:
		n+=1
		organisms = generate_cells(n_species-2, len(totals), silent=True) + old_organisms[:2]
		abundances = generate_abundances(len(organisms), silent=True)
		abundances[-3] += abundances[-1]-3
		abundances[-4] += abundances[-2]-4
		abundances[-1] -= abundances[-1]-3
		abundances[-2] -= abundances[-2]-4
		
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
			if best_dist<distance_max:
				break
		if n==500000:
			break
	return best_organisms, best_abundances
		

def get_rearrangement_coordinates(new_organisms, new_abundances, max_y):
	print "calculating coordinates for lines representing rearrangement adundance changes..."
	rear_before=[]
	rear_after=[]
	y=0
	for i,org in enumerate(new_organisms):
		rear_before.append([float(i)/len(new_organisms)+1,float(i+1)/len(new_organisms)+1])
		rear_after.append([y, y+new_abundances[i]])
		y+=new_abundances[i]
	
	scaling = float(rear_after[-1][-1])/(max_y-1)
	scaling=1
	for i in range(len(rear_after)):
		rear_after[i][0] = int(rear_after[i][0]/scaling+1)
		rear_after[i][1] = int(rear_after[i][1]/scaling+1)
	#rear_after[-1][-1]-=1
	return rear_before, rear_after


def draw_label(ax, label):
	ax.text(0.5, 0.5, label, ha="center", fontsize=24)
	ax.axis("off")

def draw_grid(ax):
	ax.set_xlim([0, n_functions+1])
	for xmaj in range(n_functions+2): ax.axvline(x=xmaj-0.5, ls='--', linewidth=1, alpha=0.12, c="k", zorder=-1)
	#ax.get_xaxis().set_visible(False)
	#ax.get_yaxis().set_ticks([])
	plt.axis("off")


def calculate_ri (organisms, before, after):
	functions={}
	abundances={}
	for i,a in enumerate(before):
		b=after[i]
		org_name="spp_"+str(i)
		abundances[org_name]=[a,b]

	for o,org in enumerate(organisms):
		org_name="spp_"+str(o)
		for i in org:
			function="func_"+str(i)
			if function not in functions:
				functions[function] = {}
			functions[function][org_name] = []

	RIs=[]
	for func in functions:
		numerator=0
		denominator=0
		for org in functions[func]:
			a = abundances[org][0]
			b = abundances[org][1]
			numerator+=abs(b-a)
			denominator+=(a+b)
		RI = numerator*1.0/denominator
		RIs.append(RI)
	print "\nGradual adjustement RI="+str(np.mean(RIs))[:4] + "+/-" + str(np.std(RIs))[:4] + "\n"
	return np.mean(RIs), np.std(RIs)


def calculate_rear_ri (before_organisms, before_abundances, after_organisms, after_abundances):
	functions={}
	abundances={}
	for i,a in enumerate(before_abundances):
		org_name="spp_before_"+str(i)
		abundances[org_name]=[a,0]
	for i,b in enumerate(after_abundances):
		if i==0 or i==1:
			org_name="spp_before_"+str(i)
			abundances[org_name][1] = b
		else:
			org_name="spp_after_"+str(i)
			abundances[org_name]=[0,b]
	for o,org in enumerate(before_organisms):
		org_name="spp_before_"+str(o)
		for i in org:
			function="func_"+str(i)
			if function not in functions:
				functions[function] = {}
			functions[function][org_name] = []
	for o,org in enumerate(after_organisms):
		if o==0 or o==1:
			org_name="spp_before_"+str(o)
			continue
		org_name="spp_after_"+str(o)
		for i in org:
			function="func_"+str(i)
			if function not in functions:
				functions[function] = {}
			functions[function][org_name] = []
	
	RIs=[]
	for func in functions:
		numerator=0
		denominator=0
		for org in functions[func]:
			a = abundances[org][0]
			b = abundances[org][1]	
			numerator+=abs(b-a)
			denominator+=(a+b)
		RI = numerator*1.0/denominator
		RIs.append(RI)
	print "\nRapid rearrangement RI="+str(np.mean(RIs))[:4] + "+/-" + str(np.std(RIs))[:4] + "\n"
	return np.mean(RIs), np.std(RIs)
			

print "\nGENERATING DATA"

# random generation seed
seed = None
seed = 4013758737

#core parameters
n_species=8; n_functions=8
fig = plt.figure(figsize=(10, 10))
font = {'family': 'arial', 'weight': 'normal', 'size': 12}
plt.rc('font', **font)


left=0.1
middle=0.4
right=0.6
width=0.3
height=0.3
landscape_height=0.08
#top=0.53
#bot=0.05
bot=0.53
top=0.05


print '\nPLOTTING THE "BEFORE" STATES'
while True:
	if seed==None:
		s = random.randint(2**32-1)
	else:
		s=seed
	random.seed(s)
	plt.clf()
	organisms = generate_cells(n_species, n_functions)
	n_species=len(organisms)
	abundances = generate_abundances(len(organisms))
	colors = generate_colors(len(organisms))
	
	ax = fig.add_axes([left, top, width, height])
	before_coord, before_totals = draw_dots(ax, organisms, abundances, colors, n_functions, labelside="left", seedbank="bot")

	ax = fig.add_axes([left, top+height*0.95, width, landscape_height])
	state = draw_landscape(ax, before_totals, "Transient functional landscape")

	ax = fig.add_axes([right, bot, width, height])
	before_coord, before_totals = draw_dots(ax, organisms, abundances, colors, n_functions, labelside="right", seedbank="bot")
	ax = fig.add_axes([right, bot+height*0.95, width, landscape_height])
	state = draw_landscape(ax, before_totals, "Transient functional landscape")


	if state=="bad":
		print "it would have been an ugly plot anyways... retry\n"
		continue
	
	transient_abundances=abundances[:]
	random.shuffle(abundances)
	ax = fig.add_axes([right, top, width, height])
	after_coord, after_totals = draw_dots(ax, organisms, abundances, colors, n_functions, labelside="right", seedbank="bot")
	ax = fig.add_axes([right, top+height*0.95, width, landscape_height])
	state = draw_landscape(ax, after_totals, "Recovered functional landscape")

	change = compute_distance(before_totals, after_totals)
	if change<40:
		print "the functional shift was boring (diff="+str(change)+")... retry\n"
		continue
	if state=="bad":
		print "it would have been an ugly plot anyways... retry\n"
		continue
	
	ax = fig.add_axes([middle, top, 0.2, height])
	draw_connections(ax, before_coord, after_coord, colors)
	ax.set_ylim(0, sum(abundances)+3)


	print "\nSattisfactory solution found! (change="+str(change)+")"
	print("Random seed was:", s)

	adjustment_ri, adjustment_ri_std = calculate_ri (organisms, transient_abundances, abundances)
	break



print "\nPLOTTING RAPID REARRANGEMENT"
new_organisms, new_abundances = brute_force_rearrangement(after_totals, n_species, organisms, 20)
new_colors = generate_colors(len(new_organisms)-2, new=True) + colors[:2]
rearrangement_ri, rearrangement_ri_std = calculate_rear_ri (new_organisms, new_abundances, organisms, transient_abundances)

ax = fig.add_axes([left, bot, width, height])
after_coord, new_totals = draw_dots(ax, new_organisms, new_abundances, new_colors, n_functions, labelside="left", seedbank="top")
ax = fig.add_axes([left, bot+height*0.95, width, landscape_height])
state = draw_landscape(ax, new_totals, "Initial functional landscape")


# plot joining lines
rear_before, rear_after = get_rearrangement_coordinates(new_organisms, new_abundances, before_coord[-1][-1])

after_coord=[]
for i in range(len(before_coord)):
	if i<2:
		after_coord.append(rear_after[-2+i])
	else:
		after_coord.append([before_coord[-1][-1]+float(i)/len(before_coord)-1, before_coord[-1][-1]+float(i+1)/len(before_coord)-1])

rear_before=rear_before[:-2]
rear_after=rear_after[:-2]
	
ax = fig.add_axes([middle, bot, 0.2, height])
draw_connections(ax, after_coord, before_coord, colors)
draw_connections(ax, rear_after, rear_before, new_colors)
ax.set_ylim(0, sum(abundances)+3)


# plot backgorund grid lines
ax = fig.add_axes([left, top+height*0.95, width, landscape_height])
draw_grid(ax)
ax = fig.add_axes([right, top+height*0.95, width, landscape_height])
draw_grid(ax)
ax = fig.add_axes([left, bot+height*0.95, width, landscape_height])
draw_grid(ax)
ax = fig.add_axes([right, bot+height*0.95, width, landscape_height])
draw_grid(ax)


# draw panel labels
ax = fig.add_axes([0, top+height+landscape_height, 0.05, 0.05])
draw_label(ax, "B")
ax = fig.add_axes([middle, top+height+landscape_height, 0.2, 0.05])
draw_label(ax, "Slow adjustment ($\it{Type}$ $\it{II}$)")

ax = fig.add_axes([0, bot+height+landscape_height, 0.05, 0.05])
draw_label(ax, "A")
ax = fig.add_axes([middle, bot+height+landscape_height, 0.2, 0.05])
draw_label(ax, "Rapid turnover ($\it{Type}$ $\it{I}$)")


# draw arrows and labels
ax = fig.add_axes([middle, bot+height, 0.2, landscape_height])
ax.arrow(0.1, 0.5, 0.8, 0, head_width=0.1, head_length=0.05, fc='k', ec='k')
ax.text(0.5, 0.6, r"$TTI$=" + str(rearrangement_ri)[:4] + r"$\pm$" + str(rearrangement_ri_std)[:4], horizontalalignment='center')
plt.axis("off")

ax = fig.add_axes([middle, top+height, 0.2, landscape_height])
ax.arrow(0.1, 0.5, 0.8, 0, head_width=0.1, head_length=0.05, fc='k', ec='k')
ax.text(0.5, 0.6, r"$TTI$=" + str(adjustment_ri)[:4] + r"$\pm$" + str(adjustment_ri_std)[:4], horizontalalignment='center')
plt.axis("off")


print "\nFINISHED"
#plt.savefig("figure_3e.png", dpi=300)
plt.savefig("random_"+str(s)+".png", dpi=300)
plt.savefig("figure_5.png", dpi=300)
#plt.savefig("figure_5.eps", dpi=600)
#plt.show()






