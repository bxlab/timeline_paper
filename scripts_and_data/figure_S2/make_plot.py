#!/usr/bin/env python
print "loading python packages..."
import sys, getopt, os
import matplotlib.pyplot as plt
import numpy as np
import math


def inches_to_cm(inch):
	return inch * 25.4

def F_to_C (temp):
	return (temp-32.0)*5.0/9.9


def load_data(filename):
	points=[]
	years=[]
	months=[]
	dates=[]
	max_temps=[]
	min_temps=[]
	max_humis=[]
	min_humis=[]
	precips=[]
	
	ct=0
	for line in open(filename):
		cut=line.strip().split("\t")
		if line.startswith("20") and len(cut[0])==4:
			year = int(cut[0])
			continue
		if cut[1]=="Max":
			month = cut[0]
			continue
		else:
			date = int(cut[0])
			ct+=1
	
			temp_max = F_to_C(float(cut[1]))
			temp_min = F_to_C(float(cut[3]))
		
			humi_max = float(cut[7])
			humi_min = float(cut[9])
		
			precip = inches_to_cm(float(cut[17]))
	
			points.append(ct)		
			years.append(year)
			months.append(month)
			dates.append(date)
			max_temps.append(temp_max)
			min_temps.append(temp_min)
			max_humis.append(humi_max)
			min_humis.append(humi_min)
			precips.append(precip)
			if precip>0:
				print line.strip()
	return points, years, months, dates, max_temps, min_temps, max_humis, min_humis, precips


def draw_time_axis(points, years, months, dates, ax):
	ax.get_yaxis().set_visible(False)
	ax.set_xticks(points)
	ax.tick_params(axis=u'both', which=u'both',length=0)
	ax.set_xticklabels([""]*len(points))
	
	colors=["gold", "cyan", "royalblue", "magenta"]
	for x,year in enumerate(years):
		c=colors[year-2014]
		ax.axvline(x=x, c=c)
	
	labels={"Jan":None, "Apr":None, "Jul":None, "Oct":None}
	for x,month in enumerate(months):
		if dates[x]==1:
			ax.axvline(x=x, ymin=0.9, ymax=1, c='k')
			#if month in labels:
			#	ax.text(x, 0.6, month, rotation=-45)
		if month=="May" and dates[x]==15:
			ax.text(x, 0.2, years[x], fontsize=18)
	ax.text(1100, 0.50, "2017", rotation=90, fontsize=18)


def draw_precipitation(points, months, precips, ax):
	ax.set_xlim(0,len(points))
	ax.bar(points, precips, width=10)
	for x,month in enumerate(months):
		if dates[x]==1:
			ax.axvline(x=x, c='k', linestyle='--', alpha=0.2, linewidth=0.5)
	ax.get_xaxis().set_visible(False)
	ax.set_ylabel("Precipitation (mm)", fontsize=12)
	ax.set_ylim(0,24)
	plt.grid(axis="y", linestyle='--', alpha=0.3)


def draw_humidity(points, months, max_humis, min_humis, ax):
	ax.set_xlim(0,len(points))
	ax.scatter(points, max_humis, s=10, alpha=0.1, c="r")
	ax.scatter(points, min_humis, s=10, alpha=0.1, c="b")
	for x,month in enumerate(months):
		if dates[x]==1:
			ax.axvline(x=x, c='k', linestyle='--', alpha=0.2, linewidth=0.5)
	ax.get_xaxis().set_visible(False)
	ax.set_ylabel("Relative humidity (%)", fontsize=12)
	ax.set_ylim(0,105)
	plt.grid(axis="y", linestyle='--', alpha=0.3)


def draw_temperature(points, months, max_temps, min_temps, ax):
	ax.set_xlim(0,len(points))
	ax.scatter(points, max_temps, s=10, alpha=0.1, c="r")
	ax.scatter(points, min_temps, s=10, alpha=0.1, c="b")
	for x,month in enumerate(months):
		if dates[x]==1:
			ax.axvline(x=x, c='k', linestyle='--', alpha=0.2, linewidth=0.5)
	ax.get_xaxis().set_visible(False)
	ax.set_ylabel("Temperature (Celcius)", fontsize=12)
	ax.set_ylim(0,34)
	plt.grid(axis="y", linestyle='--', alpha=0.3)

def add_expeditions(points, days, months, years, ax1, ax2):
	S1=[]
	S2=[]
	for i,x in enumerate(points):
		day = days[i]
		month = months[i]
		year = years[i]
		if month=="Sep" and year==2014 and day==1:
			S1.append(x)
		if month=="Jun" and year==2015 and day==1:
			S1.append(x)
		if month=="Feb" and year==2016 and day==8:
			S1.append(x)
			S2.append(x)
		if month=="Feb" and year==2017 and day==20:
			S1.append(x)
			S2.append(x)
		if month=="Jul" and year==2016 and day==11:
			S2.append(x)
		if month=="Oct" and year==2016 and day==20:
			S2.append(x)

	for x in S1:
		ax1.arrow( x, 0.5, 0.0, 0.2, fc="k", ec="k", width=5, head_width=15, head_length=0.1 , zorder=100)
	for x in S2:
		ax2.arrow( x, 5, 0.0, -3, fc="w", ec="k", width=5, head_width=15, head_length=1 , zorder=100)
			





# MAIN START
font = {'family': 'arial', 'weight': 'normal', 'size': 12}
plt.rc('font', **font)

print "loading weather data..."
points, years, months, dates, max_temps, min_temps, max_humis, min_humis, precips = load_data("2014-2017_weather.tab")


print "plotting..."
fig = plt.figure(figsize=(10, 10))

print "drawing axis..."
ax1 = fig.add_axes([0.1, 0.05, 0.85, 0.1])
draw_time_axis(points, years, months, dates, ax1)


print "drawing precipitation..."
ax2 = fig.add_axes([0.1, 0.15, 0.85, 0.2])
draw_precipitation(points, months, precips, ax2)
ax2.annotate("C", xy=(-0.08, 0.95), xycoords="axes fraction", fontsize=26)


print "drawing humidity..."
ax3 = fig.add_axes([0.1, 0.35, 0.85, 0.3])
draw_humidity(points, months, max_humis, min_humis, ax3)
ax3.annotate("B", xy=(-0.08, 0.95), xycoords="axes fraction", fontsize=26)


print "drawing temperature..."
ax4 = fig.add_axes([0.1, 0.65, 0.85, 0.3])
draw_temperature(points, months, max_temps, min_temps, ax4)
ax4.annotate("A", xy=(-0.08, 0.95), xycoords="axes fraction", fontsize=26)


print "adding sampling dates..."
add_expeditions(points, dates, months, years, ax1, ax2)


plt.subplots_adjust(left=0.1, right=0.95, top=0.9, bottom=0.1)
#plt.savefig("figure_2.eps", dpi=600)
plt.savefig("figure_S2.png", dpi=600)
#plt.show()


