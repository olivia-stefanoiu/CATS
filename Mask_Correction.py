#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import ROOT
import math
from sklearn.linear_model import LinearRegression
from collections import Counter
from numpy import arctan, degrees, radians, cos, sin
from sklearn.linear_model import RANSACRegressor
from array import array
from numpy import arctan, degrees
from ROOT import TGraph, TLine, kRed, kAzure


# In[ ]:


def RANSAC_fit_line(x, y, threshold, show):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)

    x_2d = x.reshape(-1, 1)

    model = RANSACRegressor(residual_threshold=threshold)
    model.fit(x_2d, y)

    inliers = model.inlier_mask_
    slope = float(model.estimator_.coef_[0])
    intercept = float(model.estimator_.intercept_)

    if show:
        x_min, x_max = float(x.min()), float(x.max())
        y_min, y_max = float(y.min()), float(y.max())

        hist = ROOT.TH2F("hist_ransac", "", 200, x_min, x_max, 200, y_min, y_max)
        for index in range(len(x)):
            hist.Fill(float(x[index]), float(y[index]))

        canvas = ROOT.TCanvas("c_ransac", "RANSAC Inliers", 800, 600)
        hist.Draw("COLZ")

        x_inliers = x[inliers]
        y_inliers = y[inliers]

        graph = ROOT.TGraph(len(x_inliers), x_inliers, y_inliers)
        graph.SetMarkerColor(ROOT.kRed)
        graph.SetMarkerStyle(20)
        graph.SetMarkerSize(0.3)
        graph.Draw("P SAME")

        y1 = slope * x_min + intercept
        y2 = slope * x_max + intercept
        line = ROOT.TLine(x_min, y1, x_max, y2)
        line.SetLineColor(ROOT.kAzure + 7)
        line.SetLineWidth(2)
        line.Draw("SAME")

        canvas.Update()
        input("Press Enter to continue...")

    return inliers, slope, intercept


# In[ ]:


x_vals = []
y_vals = []

text_file_handle = open("/home/olivia/Desktop/scripts/CATS/CATS2_centroid.txt", "r")
for line in text_file_handle:
    x_value, y_value = line.split()
    x_vals.append(float(x_value))
    y_vals.append(float(y_value))
text_file_handle.close()

#pentru CATS1 
#selectam punctele pentru corectia verticala
selected_pairs = [(x_value, y_value) for x_value, y_value in zip(x_vals, y_vals) if (-20 <= x_value <= 7) and (0.6 <= y_value <= 5)]
x_selected = [pair[0] for pair in selected_pairs]
y_selected = [pair[1] for pair in selected_pairs]
#selectam punctele pentru corectia orizontala 
x_selected2 = [x for x, y in zip(x_vals, y_vals) if -8.5 <= x <= -5.7]
x_mean = sum(x_selected2) / len(x_selected2)
inliers,slope,intercept =RANSAC_fit_line(x_selected,y_selected,0.0001,True)
# slope =0.02250442496464144 
# intercept = 2.7132336197449174
# x_mean =-7.046736916766848

print(slope,intercept)
print(x_mean)

#CATS2
# selected_pairs = [(x_value, y_value) for x_value, y_value in zip(x_vals, y_vals) if (-15 <= x_value <= 14) and (-2.5 <= y_value <= 1)]
# x_selected = [pair[0] for pair in selected_pairs]
# y_selected = [pair[1] for pair in selected_pairs]
# #selectam punctele pentru corectia orizontala 
# x_selected2 = [x for x, y in zip(x_vals, y_vals) if -1.4 <= x <= -0.2]
# x_mean = sum(x_selected2) / len(x_selected2)
# inliers,slope,intercept =RANSAC_fit_line(x_selected,y_selected,0.0001,True)
# slope =0.005910027638249295
# intercept = -0.7013276374635101
# x_mean = -0.8083158756358418
# print(slope,intercept)
# print(x_mean)

x_center = x_mean

x_shifted = []
y_shifted = []

rotation_angle = math.atan(slope)
cos_angle = math.cos(-rotation_angle)
sin_angle = math.sin(-rotation_angle)

x_shifted = []
y_shifted = []

for x_value, y_value in zip(x_vals, y_vals):
    x_new = x_value + abs(x_center)
    #-9.5 cats2 -6.5 cats1
    y_new = y_value - (slope * (-6.5) + intercept)

    x_rotated = x_new * cos_angle - y_new * sin_angle
    y_rotated = x_new * sin_angle + y_new * cos_angle

    x_shifted.append(x_rotated)
    y_shifted.append(y_rotated)

h_shifted = ROOT.TH2F("h", "h;X shifted;Y shifted", 2000, min(x_shifted), max(x_shifted), 2000, min(y_shifted), max(y_shifted))
for point_index in range(len(x_shifted)):
    h_shifted.Fill(x_shifted[point_index], y_shifted[point_index])
c = ROOT.TCanvas("c", "c", 900, 700)
h_shifted.Draw("COLZ")
c.Update()
input()




