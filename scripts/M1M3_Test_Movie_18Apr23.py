#!/usr/bin/env python
# coding: utf-8
"""
M1M3 cell learning
Craig Lage - 18-Apr-23 \
The 17 tons of mirror are supported by 156 pneumatic actuators where 44 are single-axis and provide support only on the axial direction, 100 are dual-axis providing support in the axial and lateral direction, and 12 are dual-axis providing support in the axial and cross lateral directions. \
Positioning is provided by 6 hard points in a hexapod configuration which moves the mirror to a fixed operational position that shall be maintained during telescope operations. The remaining optical elements will be moved relative to this position in order to align the telescope optics. Support and optical figure correction is provided by 112 dual axis and 44 single axis pneumatic actuators. 
"""


import sys, time, os, asyncio, glob
from datetime import datetime
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.colors import LightSource as LS
import pickle as pkl
from astropy.time import Time, TimeDelta
import lsst.ts.cRIOpy.M1M3FATable as M1M3FATable
from lsst_efd_client import EfdClient

########################## SUBROUTINES #######################################

def actuatorLayout(ax, FATABLE):
    ax.set_xlabel("X position (m)")
    ax.set_ylabel("Y position (m)")
    ax.set_title("M1M3 Actuator positions and type\nHardpoints are approximate", fontsize=18)
    types = [['SAA','NA', 'o', 'Z', 'b'], ['DAA','+Y', '^', '+Y','g'], ['DAA','-Y', 'v', '-Y', 'cyan'], \
             ['DAA','+X', '>', '+X', 'r'], ['DAA','-X', '<', '-X', 'r']]
    for [type, orient, marker, label, color] in types: 
        xs = []
        ys = []
        for i in range(len(FATABLE)):
            x = FATABLE[i][M1M3FATable.FATABLE_XPOSITION]
            y = FATABLE[i][M1M3FATable.FATABLE_YPOSITION]
            if FATABLE[i][M1M3FATable.FATABLE_TYPE] == type and FATABLE[i][M1M3FATable.FATABLE_ORIENTATION] == orient:
                xs.append(x)
                ys.append(y)
            else:
                continue
        ax.scatter(xs, ys, marker=marker, color=color, s=200, label=label)        

    # Now plot approximate hardpoint location
    Rhp = 3.1 # Radius in meters
    for i in range(6):
        theta = 2.0 * np.pi / 6.0 * float(i)
        if i == 0:
            ax.scatter(Rhp * np.cos(theta), Rhp * np.sin(theta), marker='o', color='magenta', s=200, label='HP')
        else:
            ax.scatter(Rhp * np.cos(theta), Rhp * np.sin(theta), marker='o', color='magenta', s=200, label='_nolegend_')
    ax.legend(loc='lower left', fontsize=9)


def barChartZ(df, ax, FATABLE, index):
    ax.set_xlabel("X position (m)")
    ax.set_ylabel("Y position (m)")
    ax.set_zlabel("Force (nt)")
    ax.set_title("M1M3 Actuator Z forces", fontsize=18)

    lightsource = LS(azdeg=180, altdeg=78)
    greyColor = '0.9'
    colors = []
    xs = []
    ys = []
    for i in range(len(FATABLE)):
        x = FATABLE[i][M1M3FATable.FATABLE_XPOSITION]
        y = FATABLE[i][M1M3FATable.FATABLE_YPOSITION]
        xs.append(x)
        ys.append(y)
        if FATABLE[i][M1M3FATable.FATABLE_TYPE] == 'SAA':
            colors.append('blue'); colors.append('blue')
            colors.append(greyColor); colors.append(greyColor)
            colors.append(greyColor); colors.append(greyColor)
        else:
            if FATABLE[i][M1M3FATable.FATABLE_ORIENTATION] in ['+Y', '-Y']:
                colors.append('green'); colors.append('green')
                colors.append(greyColor); colors.append(greyColor)
                colors.append(greyColor); colors.append(greyColor)
            if FATABLE[i][M1M3FATable.FATABLE_ORIENTATION] in ['+X', '-X']:
                colors.append('red'); colors.append('red')
                colors.append(greyColor); colors.append(greyColor)
                colors.append(greyColor); colors.append(greyColor)

    zs = np.zeros([len(FATABLE)])
    for i in range(len(FATABLE)):
        name=f"zForces{i}"
        zs[i] = df.iloc[index][name]

    dxs = 0.2 * np.ones([len(FATABLE)])
    dys = 0.2 * np.ones([len(FATABLE)])
    bottom = np.zeros([len(FATABLE)])
    ax.bar3d(xs, ys, bottom, dxs, dys, zs, shade=True, alpha=0.5, lightsource=lightsource, color=colors)

    ax.set_zlim(0, 1500)
    ax.view_init(elev=30., azim=225)


def heatMapZ(df, ax, FATABLE, index):
    ax.set_xlabel("X position (m)")
    ax.set_ylabel("Y position (m)")
    ax.set_title("M1M3 Actuator Z forces (nt)", fontsize=18)

    types = [['SAA','NA', 'o', 'Z'], ['DAA','+Y', '^', '+Y'], ['DAA','-Y', 'v', '-Y'], ['DAA','+X', '>', '+X'], ['DAA','-X', '<', '-X']]

    for [type, orient, marker, label] in types: 
        xs = []
        ys = []
        zs = []
        for i in range(len(FATABLE)):
            x = FATABLE[i][M1M3FATable.FATABLE_XPOSITION]
            y = FATABLE[i][M1M3FATable.FATABLE_YPOSITION]
            if FATABLE[i][M1M3FATable.FATABLE_TYPE] == type and FATABLE[i][M1M3FATable.FATABLE_ORIENTATION] == orient:
                xs.append(x)
                ys.append(y)
                name=f"zForces{i}"
                zs.append(df.iloc[index][name])
            else:
                continue
        im = ax.scatter(xs, ys, marker=marker, c=zs, cmap='RdBu_r', vmin=800.0, vmax=1500, s=200, label=label) 
    plt.colorbar(im, ax=ax,fraction=0.055, pad=0.02, cmap='RdBu_r')


def lateralForces(df, ax, FATABLE, index, forceMax=1500):
    ax.set_xlabel("X position (m)")
    ax.set_ylabel("Y position (m)")
    ax.set_title("M1M3 lateral forces (nt)", fontsize=18)
    ax.set_xlim(-4.5,4.5)
    ax.set_ylim(-4.5,4.5)
    types = [['DAA','+Y', '^', '+Y','g'], ['DAA','-Y', 'v', '-Y', 'cyan'], \
             ['DAA','+X', '>', '+X', 'r'], ['DAA','-X', '<', '-X', 'r']]
    for [type, orient, marker, label, color] in types: 
        xs = []
        ys = []
        arrowXs = []
        arrowYs = []
        for i in range(len(FATABLE)):
            x = FATABLE[i][M1M3FATable.FATABLE_XPOSITION]
            y = FATABLE[i][M1M3FATable.FATABLE_YPOSITION]
            if FATABLE[i][M1M3FATable.FATABLE_TYPE] == type and FATABLE[i][M1M3FATable.FATABLE_ORIENTATION] == orient:
                xs.append(x)
                ys.append(y)
                if orient == '+X':
                    name = f"xForces{FATABLE[i][M1M3FATable.FATABLE_XINDEX]}"
                    arrowXs.append(df.iloc[index][name] / forceMax)
                    arrowYs.append(0.0)
                if orient == '-X':
                    name = f"xForces{FATABLE[i][M1M3FATable.FATABLE_XINDEX]}"
                    arrowXs.append(-df.iloc[index][name] / forceMax)
                    arrowYs.append(0.0)
                if orient == '+Y':
                    name = f"yForces{FATABLE[i][M1M3FATable.FATABLE_YINDEX]}"
                    arrowXs.append(0.0)
                    arrowYs.append(df.iloc[index][name] / forceMax)
                if orient == '-Y':
                    name = f"yForces{FATABLE[i][M1M3FATable.FATABLE_YINDEX]}"
                    arrowXs.append(0.0)
                    arrowYs.append(-df.iloc[index][name] / forceMax)
            else:
                continue
        ax.scatter(xs, ys, marker=marker, color=color, s=50, label=label) 
        for ii in range(len(xs)):
            ax.arrow(xs[ii], ys[ii], arrowXs[ii], arrowYs[ii], color=color)

    ax.plot([-4.0,-3.0], [-4.0,-4.0], color='g')
    ax.text(-4.0, -4.3, f"{forceMax} nt")

def buildMovie(df, FATABLE):
    for n in range(0, len(df), 5):
        fig = plt.figure(figsize=(16,16))
        ax1 = fig.add_subplot(2,2,1)
        actuatorLayout(ax1, FATABLE)
        ax2 = fig.add_subplot(2,2,2, projection='3d')
        barChartZ(df, ax2, FATABLE, n)
        ax3 = fig.add_subplot(2,2,3)
        lateralForces(df, ax3, FATABLE, n, forceMax=25)
        ax4 = fig.add_subplot(2,2,4)
        heatMapZ(df, ax4, FATABLE, n)
        plt.savefig(f"/home/craiglagegit/movies/m1m3_test_18apr23/Frame_{n:05d}.png")
        plt.close(fig)


####################### MAIN PROGRAM ###################################



async def main():
    client = EfdClient('summit_efd')
    FATABLE = M1M3FATable.FATABLE
    # Times to start looking at encoder values
    start = Time("2023-04-18 16:10:00Z", scale='utc')
    end = Time("2023-04-18 16:15:00Z", scale='utc')
    forces = await client.select_time_series("lsst.sal.MTM1M3.appliedForces", "*", start, end)
    buildMovie(forces, FATABLE)

asyncio.run(main())



