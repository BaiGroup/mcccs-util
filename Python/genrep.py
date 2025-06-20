#!/usr/bin/env python

import argparse
import os
import math
import simutil as su
import statistics as stats


def parse_args():
    parser = argparse.ArgumentParser(description='MCCCS data extractor.')
    parser.add_argument('-b', '--nbox', action='store', default=2, type=int)
    parser.add_argument('-e', '--equilibration',
                        action='store', default=2, type=int)
    parser.add_argument('-f', '--folder-pattern',
                        action='store', required=True)
    parser.add_argument('-i', '--ignore-case', action='store_true',
                        default=True)
    parser.add_argument('--ignore-error', action='store_true')
    parser.add_argument('-m', '--nmolty', action='store', default=1, type=int)
    parser.add_argument('-o', '--output', action='store', default='rep.txt')
    parser.add_argument('-p', '--plot-trajectory', action='store_true')
    parser.add_argument('-r', '--run-pattern', action='store',
                        default=r'RUN(\d+)$')
    args, unknown = parser.parse_known_args()
    return args


def traj(ax, result):
    if not hasattr(traj, "firstBlock"):
        traj.firstBlock = 0
    block = np.arange(traj.firstBlock, traj.firstBlock +
                      len(result.energyTrajectory[1]))
    traj.firstBlock += len(result.energyTrajectory[1])
    for ibox in range(result.nbox):
        ax[ibox].plot(block, result.energyTrajectory[ibox+1],
                      'o:', label='energyTrajectory')
        ax[2+ibox].plot(block, result.densityTrajectory[ibox+1],
                        'o:', label='densityTrajectory')
        ax[4+ibox].plot(block, result.pressureTrajectory[ibox+1],
                        'o:', label='pressureTrajectory')
        return


if __name__ == '__main__':
    args = parse_args()
    print(args)

    plot_trajectory = args.plot_trajectory
    nbox = args.nbox
    nmolty = args.nmolty

    # X1: mole fraction for box 1
    # X2: mole fraction for box 2
    # X3: mole fraction for box 3
    fou = open(args.output, 'w')
    fou.write('#Name NCycles T[K] P_vap[kPa] E')
    for imolty in range(nmolty):
        fou.write(' X1_{} E X2_{} E X3_{} E'.format(imolty+1, imolty+1,
                                                    imolty+1))
    for ibox in range(nbox):
        fou.write(' rho_{}[g/mL] E Z_{} E'.format(ibox+1, ibox+1))
    fou.write('\n')

    avg = stats.Stat('Press', X1=nmolty, X2=nmolty, X3=nmolty, Dens=nbox,
                     Z=nbox)
    avg.register('Avg', stats.mean)
    avg.register('SE', stats.stderr)

    if plot_trajectory:
        import numpy as np
        import matplotlib as mpl
        import matplotlib.pyplot as plt
        mpl.rcParams['figure.figsize'] = (3.25, 13)

    sim = su.ProcDir(args.folder_pattern, args.run_pattern, args.ignore_case)

    for curSim in sim:
        print(curSim)
        # averages for each run
        RunCycles = []  # total number of cycles performed
        avg.clear()

        if plot_trajectory:
            fig, ax = mpl.pyplot.subplots(nrows=6, ncols=1, sharex=True)
            fig.subplots_adjust(left=0.15, right=0.95,
                                bottom=0.05, top=0.95, hspace=0)
            for oneax in ax:
                oneax.ticklabel_format(
                    style='scientific', scilimits=(0, 0), axis='y')

        run = sim.getRuns()
        for pathname in run:
            os.chdir(pathname)

            (suffix, last, cur) = su.getSuffix('log')
            print("%s -> %s%d" % (pathname, suffix, last))

            RunCycles.append(0)
            avg.Press.append(0.0)
            for ibox in range(nbox):
                avg.Dens[ibox].append(0.0)
                avg.Z[ibox].append(0.0)

            for imol in range(nmolty):
                avg.X1[imol].append(0.0)
                avg.X2[imol].append(0.0)
                avg.X3[imol].append(0.0)

            if suffix == 'equil':
                st = last-args.equilibration+1
                last = args.equilibration
            else:
                st = 1

            for ifile in range(st, last+st):
                try:
                    result = su.procRun("run."+suffix+str(ifile))
                    temp = result.temp
                    if plot_trajectory:
                        traj(ax, result)

                    RunCycles[run.count-1] += result.ncycles
                    avg.Press[run.count-1] += result.press[result.nbox] * \
                        result.ncycles

                    for ibox in range(result.nbox):
                        avg.Dens[ibox][run.count -
                                       1] += result.dens[ibox+1]*result.ncycles
                        avg.Z[ibox][run.count -
                                    1] += result.Z[ibox+1]*result.ncycles

                    for imol in range(nmolty):
                        X1 = result.molfrac[1][imol+1]
                        if nbox > 1:
                            X2 = result.molfrac[2][imol+1]
                        if nbox > 2:
                            X3 = result.molfrac[3][imol+1]
                        avg.X1[imol][run.count-1] += X1*result.ncycles
                        if nbox > 1:
                            avg.X2[imol][run.count-1] += X2*result.ncycles
                        if nbox > 2:
                            avg.X3[imol][run.count-1] += X3*result.ncycles
                except IOError as e:
                    print(e)
                    if args.ignore_error:
                        continue

            avg.Press[run.count-1] /= RunCycles[run.count-1]

            for ibox in range(nbox):
                avg.Dens[ibox][run.count-1] /= RunCycles[run.count-1]
                avg.Z[ibox][run.count-1] /= RunCycles[run.count-1]

            for imol in range(nmolty):
                avg.X1[imol][run.count-1] /= RunCycles[run.count-1]
                avg.X2[imol][run.count-1] /= RunCycles[run.count-1]
                avg.X3[imol][run.count-1] /= RunCycles[run.count-1]

        os.chdir(curSim)
        if plot_trajectory:
            fig.savefig('trajectory.pdf')

        # averages
        avg.update(RunCycles)

        fou.write('%s %.2f %3.2f %4.3f %1.3f' % (sim.simMatch.group(1), float(
            sum(RunCycles))/run.count, temp, avg.AvgPress, avg.SEPress))

        for imol in range(nmolty):
            fou.write(' %.5G %.9G %.5G %.9G %.5G %.9G' % (
                avg.AvgX1[imol], avg.SEX1[imol], avg.AvgX2[imol], avg.SEX2[imol], avg.AvgX3[imol], avg.SEX3[imol]))

        for ibox in range(nbox):
            fou.write(' %1.5G %1.5G %1.5G %1.5G' % (avg.AvgDens[ibox], avg.SEDens[ibox], avg.AvgZ[ibox], avg.SEZ[ibox]))
        fou.write('\n')

    fou.close()
