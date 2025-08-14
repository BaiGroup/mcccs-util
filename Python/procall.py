#!/usr/bin/env python

''' A program to set up input files and optimize simulation parameters

You need to modify the code according to the particular circumstances of your simulations.
Be ultra cautious when you use this kind of script: go through the file and make sure you understand each and every line.

Examples: For a directory hierarchy named like `~/projects/adsorption/Polyols/MFI/13PrDO/333K/512_1/RUN1`
1) During equilibration:
    /path/to/procall.py -b 3 -m 4 --swatch-pair 2 --bias-equal '2 3' --bias-nmol 1 -t 24 -f '/Polyols/MFI/.*/([\.\d]+)K/(\d+)_(\d+)$'
2) Transition from equilibration to production:
    /path/to/procall.py -b 3 -m 4 --swatch-pair 2 --bias-equal '2 3' --bias-nmol 1 -t 24 -f '/Polyols/MFI/.*/([\.\d]+)K/(\d+)_(\d+)$' --average
'''

import argparse
import math
import os
import random
import shutil
import simutil as su


def parse_args():
    parser = argparse.ArgumentParser(description='MCCCS input adjustor.')
    parser.add_argument('--average', action='store_true')
    parser.add_argument('-b', '--nbox', action='store', default=2, type=int)
    parser.add_argument('--bias-equal', action='store', default='')
    parser.add_argument('--bias-nmol', action='store', default=1, type=float)
    parser.add_argument('-f', '--folder-pattern',
                        action='store', required=True)
    parser.add_argument('-i', '--ignore-case', action='store_true',
                        default=True)
    parser.add_argument('-m', '--nmolty', action='store', default=1, type=int)
    parser.add_argument('--modify', action='store', default='equil1 prod1')
    parser.add_argument('-n', '--ncycles', action='store', type=int)
    parser.add_argument('--reset-swatch-pair', action='store_true')
    parser.add_argument('-r', '--run-pattern', action='store',
                        default=r'RUN(\d+)$')
    parser.add_argument('--swatch-pair', action='store', default=1, type=int)
    parser.add_argument('-t', '--time', action='store', type=float)
    parser.add_argument('--time-buffer-ratio', action='store',
                        default=0.92, type=float)
    args, unknown = parser.parse_known_args()
    return args


class TopmonInput():
    def __init__(self, args):
        self.args = args
        self.nmolty = args.nmolty
        self.nbox = args.nbox
        self.nswaty = args.swatch_pair
        self.bias_equal = [int(x) for x in args.bias_equal.split()]
        self.bias_nmol = args.bias_nmol
        self.nboxpair = int(self.nbox*(self.nbox-1)/2)
        if self.nboxpair < 1:
            self.nboxpair = 1
        self.avgParam = su.avgResult(nbox=self.nbox, nmolty=self.nmolty, nmoltyeta=self.nmolty,
                                     nswapbp=self.nmolty, nswatbp=args.swatch_pair, npairbox=self.nboxpair)

    def average(self):
        avgParam = self.avgParam
        avgParam.average()

        self.ncycles = avgParam.ncycles
        self.pmswat = avgParam.pmswat
        self.pmswap = avgParam.pmswap
        self.pmcb = avgParam.pmcb
        self.pmtra = avgParam.pmtra
        self.swatch_pair_prob = avgParam.pmsatc
        self.swatch_box_prob = avgParam.pmswtcb
        self.swap_box_prob = avgParam.pmswapb
        self.swap_mol_prob = avgParam.pmswmt
        self.eta = avgParam.eta

    def modify(self, stage=None, last=None, unfinished=False):
        if self.args.time is None:
            if unfinished:
                print(
                    'Provide the walltime used for the unfinished run via the "-t" option')
                raise RuntimeError
            else:
                self.args.time = 24.0

        time_target = self.args.time*self.args.time_buffer_ratio
        if unfinished:
            ncyclesn = su.procCFG("save-config.1", "fort.77")
            ncyclesn = su.estimateRuntime(
                timeTarget=time_target,
                ncycles=ncyclesn,
                timeElapsed=self.args.time)
        else:
            result = su.procRun("run."+stage+str(last))
            assert result.nbox == self.nbox
            assert result.nmolty == self.nmolty
            assert len(result.pmswtcb) == self.nswaty + 1
            timeElapsed = su.getTime("log")
            ncyclesn = su.estimateRuntime(timeTarget=time_target,
                                          ncycles=result.ncycles,
                                          timeElapsed=timeElapsed/3600)
        if self.args.ncycles is not None:
            if ncyclesn < self.args.ncycles:
                print('Requested ncycles = {} will likely run over walltime. ncycles = {} correspond to {} hours'.format(
                    self.args.ncycles, ncyclesn, time_target))
            ncyclesn = self.args.ncycles

        self.ncycles = ncyclesn
        self.avgParam.ncycles.append(ncyclesn)

        if not unfinished and stage == 'equil':
            self.pmvol = result.pmvol

            nboxpair = self.nboxpair

            if result.pmswat != 0:
                self.pmswat = (result.pmswat*result.ncycles /
                               result.ntswat + self.pmvol)

                fcfg_in = open("fort.4."+self.args.modify.split()[0])
                fdat = fcfg_in.readlines()
                fcfg_in.close()

                (lineNo, colNo) = su.get_namelist_keyword(fdat, "mc_swatch",
                                                          "nswaty")
                nswaty = int(fdat[lineNo][colNo+1:].strip())
                assert nswaty == self.nswaty

                if self.args.reset_swatch_pair:
                    swatch_pair_prob = [1.0/self.nswaty] * self.nswaty
                    swatch_pair_prob = [sum(swatch_pair_prob[:i+1])
                                        for i in range(len(swatch_pair_prob))]
                else:
                    swatch_pair_prob = su.estimateSwatchPairProb(range(self.nswaty), result.nswat,
                                                                 result.pmsatc, nboxpair=nboxpair)
                self.swatch_pair_prob = swatch_pair_prob
                for i, p in enumerate(swatch_pair_prob):
                    self.avgParam.pmsatc[i].append(p)

                self.swatch_box_prob = []
                for ipair in range(self.nswaty):
                    self.swatch_box_prob.append(su.estimateSwatchBoxProb(
                        result.nswat[ipair*nboxpair*2:(ipair+1)*nboxpair*2+1],
                        result.pmswtcb[ipair+1], boxpair=range(nboxpair)
                    ))
                    for i, p in enumerate(self.swatch_box_prob[ipair]):
                        self.avgParam.pmswtcb[ipair][i].append(p)
            else:
                self.pmswat = result.pmvol

            if result.pmswap != 0:
                try:
                    self.pmswap = (result.pmswap*result.ncycles /
                                   result.ntswap + self.pmswat)
                except ZeroDivisionError:
                    self.pmswap = 0.6

                if self.pmswap > 0.6:
                    print('Swap probabilities capped: pmswat = {:.5G}, pmswap = {:.5G}'.format(self.pmswat, self.pmswap))
                    self.pmswat = (self.pmvol +
                                   (self.pmswat-self.pmvol) *
                                   (0.6-self.pmvol) /
                                   (self.pmswap-self.pmvol))
                    self.pmswap = 0.6

                self.swap_box_prob = []
                for imol in range(result.nmolty):
                    nmol = sum([x[imol+1] for x in result.nmol[1:]])
                    self.swap_box_prob.append(
                        su.estimateSwapBoxProb(result.nswap[imol*nboxpair*2:(imol+1)*nboxpair*2+1],
                                               result.pmswapb[imol+1],
                                               boxpair=range(nboxpair)))
                    for i, p in enumerate(self.swap_box_prob[imol]):
                        self.avgParam.pmswapb[imol][i].append(p)
                    result.nswap[imol*nboxpair*2+1:(imol+1)*nboxpair*2+1] = [
                        float(x) / nmol
                        for x in result.nswap[imol*nboxpair*2+1:(imol+1)*nboxpair*2+1]
                    ]

                self.swap_mol_prob = su.estimateMolSwapProb(range(result.nmolty), result.nswap,
                                                            result.pmswmt, nboxpair=nboxpair)
                for i, p in enumerate(self.swap_mol_prob):
                    self.avgParam.pmswmt[i].append(p)
            else:
                self.pmswap = self.pmswat

            if result.pmswat != 0 or result.pmswap != 0:
                # adjust biasing potentials
                self.eta = [[] for ibox in range(self.nbox)]
                for imol in range(1, result.nmolty+1):
                    ref_box = 0
                    if imol in self.bias_equal:
                        # ref_box is the simulation box containing imol but no biasing potentials
                        # ref_box will be zero if all imol molecules are in another box, e.g., from too negative biasing potentials
                        for ibox in range(1, result.nbox+1):
                            if result.nmol[ibox][imol] != 0 and result.eta[ibox][imol] == 0:
                                ref_box = ibox
                                break

                    for ibox in range(1, result.nbox+1):
                        if result.eta[ibox][imol] == 0:
                            # adjust only if the simulations ran with a biasing potential for imol in ibox
                            # set eta to some small value (e.g., 0.1) to signal the desire for adjustment
                            self.eta[ibox-1].append(0.0)
                        elif result.nmol[ibox][imol] == 0:
                            # cannot use the formula below
                            delta = abs(result.eta[ibox][imol]) / 2
                            if delta < result.temp:
                                delta = result.temp
                            self.eta[ibox-1].append(result.eta[ibox][imol] - delta)
                        else:
                            if imol in self.bias_equal:
                                if ref_box == 0:
                                    # cannot use the formula below
                                    delta = abs(result.eta[ibox][imol]) / 2
                                    if delta < result.temp:
                                        delta = result.temp
                                    self.eta[ibox-1].append(result.eta[ibox][imol] + delta)
                                else:
                                    # \Delta \mu_{ri}^{new} = \mu_{i} - \mu_{r} = \Delta \mu_{ri}^{old} + T \ln (N_i / N_r)
                                    # \mu_{r} will be set to zero from the logic earlier
                                    self.eta[ibox-1].append(
                                        result.eta[ibox][imol] - result.eta[ref_box][imol] +
                                        result.temp * math.log(result.nmol[ibox][imol] / result.nmol[ref_box][imol]))
                            else:
                                # \Delta \mu_{i}^{new} = \Delta \mu_{i}^{old} + T \ln (N_i / N_target)
                                self.eta[ibox-1].append(result.eta[ibox][imol] + result.temp * math.log(result.nmol[ibox][imol] / self.bias_nmol))

                        self.avgParam.eta[ibox-1][imol-1].append(self.eta[ibox-1][imol-1])

            if result.pmcb != 0:
                self.pmcb = (1+2*self.pmswap) / 3
            else:
                self.pmcb = self.pmswap

            self.pmtra = (1+self.pmcb) / 2
            self.avgParam.pmswat.append(self.pmswat)
            self.avgParam.pmswap.append(self.pmswap)
            self.avgParam.pmcb.append(self.pmcb)
            self.avgParam.pmtra.append(self.pmtra)

        self.save(write_mc_param=(not unfinished and stage == 'equil'))

    def save(self, write_mc_param=True):
        suffices = []
        for s in self.args.modify.split():
            suffices.append(s)

        for suffix in suffices:
            fcfg_in = open("fort.4."+suffix)
            fdat = fcfg_in.readlines()
            fcfg_in.close()

            su.set_namelist_keyword(fdat, "mc_shared", "seed",
                                    int(random.random()*1e6))
            su.set_namelist_keyword(fdat, "mc_shared", "nstep", self.ncycles)
            ncycles_per_block = int(self.ncycles/10)
            su.set_namelist_keyword(
                fdat, "mc_shared", "iprint", ncycles_per_block)
            su.set_namelist_keyword(
                fdat, "mc_shared", "iblock", ncycles_per_block)
            su.set_namelist_keyword(
                fdat, "mc_shared", "imv", ncycles_per_block)

            if write_mc_param:
                if self.pmswat > self.pmvol:
                    su.set_namelist_keyword(fdat, "mc_swatch", "pmswat", "{:.5e}".
                                            format(self.pmswat))  # swatch prob
                    su.set_namelist_keyword(fdat, "mc_swatch", "pmsatc",
                                            (" {:.5e}"*self.nswaty).format(*self.swatch_pair_prob))
                    for ipair in range(self.nswaty):
                        # swatch box pair prob
                        su.set_mcswatch(
                            fdat, ipair+1, self.swatch_box_prob[ipair])

                if self.pmswap > self.pmswat:
                    su.set_namelist_keyword(fdat, "mc_swap", "pmswap", "{:.5e}".
                                            format(self.pmswap))  # swap prob
                    for imol in range(self.nmolty):
                        # swap box pair prob
                        su.set_mcswap(fdat, imol+1, self.swap_box_prob[imol])

                    su.set_namelist_keyword(fdat, "mc_swap", "pmswmt",
                                            (" {:.5e}"*self.nmolty).format(*self.swap_mol_prob))

                if self.pmcb > self.pmswap:
                    su.set_namelist_keyword(fdat, "mc_cbmc", "pmcb", "{:.5e}".
                                            format(self.pmcb))  # cbmc prob

                su.set_namelist_keyword(fdat, "mc_simple", "pmtra", "{:.5e}".
                                        format(self.pmtra))  # translation prob

                if self.pmswat > self.pmvol or self.pmswap > self.pmswat:
                    # biasing potentials
                    line_start = su.get_section(fdat, 'MOLECULE_TYPE')
                    line_end = su.get_section(fdat, 'END MOLECULE_TYPE')
                    for imol in range(self.nmolty):
                        line_start = su.set_next_eta(fdat, [str(self.eta[ibox][imol]) for ibox in range(self.nbox)], line_start+1, line_end)

            fcfg_out = open("fort.4."+suffix, 'w')
            fcfg_out.writelines(fdat)
            fcfg_out.close()


if __name__ == '__main__':
    args = parse_args()
    print(args)

    sim = su.ProcDir(args.folder_pattern, args.run_pattern, args.ignore_case)
    for curSim in sim:
        topmon_input = TopmonInput(args)

        run = sim.getRuns()
        for pathname in run:
            os.chdir(pathname)

            (suffix, last, cur) = su.getSuffix('log')
            print("%s -> %s%d -> %d" % (pathname, suffix, last, cur))

            if os.path.exists('final-config'):
                print("volume too small!")
                last -= 1
                shutil.move('final-config', 'fort.77')
                su.batchMove(suffix, last, nbox=args.nbox)
            elif os.path.exists('run1a.dat'):
                print("not finished!")
                topmon_input.modify(suffix, last, unfinished=True)
                if suffix == 'equil':
                    last = cur
                    shutil.move('save-config.1', "fort.77")
                    su.batchMove(suffix, last, nbox=args.nbox)
            else:
                topmon_input.modify(suffix, last)

            if args.average:
                topmon_input.avgParam.readMaximumDisplacement('fort.77')

        if args.average:
            print('#####Averaging#####')
            print(topmon_input.avgParam)
            try:
                topmon_input.average()

                run = sim.getRuns()
                for pathname in run:
                    os.chdir(pathname)

                    (suffix, last, cur) = su.getSuffix('log')
                    print("%s -> %s%d -> %d" % (pathname, suffix, last, cur))

                    topmon_input.save(write_mc_param=True)
                    topmon_input.avgParam.writeMaximumDisplacement('fort.77')
            except (TypeError, ValueError) as e:
                print(e)
                print('skip averaging..')
                continue
            finally:
                print("###################\n")
