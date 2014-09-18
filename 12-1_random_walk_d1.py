#! /usr/bin/env python
# -*- coding:utf-8 -*-
#
# written by Shotaro Fujimoto, May 2014.

import numpy as np
import matplotlib.pyplot as plt


class RandomWalk1():

    def __init__(self, prob=0.7, l=1, nwalkers=1000, x0=0):
        """ Initial function in RandomWalk1.

        prob     : probability that a particle moves right
        l        : step length
        nwalkers : number of trials
        x0       : initial position
        """
        self.prob = prob
        self.l = l
        self.nwalkers = nwalkers
        self.x0 = x0

    def random_walk_d1(self, N):
        """ Caluculate the displacements of each walkers.

        N : A list of walk steps
        """
        x = np.zeros([self.nwalkers, max(N)], 'i')
        p = np.random.random([self.nwalkers, max(N) - 1])
                             # generate random number in [0,1)
        prob = self.prob
        l = self.l
        x0 = self.x0

        for n in range(self.nwalkers):
            x[n][0] = x0
            for i in range(1, max(N)):
                d = +l if p[n][i - 1] < prob else -l
                x[n][i] = x[n][i - 1] + d
        self.x = x
        self.N = N

    def calc_ave(self):
        """ Caluculate the average of displacements after max(N) steps.

        You can call the results by "self.N", "self.x_ave", and "self.x_2_ave"
        """
        x = self.x
        x_ave = np.zeros(len(self.N))
        for i, nvalue in enumerate(self.N):
            x_ave[i] = sum(
                [x[n][nvalue - 1] * 1. for n in xrange(self.nwalkers)]) / self.nwalkers

        x_2_ave = np.zeros(len(self.N))
        for i, nvalue in enumerate(self.N):
            x_2_ave[i] = sum(
                [x[n][nvalue - 1] ** 2. for n in xrange(self.nwalkers)]) / self.nwalkers

        self.x_ave = x_ave
        self.x_2_ave = x_2_ave

    def show(self):
        """ Show the graph.
        """
        fig = plt.figure('random walk', figsize=(8, 8))

        ax1 = fig.add_subplot(311)
        ax1.plot([n for n in self.N], self.x_ave)
        ax1.set_ylabel(r'$<x(N)>$', fontsize=16)

        ax2 = fig.add_subplot(312)
        ax2.plot([n for n in self.N], self.x_2_ave)
        ax2.set_ylabel(r'$<x^{2}(N)>$', fontsize=16)

        ax3 = fig.add_subplot(313)
        ax3.set_ylabel(r'$<\Delta x^{2}(N)>$', fontsize=16)
        ax3.plot([n for n in self.N], self.x_2_ave - self.x_ave ** 2)
        ax3.set_xlabel(r'$N$')

        plt.show()

    def caluculate_error(self, N):
        """ Caluculate the error of <\Delta x^{2}(N)> and preview.

        N : (int)
        """
        resN_0 = 4. * self.prob * (1. - self.prob) * (self.l ** 2) * N
        _N = range(1, N + 1)

        M = 2
        count = 0
        while count < 15:
            resN = np.zeros(M, 'f')
            for m in range(M):
                self.random_walk_d1(_N)
                t = self.calc_ave()
                resN[m] = self.x_2_ave[N - 1] - self.x_ave[N - 1] ** 2
            std_resN = np.std(resN)

            if M > (std_resN * 100. / resN_0) ** 2:
                print str(M) + " & $>$ & " + str((std_resN / (0.01 * resN_0)) ** 2) \
                    + " & " + str(count + 1) + " \\\\"
                count += 1
            else:
                print str(M) + " & $<$ & " + str((std_resN / (0.01 * resN_0)) ** 2) + " & \\\\"
            M += 1
        return None

if __name__ == '__main__':

    rw1 = RandomWalk1()
    # --- 問題a ---
    N = [4, 8, 16, 32]  # caluculate when N = *
    rw1.random_walk_d1(N)
    rw1.calc_ave()
    rw1.show()

    # --- 問題b ---
# rw1.caluculate_error(8) # 8 or 32
