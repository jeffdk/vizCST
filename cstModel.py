import consts
import numpy
from fortranFile import FortranFile
import matplotlib.pyplot as plt


class cstModel(object):

    # dictionary variables defined in init
    # See CST paper for more details on variables:
    #  Cook, G. B., Shapiro, S. L., & Teukolsky, S. A.
    #  1992, Astrophys. J., 398, 203
    dataDict = None

    #rs are the values of the radial grid-points
    rs = None
    NR = None
    NSurf = None
    #mus are the cos(\thetas) of the polar angle grid-points
    mus = None
    NU = None

    # CST output is in polytropic 'scaled' units with
    polyK = 617.714470405639
    polyGamma = 2.0
    polyN = 1.0 / (polyGamma - 1.0)

    def __init__(self, filename):
        f = FortranFile(filename)

        # Line 1:
        # NSI, NUI, NSURF
        NS, NU, NSurf = f.readArray(numpy.int32)

        # Line 2:
        # rpoe, RADEQUAT, COMEGA
        rpoe, RADEQUAT, COMEGA = f.readArray(numpy.float64)
        # Line 3:
        # Etot, BaryM, AngMom, snp, ROTFUNC
        Etot, BaryM, AngMom, snp, ROTFUNC = f.readArray(numpy.float64)

        # Data
        #  rho, gamma, alpha, omega, ed, angv
        data = f.readArray(numpy.float64)
        data = data.reshape((6, NS + 1, NU + 1))

        self.dataDict = {'rho': data[0],     # CST metric variable \rho
                         'gamma': data[1],   # CST metric variable \gamma
                         'alpha': data[2],   # CST metric variable \alpha
                         'omega': data[3],   # CST metric variable \omega
                         'ed': data[4],      # matter energy density
                         'angv': data[5]}    # matter angular velocity

        # s is the compactified CST coordinate spaced uniformly in 0 - 1
        # transform to get rs
        # note, arange 0.0 to 1.0 omits the point at 1.0
        self.rs = numpy.array([RADEQUAT * s / (1.0 - s)
                               for s in numpy.arange(0.0, 1.00, 1.0 / NS)])

        #scale rs to proper units: kilometers
        self.rs *= self.polyK ** (self.polyN / 2.0) * consts.AGEO_LENGTH_IN_CM
        #scale ed to proper units: g/cm^3
        self.dataDict['ed'] *= self.polyK ** (-self.polyN) * consts.AGEO_DENSITY_IN_CGS

        self.mus = numpy.array([float(j) / float(NU) for j in range(0, NU + 1)])

        #Last point in rs is at infinity so we must trim it off
        # (rs point at infinity is omitted by omitting s = 1.0 above)
        for key in self.dataDict.keys():
            self.dataDict[key] = self.dataDict[key][:-1, :]
            self.dataDict[key] = self.dataDict[key][:-(NSurf), :]
        self.rs = self.rs[:-(NSurf)]
        self.NR = len(self.rs)
        self.NU = len(self.mus)

    def rawPlot(self, var, title=None, rmax=25.0, func=lambda x:x*x):
        rs, mus = numpy.meshgrid(self.rs, self.mus)
        ys = rs * mus
        xs = rs * numpy.sqrt(1.0 - mus * mus)

        print xs
        print ys
        plt.xlim([0, rmax])
        plt.ylim([0, rmax])
        #plt.pcolor(xs, ys, xs * func(self.dataDict[var].transpose()))  # , vmin=10, vmax=14.5)
        plt.pcolor(xs, ys, func(self.dataDict[var].transpose()), vmin=10, vmax=14.5)
        plt.colorbar()
        plt.title(title)
        plt.show()

    def integralOverStar(self, vars, func, edMin):
        """
        Does a volume integral over the star integrating over r, until
        reaching edMin, then integrating over theta.
        Integrates func(vars[0], vars[1], ...)
        Therefore func must be a function of len(vars) variables
        """

        sum = 0.0
        previousIntegralOverR = 0.0
        for muInd in range(0, self.NU):
            currentIntegralOverR = 0.0
            for rInd in range(0, self.NR - 1):
                if self.dataDict['ed'][rInd, muInd] < edMin:
                    #print "RSURF: ", self.rs[rInd]
                    break
                r_i = self.rs[rInd]
                r_iplus1 = self.rs[rInd + 1]
                drCubed = r_iplus1 ** 3 - r_i ** 3
                varsValues_i = [self.dataDict[var][rInd, muInd] for var in vars]
                varsValues_iplus1 = [self.dataDict[var][rInd + 1, muInd] for var in vars]
                averageValue = 0.5 *(func(*varsValues_iplus1) + func(*varsValues_i))
                currentIntegralOverR += drCubed * averageValue
            if muInd > 0:
                dmu = self.mus[muInd] - self.mus[muInd - 1]
                sum += dmu * 0.5 * (currentIntegralOverR + previousIntegralOverR)
            previousIntegralOverR = currentIntegralOverR

        return sum * 4.0 / 3.0 * numpy.pi
