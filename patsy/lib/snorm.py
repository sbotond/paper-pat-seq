#
# Copyright (C) 2013 EMBL - European Bioinformatics Institute
#
# This program is free software: you can redistribute it
# and/or modify it under the terms of the GNU General
# Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your
# option) any later version.
#
# This program is distributed in the hope that it will be
# useful, but WITHOUT ANY WARRANTY; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A
# PARTICULAR PURPOSE. See the GNU General Public License
# for more details.
#
# Neither the institution name nor the name patsy
# can be used to endorse or promote products derived from
# this software without prior written permission. For
# written permission, please contact <sbotond@ebi.ac.uk>.

# Products derived from this software may not be called
# patsy nor may patsy appear in their
# names without prior written permission of the developers.
# You should have received a copy of the GNU General Public
# License along with this program. If not, see
# <http://www.gnu.org/licenses/>.

import  numpy   as np
import  scipy
from    scipy   import  stats
from    scipy   import  optimize
import      matplotlib                                                                                                                                         
from        matplotlib.backends.backend_pdf import  PdfPages
from        matplotlib                      import  pyplot          as  plt

class FitSkewNormal:
    """ Fit skew normal distribution on the fragment size distribution. """
    def __init__(self, frags,w, log, report):
        self.frags          = frags
        self.log            = log
        self.report         = report
        self.low            = 0
        self.high           = len(frags) - 1
        self.w              = w

    def snorm_pdf(self, x, loc, scale, shape):
        """ Calculate skew normal density. """
        x   = (x - loc)/float(scale)
        return 2.0 / scale * stats.norm.pdf(x) * stats.norm.cdf(shape*x)

    def pdf(self,x):
        return self.snorm_pdf(x,self.loc,self.scale,self.shape)
    
    def init_params(self):
        """ Initialize parameters based on sample mean and sd. """
        fd  = self.frags

        m   = 0.0 
        for i in xrange(self.low, self.high+1):
            m += fd[i] * i

        s   = 0.0
        for i in xrange(self.low, self.high+1):
            s += fd[i] * (i - m) ** 2
        s = np.sqrt(s)
        return np.array([m,s,1.0],dtype=float) 

    def snorm_var(self, loc,scale,shape):
        d   = shape/np.sqrt(1+shape**2)
        var = ( (scale**2) * (1 - 2*(d**2)/np.pi) )
        return var

    def log_lik(self, p, x, d):
        """ Calculate log likelihood """
        pdf = self.snorm_pdf(x, p[0], p[1], p[2])
        i   = pdf > 0.0
        L   = (1-self.w) * np.sum(np.log(pdf[i]) * d[i]) - self.w * self.snorm_var(p[0],p[1],p[2])
        return L
        
    def fit(self):
        """ Fit skew normal distribution. """
        f           = self.frags
        init_params = self.init_params()
        
        x       = np.arange(self.low, self.high+1)
        d       = self.frags

        cost_func           = lambda p: -self.log_lik(p,x,d)
        fitted_params  = optimize.fmin(func=cost_func, x0=init_params, maxiter=10**6, maxfun=10**6, disp=False)
        #fitted_params  = optimize.fmin(func=cost_func, x0=init_params, maxiter=10**6, maxfun=10**6, disp=True)
        #fitted_params  = optimize.fmin_tnc(func=cost_func, approx_grad=True, bounds=[ (0,float('inf')), (0,float('inf')), (float('-inf'),float('inf'))], x0=init_params, disp=0)[0]
        self.L          = self.log_lik(fitted_params, x, d)

        self.loc        = fitted_params[0]
        self.scale      = fitted_params[1]
        self.shape      = fitted_params[2]

        self.log.log("Estimated location parameter: %d" % self.loc)
        self.log.log("Estimated scale parameter: %d" % self.scale)
        self.log.log("Estimated shape parameter: %g" % self.shape)

    def plot_model(self, name):
        """ Plot fitted skew normal distribution. """
        fig                 = plt.figure()
        f                   = self.frags/np.sum(self.frags)
        # Plot empirical distribution:
        plt.bar(np.arange(len(f)),f,width=0.1)
        # Plot fitted curve:
        #y = stats.truncnorm.pdf(np.arange(len(f)+1),self.low, self.high, loc=self.mean, scale=self.sd) 
        y = self.snorm_pdf(x=np.arange(self.high+1), loc=self.loc, scale=self.scale, shape=self.shape) 
        y = y/np.sum(y)
        plt.plot(np.arange(len(y)),y)
        plt.xlim((self.low-100, self.high+100))
        plt.title("%s - Fitted fragment size distribution (skew normal)"  % name) 
        self.report.pages.savefig(fig)

