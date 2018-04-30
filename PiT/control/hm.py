# Solver for damped harmonic motion

import numpy as np
import scipy.integrate
import yaml

from interval import Interval

class Integrator(object):
    """Our system is:
    d^2x / dt^2 + b dx/dt + omega^2 x = 0

    Transform to 2 coupled equations:

    0: dx/dt = v
    1: dv/dt = -b v - omega^2 x

    Define y = (x, v)
    """
    
    def __init__(self, params):
        self.b = params['b']
        self.omega = params['omega']
        # The system coefficients
        self.mat = np.array([
            [0.0, 1.0],
            [-self.omega**2, -self.b]], dtype=float)

        return
        
    def rhs(self, y, t):
        """Evaluates the rhs of the system."""
        
        return np.matmul(self.mat, y)

    def jac(self, y, t):
        return self.mat

    pass

class CoarseInt(Integrator):
    def __call__(self, y0, iv):
        ts = iv.generate()
        vals, info = scipy.integrate.odeint(
            self.rhs,
            y0,
            ts,
            rtol=1e-2,
            atol=1e-3,
            full_output=True)

        return (ts, vals)
    pass
class FineInt(Integrator):
    def __call__(self, y0, iv):
        ts = iv.generate()
        vals, info = scipy.integrate.odeint(
            self.rhs,
            y0,
            ts,
            Dfun=self.jac,
            rtol=1e-4,
            atol=1e-6,
            full_output=True)
        return (ts, vals)
    pass

class Problem(object):
    def __init__(self, params, ic, iv, solver):
        self.yi = np.array((ic['x'], ic['v']), dtype=float)
        self.iv = iv

        IntCls = {'fine': FineInt, 'coarse': CoarseInt}[solver]
        self.integrator = IntCls(params)
        return
    
    @classmethod
    def from_file(cls, fn):
        with open(fn) as f:
            return cls.from_stream(f)
        
    @classmethod
    def from_stream(cls, stream):
        state = yaml.load(stream)
        return cls.from_dict(state)
    
    @classmethod
    def from_dict(cls, state):
        params = state['params']
        time = state['time']
        iv = Interval(time['dt'], time['start'], time['n'])
        return cls(
            state['params'],
            state['initial'], iv,
            state['solver'])
    
    def run(self, outfile):
        ts, vals = self.integrator(self.yi, self.iv)
        data = np.concatenate((ts[:, np.newaxis], vals), axis=1)
        np.savetxt(outfile, data)
        
if __name__ == '__main__':
    import sys
    prob = Problem.from_file(sys.argv[1])
    prob.run(sys.argv[2])
    
    
