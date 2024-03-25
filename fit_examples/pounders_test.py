import numpy as np
from pounders import pounders_mpi_wrapper, stat_analysis

class linear_model(object):

    r"""
    Linear model y = a*x + b
    """

    def __init__(self):
        # initialize data and attributes
        coeffs = [0.6, 0.4]
        n_obs = 10
        n_params = len(coeffs)
        np.random.seed(123)
        x = np.random.rand(n_obs)
        error = np.random.normal(0.0, 0.1, n_obs)
        y = coeffs[0]*x + coeffs[1] + error

        self.n_obs = n_obs
        self.n_params = n_params
        self.x = x
        self.y = y

    def formResidual(self, X):
        a, b = X
        return (a * self.x + b) - self.y

# A test example for POUNDerS.
if __name__ == '__main__':
    model = linear_model()
    max_iter = None
    sol, summary = pounders_mpi_wrapper(model.formResidual, (0.5,0.5), model.n_obs, max_iter=max_iter)
    print('POUNDerS', sol, summary)
    stat_result = stat_analysis(model.formResidual, sol, 1e-3, summary['residual'])
    print('stat_analysis', *stat_result)

    # Compare with numpy.polyfit. 
    p, v = np.polyfit(model.x, model.y, 1, cov=True)
    print('Numpy Polyfit', p, v)