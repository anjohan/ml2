import numpy as np
import sys
from tqdm import tqdm
from sklearn.linear_model import Lasso
from sklearn.metrics import mean_squared_error, r2_score


def energy(s):
    return -(s[0] * s[-1] + np.dot(s[1:], s[:-1]))


num_states = int(sys.argv[1])
num_bootstraps = 20

L = 40

# create 10000 random Ising states
states = np.random.choice([-1, 1], size=(2 * num_states, L))


def ising_energies(states, L):
    """
    This function calculates the energies of the states in the nn Ising Hamiltonian
    """
    J = np.zeros((L, L))
    for i in range(L):
        J[i, (i + 1) % L] -= 1.0
    # compute energies
    E = np.einsum("...i,ij,...j->...", states, J, states)

    return E


# calculate Ising energies
energies = ising_energies(states, L)

states = -np.einsum("...i,...j->...ij", states, states)
shape = states.shape
states = states.reshape((shape[0], shape[1] * shape[2]))

X_train = states[:num_states]
y_train = energies[:num_states]
X_test = states[num_states:]
y_test = energies[num_states:]


lambdas = 10.0 ** np.arange(-4, 5)
betas = np.zeros((num_bootstraps, L ** 2))
predictions = np.zeros((num_bootstraps, num_states))
mses = np.zeros(num_bootstraps)

with open("data/J_LASSO_" + str(num_states) + ".dat", "w") as f_J, open(
    "data/mse_LASSO_" + str(num_states) + ".dat", "w"
) as f_mse:
    f_J.write("%d\n" % L)
    for l in lambdas:
        fitter = Lasso(alpha=l, fit_intercept=False)
        print(l)
        fitter.fit(X_train, y_train)

        J = np.array(fitter.coef_).reshape(L, L)
        f_J.write("\n")
        print(J)
        np.savetxt(f_J, J)
        r2 = fitter.score(X_test, y_test)
        pred = fitter.predict(X_test)
        mse = mean_squared_error(pred, y_test)

        for i in tqdm(range(num_bootstraps)):
            indices = np.random.randint(0, num_states, num_states)
            X_sel = X_train[indices]
            y_sel = y_train[indices]
            fitter.fit(X_sel, y_sel)
            betas[i, :] = fitter.coef_
            predictions[i, :] = fitter.predict(X_test)
            mses[i] = mean_squared_error(predictions[i, :], y_test)

        mse = np.mean(mses)
        bias = np.mean((y_test - np.mean(predictions, axis=0, keepdims=True)) ** 2)
        variance = np.mean(
            (predictions - np.mean(predictions, axis=0, keepdims=True)) ** 2
        )

        f_mse.write("%g %g %g %g %g\n" % (l, r2, mse, bias, variance))
