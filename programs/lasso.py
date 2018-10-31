import numpy as np
from tqdm import tqdm
from sklearn.linear_model import Lasso
from sklearn.metrics import mean_squared_error, r2_score


def energy(s):
    return -(s[0] * s[-1] + np.dot(s[1:], s[:-1]))


num_states = 100
num_bootstraps = 20
L = 40

spins = np.random.choice([-1.0, 1.0], size=(2 * num_states, L))
print(spins)

couplings = np.zeros((2 * num_states, L**2))
energies = np.zeros(2 * num_states)
for i in range(2 * num_states):
    couplings[i, :] = -np.ravel(np.outer(spins[i, :], spins[i, :]))
    energies[i] = energy(spins[i, :])
print(couplings)
print(energies)

test_couplings = couplings[num_states:]
test_energies = energies[num_states:]

couplings = couplings[:num_states]
energies = energies[:num_states]
print(couplings)
print(energies)

lambdas = 10.0**np.arange(-4, 5)
betas = np.zeros((num_bootstraps, L**2))

N_train = int(round(0.8 * num_states))
N_test = num_states - N_train
X_train = couplings[:N_train, :]
y_train = energies[:N_train]
X_test = couplings[N_train:, :]
y_test = energies[N_train:]
predictions = np.zeros((num_bootstraps, N_test))

with open("data/J_LASSO.dat", "w") as f_J, open("data/mse_LASSO.dat",
                                                "w") as f_mse:
    f_J.write("%d\n" % L)
    for l in lambdas:
        fitter = Lasso(alpha=l, fit_intercept=False)
        fitter.fit(couplings, energies)

        J = np.array(fitter.coef_).reshape(L, L)
        f_J.write("\n")
        print(J)
        np.savetxt(f_J, J)

        for i in tqdm(range(num_bootstraps)):
            indices = np.random.randint(0, N_train, N_train)
            X_sel = X_train[indices]
            y_sel = y_train[indices]
            fitter.fit(X_sel, y_sel)
            betas[i, :] = fitter.coef_
            predictions[i, :] = fitter.predict(X_test)

        mean_beta = np.sum(betas, axis=0)
        pred = np.matmul(test_couplings, mean_beta)
        r2 = r2_score(test_energies, pred)
        mse = mean_squared_error(pred, test_energies)

        bias = np.mean((y_test - np.mean(predictions, axis=0, keepdims=True))
                       **2)
        variance = np.mean(
            (predictions - np.mean(predictions, axis=0, keepdims=True))**2)

        f_mse.write("%g %g %g %g %g\n" % (l, r2, mse, bias, variance))
