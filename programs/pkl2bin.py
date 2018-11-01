import numpy as np
import pickle


def write_arrays(name, states, labels):
    num_states = labels.shape[0]
    num_spins = states.shape[1]
    with open("data/" + name + "_labels.bin", "wb") as f:
        np.array([num_states]).astype("int32").tofile(f)
        labels.astype("int32").tofile(f)

    with open("data/" + name + "_states.bin", "wb") as f:
        np.array([num_states, num_spins]).astype("int32").tofile(f)
        states.astype("int32").tofile(f)


with open("data/states.pkl", "rb") as f:
    data = pickle.load(f)
with open("data/labels.pkl", "rb") as f:
    labels = pickle.load(f)

data = np.unpackbits(data).reshape(-1, 1600).astype("int32")
data[data == 0] = -1
labels = labels.astype("int32")


X_ordered = data[:70000, :]
Y_ordered = labels[:70000]

X_critical = data[70000:100000, :]
Y_critical = labels[70000:100000]

X_disordered = data[100000:, :]
Y_disordered = labels[100000:]

write_arrays("ordered", X_ordered, Y_ordered)
write_arrays("critical", X_critical, Y_critical)
write_arrays("disordered", X_disordered, Y_disordered)
