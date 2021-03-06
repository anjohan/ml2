import sys
import numpy as np
import pickle


with open("data/states.pkl", "rb") as f:
    data = pickle.load(f)
with open("data/labels.pkl", "rb") as f:
    labels = pickle.load(f)

data = np.unpackbits(data).reshape(-1, 1600).astype("int32")
data[data == 0] = -1
labels = labels.astype("int32")

skip = 1
if len(sys.argv) == 2:
    skip = int(sys.argv[1])
num_ordered = 70000 / skip
num_critical = 30000 / skip
num_disordered = 60000 / skip
num_spins = 1600

print(np.sum(data[:10, :] == -1, axis=1))

with open("data/states_%d.bin" % skip, "wb") as f:
    np.array([num_spins, num_ordered, num_critical, num_disordered]).astype(
        "int32"
    ).tofile(f)
    data[::skip].astype(np.float64).tofile(f)
    # f.write("%d %d %d %d\n" % (num_spins, num_ordered, num_critical, num_disordered))
    # np.savetxt(f, data, fmt="%d")


with open("data/labels_%d.bin" % skip, "wb") as f:
    labels[::skip].tofile(f)
