#!/usr/bin/python3

import numpy as np
DUMMY_DATA_1 = np.asarray([
    np.asarray([-7.4, 6.9, 3.7]),
    np.asarray([-6.2, -3.1, 8.9]),
    np.asarray([-12.7, -6.3, -3.4])
])

DUMMY_DATA_2 = np.asarray([
    np.asarray([-55.5, 42.9, 23.8]),
    np.asarray([-44.8, -21.1, 56.6]),
    np.asarray([-85.7, -42.4, -22.5])
])

XPC = {
    "init_conc": 0.14,
    "on": {
        "damaged": 0.008,
        "partially": 0.002,
        "fully": 0.002,
        "incised": 0.22,
        "resynth": 0.,
        "rechrom": 0.
    },
    "off": {
        "damaged": 0.061,
        "partially": 0.007,
        "fully": 0.007,
        "incised": 0.4,
        "resynth": 0.,
        "rechrom": 0.
    },
    "name": "XPC"
}

TFIIH = {
    "init_conc": 0.36,
    "on": {
        "damaged": 1.6,
        "partially": 0.26,
        "fully": 0.26,
        "incised": 0.0004,
        "resynth": 0.,
        "rechrom": 0.
    },
    "off": {
        "damaged": 0.053,
        "partially": 0.012,
        "fully": 0.012,
        "incised": 0.05,
        "resynth": 0.,
        "rechrom": 0.
    },
    "name": "TFIIH"
}

XPG = {
    "init_conc": 0.44,
    "on": {
        "damaged": 0.,
        "partially": 0.28,
        "fully": 0.28,
        "incised": 0.001,
        "resynth": 0.,
        "rechrom": 0.
    },
    "off": {
        "damaged": 0.,
        "partially": 0.015,
        "fully": 0.015,
        "incised": 0.1,
        "resynth": 0.,
        "rechrom": 0.
    },
    "name": "XPG"
}

XPA = {
    "init_conc": 1.11,
    "on": {
        "damaged": 0.,
        "partially": 0.13,
        "fully": 0.13,
        "incised": 0.004,
        "resynth": 0.054,
        "rechrom": 0.
    },
    "off": {
        "damaged": 0.,
        "partially": 1.04,
        "fully": 1.04,
        "incised": 0.06,
        "resynth": 0.004,
        "rechrom": 0.
    },
    "name": "XPA"
}

XPF = {
    "init_conc": 0.17,
    "on": {
        "damaged": 0.,
        "partially": 1.2,
        "fully": 1.2,
        "incised": 0.09,
        "resynth": 0.,
        "rechrom": 0.
    },
    "off": {
        "damaged": 0.,
        "partially": 0.01,
        "fully": 0.01,
        "incised": 0.05,
        "resynth": 0.,
        "rechrom": 0.
    },
    "name": "XPF"
}

RPA = {
    "init_conc": 1.11,
    "on": {
        "damaged": 0.,
        "partially": 0.15,
        "fully": 0.006,
        "incised": 0.006,
        "resynth": 0.08,
        "rechrom": 0.07
    },
    "off": {
        "damaged": 0.,
        "partially": 2.6,
        "fully": 0.021,
        "incised": 0.021,
        "resynth": 0.04,
        "rechrom": 0.04
    },
    "name": "RPA"
}

PCNA = {
    "init_conc": 1.11,
    "on": {
        "damaged": 0.,
        "partially": 0.,
        "fully": 0.,
        "incised": 0.001,
        "resynth": 0.01,
        "rechrom": 0.31
    },
    "off": {
        "damaged": 0.,
        "partially": 0.,
        "fully": 0.,
        "incised": 0.004,
        "resynth": 0.002,
        "rechrom": 0.05
    },
    "name": "PCNA"
}

DAMAGED = {
    "enzyme_dict": {
        "0000000": {
            "factor": [3.1, 11.],
            "inter": [1, 2]
        },
        "0100000": {
            "factor": [-0.08],
            "inter": [0]
        },
    },
    "mask": [True, True, False, False, False, False, False],
    "name": "damaged"
}

PARTIALLY = {
    "enzyme_dict": {
        "0000000": {
            "factor": [-3.1],
            "inter": [1]
        },
        "1100000": {
            "factor": [11.],
            "inter": [0]
        },
        "1111010": {
            "factor": [-0.74],
            "inter": [1]
        },
        "1111110": {
            "factor": [-0.74],
            "inter": [1]
        }
    },
    "mask": [True, True, True, True, True, True, False],
    "name": "partially"
}

FULLY = {
    "enzyme_dict": {
        "0000000": {
            "factor": [-11.],
            "inter": [2]
        },
        "0111110": {
            "factor": [-4.1],
            "inter": [2]
        },
        "1111010": {
            "factor": [0.74],
            "inter": [1]
        },
        "1111110": {
            "factor": [-4.1, 0.74],
            "inter": [2, 1]
        }
    },
    "mask": [True, True, True, True, True, True, False],
    "name": "fully"
}

INCISED = {
    "enzyme_dict": {
        "1111111": {
            "factor": [4.1],
            "inter": [2]
        }
    },
    "mask": [True, True, True, True, True, True, True],
    "name": "incised"
}

RESYNTH = {
    "enzyme_dict": {
        "0000011": {
            "factor": [0.05, -0.012],
            "inter": [3, 4]
        }
    },
    "mask": [False, False, False, True, False, True, True],
    "name": "resynth"
}

RECHROM = {
    "enzyme_dict": {
        "0000011": {
            "factor": [0.012],
            "inter": [4]
        }
    },
    "mask": [False, False, False, False, False, True, True],
    "name": "rechrom"
}