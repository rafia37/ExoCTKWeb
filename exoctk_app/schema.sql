CREATE TABLE "ldc" (
    'id'	INTEGER NOT NULL UNIQUE,
    'date'  TEXT NOT NULL,
    'n_bins'    INTEGER,
    'teff'  REAL,
    'logg' REAL,
    'feh'   REAL,
    'bandpass'  TEXT,
    'modeldir'  TEXT,
    'wave_min'  REAL,
    'mu_min'    REAL,
    'wave_max'  REAL,
    'local_files'   TEXT,
    'pixels_per_bin'    INTEGER,
    'uniform'   TEXT,
    'linear'    TEXT,
    'quadratic'    TEXT,
    'squareroot'    TEXT,
    'logarithmic'    TEXT,
    'exponential'    TEXT,
    'three_parameter'    TEXT,
    'four_parameter'    TEXT,
    PRIMARY KEY(id)
);

CREATE TABLE "tor" (
    'id'    INTEGER NOT NULL UNIQUE,
    'date'  TEXT NOT NULL,
    'ins'   TEXT,
    'mag'   REAL,
    'groups'    TEXT,
    'amps'  INTEGER,
    'subarray'  TEXT,
    'sat_lvl'   REAL,
    'sat'   TEXT,
    'T' REAL,
    'n_reset'   INTEGER,
    'band'  TEXT,
    'filt'  TEXT,
    PRIMARY KEY(id)
);

CREATE TABLE "exotransmit" (
    'id'    INTEGER NOT NULL UNIQUE,
    'date'  TEXT NOT NULL,
    'eos'   TEXT,
    'tp'    TEXT,
    'g' REAL,
    'R_p'   REAL,
    'R_s'   REAL,
    'P' REAL,
    'Rayleigh'  REAL,
    PRIMARY KEY(id)
);