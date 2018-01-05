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
    
);