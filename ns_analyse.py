#!/usr/bin/env python3

#ns_analyse code from pymatnest

from __future__ import print_function
import numpy as np, fileinput, itertools, sys, argparse

def parse_cli():
    p = argparse.ArgumentParser()
    p.add_argument('-M','--T_min', action='store', help="""Minimum temperature""",type=float,required=True)
    p.add_argument('-D','--dT', action='store', help="""Temperature step""",type=float,required=True)
    p.add_argument('-n','--n_T', action='store', help="""Number of temperatures""",type=int,required=True)
    p.add_argument('-k','--kB', action='store', help="""Boltzmann constant (defaults to eV/K)""", type=float, default=8.6173324e-5)
    p.add_argument('-a','--accurate_sum', action='store_true', help="""use more accurate sum (math.fsum)""")
    p.add_argument('-v','--verbose', action='store_true', help="""Verbose output (for debugging)""")
    p.add_argument('-p','--profile', action='store_true', help="""do profiling""")
    p.add_argument('-s','--skip', action='store', help="""number of lines to skip""", type=int, default=0)
    p.add_argument('-l','--line_end', action='store', help="""line to ened analysis (python std. zero based last line + 1)""", type=int, default=None)
    p.add_argument('-i','--interval', action='store', help="""interval between lines to use""", type=int, default=1)
    p.add_argument('-P','--delta_pressure', action='store', help="""delta pressure to use for reweighting with flat V prior""", type=float, default=0.0)
    p.add_argument('-V','--dump_terms', action='store', help="""dump weight of each sample""", type=float, default=-1.0)
    p.add_argument('-S','--entropy', action='store_true', help="""compute and print entropy""")
    p.add_argument('-K','--kolmogorov_smirnov_test', action='store_true', help="""Do Kolmogorov-Smirnov test for Gaussianity of distribution of Vs (when --delta_pressure is set)""")
    p.add_argument('-q','--quiet', action='store_true', help="""No progress output""")
    p.add_argument('files', help="""input energies files""")

    args = p.parse_args()
    return args

def read_inputs(args, line_skip=0, line_end=None, interval=1):

    inputs = fileinput.input(files=args)

    fields = inputs.readline().split()
    if len(fields) == 3:
        (n_walkers, n_cull, n_Extra_DOF) = fields
        flat_V_prior = "False"
        N_atoms = "0"
    elif len(fields) == 5:
        (n_walkers, n_cull, n_Extra_DOF, flat_V_prior, N_atoms) = fields
    else:
        raise ValueError("unknown number of fields %d (not 3 or 5) in first line of energies file" % len(fields))
    n_walkers = int(n_walkers)
    n_cull = int(n_cull)
    n_Extra_DOF = int(n_Extra_DOF)
    if flat_V_prior.lower() == "t" or flat_V_prior.lower() == "true":
        flat_V_prior = True
    elif flat_V_prior.lower() == "f" or flat_V_prior.lower() == "false":
        flat_V_prior = False
    else:
        sys.stderr.write("Unknown flat_V_prior string '%s'\n" % flat_V_prior)
        sys.exit(1)
    N_atoms = int(N_atoms)

    Es=[]
    Vs=[]
    lines = itertools.islice(inputs, line_skip, line_end, interval)
    n_fields = None
    for line in lines:
        fields = line.split()
        if n_fields is not None and n_fields != len(fields):
            sys.stderr.write(f'Mismatch field # prev {n_fields} cur {len(fields)}, skipping\n')
            continue
        if len(fields) == 3:
            try:
                E = float(fields[1])
                V = float(fields[2])
            except:
                continue
            if n_fields is None:
                n_fields = 3
            Es.append(E)
            Vs.append(V)
        elif len(fields) == 2:
            try:
                E = float(fields[1])
            except:
                continue
            if n_fields is None:
                n_fields = 2
            Es.append(E)
        else: # silently skip lines with problems
            sys.stderr.write("WARNING: input line with problem: number of fields not 2 or 3, or not floats\n")
            continue

    if len(Vs) == 0:
        Vs = None
    else:
        Vs = np.array(Vs)
    return (n_walkers, n_cull, n_Extra_DOF, flat_V_prior, N_atoms, np.array(Es), Vs)

def calc_log_a(n_Es, n_walkers, n_cull, interval=1):
    # log_a = math.log(float(n_walkers)) - math.log(float(n_walkers+n_cull))
    # From SENS paper PRX v. 4 p 031034 (2014) Eq. 3
    i_range = np.array(range(n_Es*interval))
    i_range_mod_n_cull = np.mod(i_range,n_cull)
    i_range_plus_1_mod_n_cull = np.mod(i_range+1,n_cull)
    # X_n = \prod_{i=0}^n \frac{N-i\%P}{N+1-i\%P}
    # \log X_n = \sum_{i=0}^n \log (N-i\%P) - \log(N+1-i\%P)
    log_X_n_term = np.log(n_walkers-i_range_mod_n_cull) - np.log(n_walkers+1-i_range_mod_n_cull)
    log_X_n = np.cumsum(log_X_n_term)
    # a_n = X_n - X_(n+1)
    # a_n = \prod_{i=0}^n \frac{N-i\%P}{N+1-i\%P} \left(1 - \frac{N-(n+1)\%P}{N+1-(n+1)\%P}\right)
    # \log(a_n) & = \sum_{i=0}^n \left[\log(N-i\%P)-\log(N+1-i\%P)\right] + \log\left(1-\frac{N-(n+1)\%P}{N+1-(n+1)\%P}\right) \\
    #    & = \log(X_n) + \log\left(1-\frac{N-(n+1)\%P}{N+1-(n+1)\%P}\right) \\
    #    & = \log(X_n) + \log\left( \frac{N+1-(n+1)\%P - N + (n+1)\%P}{N+1-(n+1)\%P} \right) \\
    #    & = \log(X_n) + \log \left(  \frac{1}{N+1-(n+1)\%P} \right) \\
    #    & = \log(X_n) -\log \left(  N+1-(n+1)\%P \right)
    log_a = log_X_n[0::interval] - np.log(n_walkers+1-i_range_plus_1_mod_n_cull[0::interval])
    return log_a

def calc_Z_terms(beta, log_a, Es, flat_V_prior=False, N_atoms=0, Vs=None):
    #DEBUG for i in range(len(log_a)):
        #DEBUG print "calc_Z_terms log_a ", log_a[i], beta*Es[i]
    log_Z_term = log_a[:] - beta*Es[:]
    if flat_V_prior:
        log_Z_term += float(N_atoms)*np.log(Vs[:])
    shift = np.amax(log_Z_term[:])
    Z_term = np.exp(log_Z_term[:] - shift)
    return (Z_term, shift)
