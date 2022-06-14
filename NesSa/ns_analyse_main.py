#!/usr/bin/env python3

from __future__ import print_function

import sys, math
import numpy as np
from NesSa import ns_analyse

def analyse(inputs=None):
    if inputs==None:
        args=ns_analyse.parse_cli()
    else:
        args=ns_analyse.argparse.Namespace(**inputs)
    
    T_min = args.T_min
    dT = args.dT
    n_T = args.n_T
    k_Boltzmann = args.kB
    debug = args.verbose
    profile = args.profile
    reweight_delta_P = args.delta_pressure

    if args.accurate_sum:
        sum=math.fsum
    else:
        sum=np.sum

    (n_walkers, n_cull, n_Extra_DOF, flat_V_prior, N_atoms, Es, Vs) = ns_analyse.read_inputs(args.files, line_skip=args.skip, line_end=args.line_end, interval=args.interval)

    try:
        from mpi4py import MPI
        comm_rank = MPI.COMM_WORLD.Get_rank()
        comm_size = MPI.COMM_WORLD.Get_size()
        sys.stderr.write("Using mpi %d/%d\n" % (comm_rank, comm_size))
    except:
        comm_rank = 0
        comm_size = 1


    if reweight_delta_P != 0.0:
        Es += reweight_delta_P*Vs

    E_min_i = np.argmin(Es)
    E_min = Es[E_min_i]
    Es -= E_min

    if args.kolmogorov_smirnov_test:
        from scipy import stats
        def ks_test_gaussianity_histogram(histogram):
            dev_max = 0.0
            histo_mean =    np.sum(histogram[0] * (histogram[1][1:]-histogram[1][0:-1]) * (histogram[1][1:]+histogram[1][0:-1])/2.0 )
            histo_2nd_mom = np.sum(histogram[0] * (histogram[1][1:]-histogram[1][0:-1]) * ((histogram[1][1:]+histogram[1][0:-1])/2.0)**2 )
            histo_std_dev = np.sqrt(histo_2nd_mom - histo_mean**2)
            for ibin in range(1,len(histogram[0])):
                numerical_cumul = np.sum(histogram[0][0:ibin] * (histogram[1][1:ibin+1]-histogram[1][0:ibin]))
                analytical_cumul = stats.norm.cdf(histogram[1][ibin], histo_mean, histo_std_dev)
                dev_max = max(dev_max, np.abs( numerical_cumul - analytical_cumul))
            return dev_max

    def analyse(log_a, T, interval=1):
        beta = 1.0/(k_Boltzmann*T)

        (Z_term, shift) = ns_analyse.calc_Z_terms(beta, log_a, Es, flat_V_prior, N_atoms, Vs)

        if debug:
            for i in range(len(Z_term)):
                print("Z_term ", T, Es[i], Z_term[i])
            print("Z_term")
            print("Z_term")

        Z = sum(Z_term)

        U_pot = sum(Z_term*Es)/Z
        U = n_Extra_DOF/(2.0*beta) + U_pot + E_min

        Cvp = n_Extra_DOF*k_Boltzmann/2.0 + k_Boltzmann *beta*beta * (sum(Z_term * Es**2)/Z - U_pot**2)

        if Vs is not None:
            V = sum(Z_term*Vs)/Z
            #thermal_exp = -1.0/V * (sum(Z_term*Vs*Vs)*(-beta)*Z - sum(Z_term*Vs)*sum(Z_term*Vs)*(-beta)) / Z**2
            thermal_exp = -1.0/V * k_Boltzmann *beta*beta * (sum(Z_term*Vs)*sum(Z_term*Es)/Z - sum(Z_term*Vs*Es)) / Z
        else:
            V = None
            thermal_exp = None

        if debug:
            sys.stderr.write("T %f extra DoF U_offset %f C_offset %f remaining Cvp %f\n" % (T, n_Extra_DOF/(2.0*beta), n_Extra_DOF*k_Boltzmann/2.0, Cvp-n_Extra_DOF*k_Boltzmann/2.0))

        log_Z = math.log(Z)+shift-beta*E_min
        # Z(T=0) contributed entirely by lowest energy, since all others are exponentially suppressed by exp(-beta E_i) term
        # Z(T=0) = a(E_min) * exp (- beta*E_min)
        # Helmholtz free energy F = - log_Z / beta
        # F(T=0) = -log_a(E_min)/beta + E_min
        # shift F so that it equals lowest internal energy (really E+PV) of all NS samples
        Helmholtz_F = -log_Z/beta + log_a[E_min_i]/beta

        Z_max = np.amax(Z_term)
        low_percentile_config = np.where(Z_term > Z_max/10.0)[0][0]*interval
        high_percentile_config = np.where(Z_term > Z_max/10.0)[0][-1]*interval
        mode_config = np.argmax(Z_term)*interval

        if args.kolmogorov_smirnov_test and reweight_delta_P != 0.0:
            very_low_percentile_config = np.where(Z_term > Z_max/100.0)[0][0]*interval
            very_high_percentile_config = np.where(Z_term > Z_max/100.0)[0][-1]*interval
            V_histogram = np.histogram(Vs[very_low_percentile_config:very_high_percentile_config], bins=30, weights=Z_term[very_low_percentile_config:very_high_percentile_config], density=True)
            ks_gaussianity = ks_test_gaussianity_histogram(V_histogram)
        else:
            ks_gaussianity = None

        Z_fract = np.sum(Z_term[low_percentile_config:high_percentile_config+1])/Z

        if args.dump_terms >= 0.0:
            n_Z_terms=0
            Z_terms_fract=0.0
            print("#Z_term i H V Z_term")
            for i in range(len(Es)):
                if Z_term[i] > Z_max*args.dump_terms:
                    Z_terms_fract += Z_term[i]
                    n_Z_terms += 1
                    print("Z_term",i, Es[i]+E_min, Vs[i], Z_term[i])
            print("#Z_term fract total ", n_Z_terms, Z_terms_fract/Z, high_percentile_config-low_percentile_config, Z_fract)

        return (log_Z, Helmholtz_F, U, Cvp, low_percentile_config, mode_config, high_percentile_config, Z_fract, V, ks_gaussianity, thermal_exp)


    log_a = ns_analyse.calc_log_a(len(Es), n_walkers, n_cull, interval=args.interval)

    if comm_rank == 0:
        # check this: a = N_w/(N_w+N_c)
        print("# n_walkers ",n_walkers," n_cull ",n_cull," simple log(alpha) ",math.log(float(n_walkers)) - math.log(float(n_walkers+n_cull)))
        n_therm = 4
        if args.entropy:
            n_therm += 1
        fmt = "{:7g} "+" ".join(["{:11g}"]*n_therm)+" "+" ".join(["{:8d}"]*3)+" "+" ".join(["{:15d}"]*3)+" {:6g}"
        header_labels=["T","log_Z(T)","F(T)"]
        if args.entropy:
            header_labels.append("S")
        header_labels.extend(["U(T)","Cvp(T)","low_pct_i(T)","mode_i(T)","high_pct_i(T)","low_pct_iter(T)","mode_iter(T)","high_pct_iter(T)","Z_fract(T)"])
        if Vs is not None:
            fmt += " {:8g} {:14g}"
            header_labels.append("V(T)")
            header_labels.append("thermal_exp(T)")
        if args.kolmogorov_smirnov_test:
            fmt += " {:13g}"
            header_labels.append("KS_Gauss_V(T)")
        header_fmt = "#"+fmt.replace("f","s").replace("g","s").replace("d","s").replace(":",":^")
        print(header_fmt.format(*header_labels))

    data = []
    for i_T in range(comm_rank, n_T, comm_size):
        if not args.quiet:
            if (i_T + 1) % 10 == 0:
                if (i_T + 1) % 100 == 0:
                    sys.stderr.write(f'{(i_T+1)//100 % 10}')
                else:
                    sys.stderr.write('.')
                sys.stderr.flush()
        T = T_min + i_T * dT
        if args.profile and i_T == 0:
            import cProfile
            cProfile.run('(log_Z, Helmholtz_F, U, Cvp, low_pct, mode, high_pct, Z_fract, V, ks_gaussianity) = analyse(log_a, T, args.interval)','ns_analyse.stats')
        else:
            (log_Z, Helmholtz_F, U, Cvp, low_pct, mode, high_pct, Z_fract, V, ks_gaussianity, thermal_exp) = analyse(log_a, T, args.interval)
        if args.entropy:
            S = (U - Helmholtz_F)/(k_Boltzmann * T)
        else:
            S = None
        data.append( (T, log_Z, Helmholtz_F, S, U, Cvp, low_pct, mode, high_pct, int(low_pct/n_cull), int(mode/n_cull), int(high_pct/n_cull), Z_fract, V, ks_gaussianity, thermal_exp) )
    if not args.quiet:
        sys.stderr.write('\n')

    try:
        data = MPI.COMM_WORLD.gather(data, root = 0)
        data = [item for sublist in data for item in sublist]
    except:
        pass

    if comm_rank == 0:
        data = sorted(data, key = lambda x: x[0])
        for row in data:
            print(fmt.format(*[r for r in row if r is not None]))

if __name__ == "__main__":
    analyse()