#!/usr/bin/python3
import sys
import pymatgen
import seekpath
import numpy
from pymatgen.core.periodic_table import get_el_sp
import json


def write_atom(f, avec, typ, nat, pos, atom):
    print("CELL_PARAMETERS angstrom", file=f)
    for ii in range(3):
        print(" %f %f %f" % (avec[ii, 0], avec[ii, 1], avec[ii, 2]), file=f)
    print("ATOMIC_SPECIES", file=f)
    for ityp in typ:
        print(" %s %f %s" % (ityp, pymatgen.Element(ityp).atomic_mass, sssp[str(ityp)]["filename"]), file=f)
    print("ATOMIC_POSITIONS crystal", file=f)
    for iat in range(nat):
        print(" %s %f %f %f" % (
            atom[iat], pos[iat][0], pos[iat][1], pos[iat][2]), file=f)


def write_middle(f, nat, ntyp, ecutwfc, ecutrho, psdir):
    print("  pseudo_dir = \'%s/\'" % psdir, file=f)
    print("/", file=f)
    print("&SYSTEM", file=f)
    print("       ibrav = 0", file=f)
    print("         nat = %d" % nat, file=f)
    print("        ntyp = %d" % ntyp, file=f)
    print("     ecutwfc = %f" % ecutwfc, file=f)
    print("     ecutrho = %f" % ecutrho, file=f)


if __name__ == '__main__':

    args = sys.argv
    if len(args) < 2:
        print("Usage:")
        print("$ cif2qe.py cif-file pseudo-directory [dk_path] [dk_grid]")
        print("Default:")
        print("$ cif2qe.py cif-file pseudo-directory 0.1 0.2")
        exit(0)
    #
    # CIF parser
    #
    structure = pymatgen.Structure.from_file(args[1])
    #
    # Default value
    #
    dk_path = 0.1
    dk_grid = 0.2
    #
    psdir = args[2]
    f = open(psdir + '/sssp_efficiency.json', 'r')
    sssp = json.load(f)
    #
    if len(args) > 3:
        dk_path = float(args[3])
        if len(args) > 4:
            dk_grid = float(args[4])
    #
    print("  dk for band : {0}".format(dk_path))
    print("  dk for grid : {0}".format(dk_grid))

    structure.remove_oxidation_states()

    #
    # Band path and primitive lattice
    #
    skp = seekpath.get_explicit_k_path((structure.lattice.matrix, structure.frac_coords,
                                        [pymatgen.Element(str(spc)).number for spc in structure.species]),
                                       reference_distance=dk_path)
    #
    # Lattice information
    #
    avec = skp["primitive_lattice"]
    bvec = skp["reciprocal_primitive_lattice"]
    pos = skp["primitive_positions"]
    nat = len(skp["primitive_types"])
    atom = [str(get_el_sp(iat)) for iat in skp["primitive_types"]]
    typ = set(atom)
    ntyp = len(typ)
    typ = set(atom)
    #
    # WFC and Rho cutoff
    #
    ecutwfc = 0.0
    ecutrho = 0.0
    for ityp in typ:
        if ecutwfc < sssp[str(ityp)]["cutoff"]:
            ecutwfc = sssp[str(ityp)]["cutoff"]
        if ecutrho < sssp[str(ityp)]["cutoff"] * sssp[str(ityp)]['dual']:
            ecutrho = sssp[str(ityp)]["cutoff"] * sssp[str(ityp)]['dual']
    #
    # k and q grid
    #
    nk = numpy.zeros(3, numpy.int_)
    for ii in range(3):
        norm = numpy.sqrt(numpy.dot(bvec[ii][:], bvec[ii][:]))
        nk[ii] = round(norm / dk_grid)
        print(norm)
    print("SCF k-grid : ", nk[0], nk[1], nk[2])
    #
    # Band path
    #
    print("Band path")
    for ipath in range(len(skp["path"])):
        start = skp["explicit_segments"][ipath][0]
        final = skp["explicit_segments"][ipath][1] - 1
        print("%5d %8s %10.5f %10.5f %10.5f %8s %10.5f %10.5f %10.5f" % (
            final - start + 1,
            skp["explicit_kpoints_labels"][start],
            skp["explicit_kpoints_rel"][start][0],
            skp["explicit_kpoints_rel"][start][1],
            skp["explicit_kpoints_rel"][start][2],
            skp["explicit_kpoints_labels"][final],
            skp["explicit_kpoints_rel"][final][0],
            skp["explicit_kpoints_rel"][final][1],
            skp["explicit_kpoints_rel"][final][2]))
    #
    # Number of electrons
    #
    nbnd = len(atom)*10
    #
    # scf.in : Charge density
    #
    with open("scf.in", 'w') as f:
        print("&CONTROL", file=f)
        print(" calculation = \'scf\'", file=f)
        write_middle(f, nat, ntyp, ecutwfc, ecutrho, psdir)
        print(" occupations = \'tetrahedra_opt\'", file=f)
        print("/", file=f)
        #
        print("&ELECTRONS", file=f)
        print(" mixing_beta = 0.3", file=f)
        print(" conv_thr = %e" % (float(nat)*1.0e-10), file=f)
        print("/", file=f)
        write_atom(f, avec, typ, nat, pos, atom)
        print("K_POINTS automatic", file=f)
        print(" %d %d %d 0 0 0" % (nk[0], nk[1], nk[2]), file=f)
    #
    # nscf.in : Dense k grid
    #
    with open("nscf.in", 'w') as f:
        print("&CONTROL", file=f)
        print(" calculation = \'nscf\'", file=f)
        write_middle(f, nat, ntyp, ecutwfc, ecutrho, psdir)
        print(" occupations = \'tetrahedra_opt\'", file=f)
        print("        nbnd = %d" % nbnd, file=f)
        print("        la2f = .true.", file=f)
        print("/", file=f)
        print("&ELECTRONS", file=f)
        print("/", file=f)
        write_atom(f, avec, typ, nat, pos, atom)
        print("K_POINTS automatic", file=f)
        print(" %d %d %d 0 0 0" % (nk[0]*2, nk[1]*2, nk[2]*2), file=f)
    #
    # band.in : Plot band
    #
    with open("band.in", 'w') as f:
        print("&CONTROL", file=f)
        print(" calculation = \'bands\'", file=f)
        write_middle(f, nat, ntyp, ecutwfc, ecutrho, psdir)
        print("        nbnd = %d" % nbnd, file=f)
        print("/", file=f)
        print("&ELECTRONS", file=f)
        print("/", file=f)
        write_atom(f, avec, typ, nat, pos, atom)
        print("K_POINTS crystal", file=f)
        print(len(skp["explicit_kpoints_rel"]), file=f)
        for ik in range(len(skp["explicit_kpoints_rel"])):
            print(" %f %f %f 1.0" % (
                skp["explicit_kpoints_rel"][ik][0],
                skp["explicit_kpoints_rel"][ik][1],
                skp["explicit_kpoints_rel"][ik][2]),
                    file=f)
    #
    # bands.in : Read by bands.x
    #
    with open("bands.in", 'w') as f:
        print("&BANDS", file=f)
        print("       lsym = .false.", file=f)
        print("/", file=f)
    #
    # proj.in : Read by projwfc.x
    #
    with open("proj.in", 'w') as f:
        print("&PROJWFC", file=f)
        print("      emin = ", file=f)
        print("      emax = ", file=f)
        print("    deltae = 0.1", file=f)
        print("/", file=f)
    #
    # pp.in : Plot Kohn-Sham orbitals
    #
    with open("pp.in", 'w') as f:
        print("&INPUTPP ", file=f)
        print(" plot_num = 7", file=f)
        print("   kpoint = 1", file=f)
        print(" kband(1) = %d" % 1, file=f)
        print(" kband(2) = %d" % nbnd, file=f)
        print("    lsign = .true.", file=f)
        print("/", file=f)
        print("&PLOT  ", file=f)
        print("         iflag = 3", file=f)
        print(" output_format = 5", file=f)
        print("       fileout = \".xsf\"", file=f)
        print("/", file=f)
    #
    # band.gp : Gnuplot script
    #
    with open("band.gp", 'w') as f:
        print("set terminal pdf color enhanced \\", file=f)
        print("dashed dl 0.5 size 8.0cm, 6.0cm", file=f)
        print("set output \"band.pdf\"", file=f)
        print("#", file=f)
        n_sym_points = 1
        final = 0
        x0 = numpy.linalg.norm(avec[0, :]) * 0.5 / numpy.pi
        print("x%d = %f" % (n_sym_points, x0*skp["explicit_kpoints_linearcoord"][final]), file=f)
        for ipath in range(len(skp["path"])):
            start = skp["explicit_segments"][ipath][0]
            if start != final:
                n_sym_points += 1
                print("x%d = %f" % (n_sym_points, x0*skp["explicit_kpoints_linearcoord"][start]), file=f)
            n_sym_points += 1
            final = skp["explicit_segments"][ipath][1] - 1
            print("x%d = %f" % (n_sym_points, x0*skp["explicit_kpoints_linearcoord"][final]), file=f)
        print("#", file=f)
        print("set border lw 2", file=f)
        print("#", file=f)
        print("set style line 1 lt 1 lw 2 lc 0 dashtype 2", file=f)
        print("set style line 2 lt 1 lw 2 lc 0", file=f)
        print("set style line 3 lt 1 lw 1 lc 1", file=f)
        print("#", file=f)
        print("set ytics scale 3.0, -0.5 font \'Cmr10,18\'", file=f)
        print("set xtics( \\", file=f)
        n_sym_points = 1
        final = 0
        label_f = skp["explicit_kpoints_labels"][final]
        if label_f == "GAMMA":
            label_f = "\\241"
        print("\"%s\" x%d" % (label_f, n_sym_points), end="", file=f)
        for ipath in range(len(skp["path"])):
            start = skp["explicit_segments"][ipath][0]
            label_s = skp["explicit_kpoints_labels"][start]
            if label_s == "GAMMA":
                label_s = "\\241"
            label_f = skp["explicit_kpoints_labels"][final]
            if label_f == "GAMMA":
                label_f = "\\241"
            if start != final:
                n_sym_points += 1
                print(", \\\n\"%s%s\" x%d" % (label_f, label_s, n_sym_points), end="", file=f)
            n_sym_points += 1
            final = skp["explicit_segments"][ipath][1] - 1
            label_f = skp["explicit_kpoints_labels"][final]
            if label_f == "GAMMA":
                label_f = "\\241"
            print(", \\\n\"%s\" x%d" % (label_f, n_sym_points), end="", file=f)
        print(") \\\noffset 0.0, 0.0 font \'Cmr10,18\'", file=f)
        print("#", file=f)
        print("set grid xtics ls 2 front", file=f)
        print("#", file=f)
        print("unset key", file=f)
        print("#", file=f)
        print("set xzeroaxis ls 1", file=f)
        print("#", file=f)
        print("set ylabel \"Energy from {/Cmmi10 \\042}_F [eV]\" offset - 0.5, 0.0 font \'Cmr10,18\'", file=f)
        print("#", file=f)
        n_sym_points = 1
        final = 0
        for ipath in range(len(skp["path"])):
            start = skp["explicit_segments"][ipath][0]
            if start == final:
                n_sym_points += 1
                final = skp["explicit_segments"][ipath][1] - 1
            else:
                break
        print("plot[:][emin:emax] \"bands.out.gnu\" u 1:($2-ef) w l ls 3", file=f)
