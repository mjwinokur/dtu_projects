import hoomd, hoomd.md
hoomd.context.initialize()
hoomd.init.read_gsd('init.gsd')
nl = hoomd.md.nlist.cell()
lj = hoomd.md.pair.lj(r_cut=2.5, nlist=nl)
lj.pair_coeff.set('A', 'A', epsilon=1.0, sigma=1.0)
hoomd.md.integrate.mode_standard(dt=0.005)
all = hoomd.group.all()
hoomd.md.integrate.langevin(group=all, kT=1.2, seed=4)
hoomd.run(1e5)
