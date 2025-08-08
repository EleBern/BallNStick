from neuron import h
from neuron.units import ms, mV
import matplotlib.pyplot as plt
import numpy as np
h.load_file('stdrun.hoc')

class BallAndStick:
    def __init__(self, gid):
        self._gid = gid
        self._setup_morphology()
        self._setup_biophysics()
    def _setup_morphology(self):
        self.soma = h.Section(name='soma', cell=self)
        self.dend = h.Section(name='dend', cell=self)
        self.all = [self.soma, self.dend]
        self.dend.connect(self.soma)
        self.soma.L = self.soma.diam = 12.6157
        self.dend.L = 200
        self.dend.diam = 1
    def _setup_biophysics(self):
        for sec in self.all:
            sec.Ra = 100    # Axial resistance in Ohm * cm
            sec.cm = 1      # Membrane capacitance in micro Farads / cm^2
        self.soma.insert('hh')                                          
        for seg in self.soma:
            seg.hh.gnabar = 0.12  # Sodium conductance in S/cm2
            seg.hh.gkbar = 0.036  # Potassium conductance in S/cm2
            seg.hh.gl = 0.0003    # Leak conductance in S/cm2
            seg.hh.el = -54.3     # Reversal potential in mV
        # Insert passive current in the dendrite                      
        self.dend.insert('pas')                                      
        for seg in self.dend:                                          
            seg.pas.g = 0.001  # Passive conductance in S/cm2        
            seg.pas.e = -65    # Leak reversal potential mV          
    def __repr__(self):
        return 'BallAndStick[{}]'.format(self._gid)
        
tstop = 500
factor = 10
text = np.arange(0,1000,100)
vexts = np.ones(len(text))
vexts[::2] = -1
vexts = vexts * factor
vextd = 0.5 * vexts
vexts = h.Vector(vexts)
vextd = h.Vector(vextd)
text = h.Vector(text)

simple_cell = BallAndStick(0)
simple_cell.soma.insert('extracellular')
simple_cell.dend.insert('extracellular')
vexts.play(simple_cell.soma(0.5)._ref_e_extracellular, text)
vextd.play(simple_cell.dend(0.5)._ref_e_extracellular, text)

# stim = h.IClamp(simple_cell.dend(1))
# stim.delay = 0
# stim.dur = tstop
# stim.amp = 0.22

# Stimulate the cell
# stim = h.NetStim() 
# syn_ = h.ExpSyn(simple_cell.dend(0.5))
# syn_.tau = 2 * ms
# stim.number = 50
# stim.start = 0
# stim.noise = 0
# ncstim = h.NetCon(stim, syn_)
# ncstim.delay = 1 * ms
# ncstim.weight[0] = 0.04 

soma_v = h.Vector().record(simple_cell.soma(0.5)._ref_v)
dend_v = h.Vector().record(simple_cell.dend(0.5)._ref_v)
t = h.Vector().record(h._ref_t)

h.finitialize(-65 * mV)
h.continuerun(tstop * ms)


t_stim = np.zeros(len(text)*2)
t_stim[::2] = text
t_stim[1::2] = text 
t_stim = t_stim[1:]
v_temp = vextd
v_stim = np.zeros(len(v_temp)*2)
v_stim[::2] = v_temp
v_stim[1::2] = v_temp
v_stim = 0.5 * v_stim[:-1] / factor
v_stim = v_stim - np.abs(max(soma_v)) + 1


plt.figure()
plt.plot(t, soma_v, label="soma")
plt.plot(t, dend_v, "r--", label="dend")
plt.plot(t_stim, v_stim, "y", label="pulse")
plt.xlabel("t [ms]")
plt.ylabel("V [mV]")
plt.xlim([0, tstop])
plt.legend(loc="upper right")
plt.show()