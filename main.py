"""import scipy"""
from scipy.integrate import odeint
from numpy import *
import matplotlib.pyplot as plt
from matplotlib import *
"""from scipy import weave"""


class IntStruct:
	def get_size(self):
		pass

	def __init__(self):
		pass

	def ddt(self, y, t):
		pass

	def augment_y0(self, y0):
		pass

class HH_RS(IntStruct):

	num_internal_vars = 4

	@staticmethod
	def m0(V):
		return 1./(.1+exp((-V - 34.5)/10.))
	@staticmethod
	def hinf(V):
		return 1./(1.+exp((V + 59.4)/10.7))
	@staticmethod
	def tauh(V):
		return .15 + 1.15/(1+exp((V+33.5)/15.))
	@staticmethod
	def minf(V):
		return 1./(1.+exp((-V-29.5)/10.))
	@staticmethod
	def taum(V):
		return .25 + 4.35*exp(-abs(V+10.)/10.)
	@staticmethod
	def m_hinf(V, V_0):
		return 1./(1.+exp((V+V_0)/5.5))
	@staticmethod
	def taum_h(V):
		return 1./(exp(-14.6 - 0.086*V) + exp(-1.87+0.07*V))

	
	def __init__(self,
		J,
		initial_index,
		V_pre_s_indices = [],
		V_pre_gj_indices = [],
		gNaF = 200.,
		gKDR = 200.,
		g_h=1.,
		V_rev_s_ = [],
		g_s_ = [],
		g_gj_ = [],
		tau_r_ = [],
		tau_d_ = [],
		V_0 = 75.,
		C = 1.,
		V0 = -75.,
		h0 = .05,
		m0 = .5,
		m_h0 = 0.,
		s0_ = []):

		"""or g_h =  25."""

		self.num_s = len(V_pre_s_indices)
		self.num_gj = len(V_pre_gj_indices)
		self.num_indices = HH_RS.num_internal_vars + self.num_s + self.num_gj

		"""or V0 = 87.5"""
		self.V_index = initial_index
		self.h_index = self.V_index+1
		self.m_index = self.h_index+1
		self.m_h_index = self.m_index+1
		self.s_indices = range(self.m_h_index+1, self.m_h_index + self.num_s+1)

		self.J = J
		self.V_pre_s_indices = V_pre_s_indices
		self.V_pre_gj_indices = V_pre_gj_indices
		self.gNaF = gNaF
		self.gKDR = gKDR
		self.g_h = g_h
		self.V_rev_s_ = array(V_rev_s_)
		self.g_s_ = array(g_s_)
		self.g_gj_ = array(g_gj_)
		self.tau_r_ = array(tau_r_)
		self.tau_d_ = array(tau_d_)
		self.V_0 =V_0
		self.C = C
		
		self.V0 = V0
		self.h0 = h0
		self.m0 = m0
		self.m_h0 = m_h0
		self.s0_ = array(s0_)

	def get_size(self):
		return self.num_indices

	def ddt(self, y, t):
		y = array(y)
		V = y[self.V_index]
		h = y[self.h_index]
		m = y[self.m_index]
		m_h = y[self.m_h_index]
		s_ = array(y[self.s_indices])

		V_pre_s_ = array(y[self.V_pre_s_indices])
		V_pre_gj_ = array(y[self.V_pre_gj_indices])

		dydt = zeros(size(y))

		dydt[self.V_index] = (1./self.C) * (-self.J - (70.+V)
			- self.gNaF * HH_RS.m0(V)**3 * h * (-50.+V)
			- self.gKDR * m**4 * (95. + V)
			- self.g_h * m_h * (35. + V)
			- sum(self.g_s_ * (V - self.V_rev_s_) * s_)
			- sum(self.g_gj_ * (V - V_pre_gj_)))
		dydt[self.h_index] = (HH_RS.hinf(V) - h)/HH_RS.tauh(V)
		dydt[self.m_index] = (HH_RS.minf(V) - m)/HH_RS.taum(V)
		dydt[self.m_h_index] = (HH_RS.m_hinf(V, self.V_0) - m_h)/HH_RS.taum_h(V)
		dydt[self.s_indices] = -s_/self.tau_d_ + (1.-s_) * (1.+tanh(V_pre_s_/10.)) / self.tau_r_

		return dydt
	def augment_y0(self, y0):
		I = [self.V_index, self.h_index, self.m_index, self.m_h_index] + list(self.s_indices)
		L = [self.V0, self.h0, self.m0, self.m_h0]+ list(self.s0_)
		y0[I] = L


class HH_FS(IntStruct):


	num_internal_vars = 3

	@staticmethod
	def m0(V):
		return 1./(1.+exp((-V - 38.)/10.))
	@staticmethod
	def hinf(V):
		return 1./(1.+exp((V + 58.3)/6.7))
	@staticmethod
	def tauh(V):
		return .225 + 1.125/(1.+exp((V+37.)/15.))
	@staticmethod
	def minf(V):
		return 1./(1.+exp((-V-27.)/11.5))
	@staticmethod
	def taum(V):
		return .25 + 4.35*exp(-abs(V+10)/10.)

	
	def __init__(self,
		initial_index,
		J,
		V_pre_s_indices = [],
		V_rev_s_ = [],
		V_pre_gj_indices = [],
		gNaF = 200.,
		gKDR = 200.,
		g_s_ = [],
		g_gj_ = [],
		tau_r_ = [],
		tau_d_ = [],
		C = 1.,
		V0 = -75.,
		h0 = .05,
		m0 = .5,
		s0_ = []):

		"""or V0 = 87.5"""

		self.num_s = len(V_pre_s_indices)
		self.num_gj = len(V_pre_gj_indices)
		self.num_indices = HH_FS.num_internal_vars + self.num_s + self.num_gj


		self.V_index = initial_index
		self.h_index = self.V_index + 1
		self.m_index = self.h_index + 1
		self.s_indices = range(self.m_index+1, self.m_index + self.num_s+1)

		self.J = J
		self.V_pre_s_indices = V_pre_s_indices
		self.V_pre_gj_indices = V_pre_gj_indices
		self.gNaF = gNaF
		self.gKDR = gKDR
		self.V_rev_s_ = array(V_rev_s_)
		self.g_s_ = array(g_s_)
		self.g_gj_ = array(g_gj_)
		self.tau_r_ = array(tau_r_)
		self.tau_d_ = array(tau_d_)
		self.C = C

		self.V0 = V0
		self.h0 = h0
		self.m0 = m0
		self.s0_ = array(s0_)

	def get_size(self):
		return self.num_indices


	def ddt(self, y, t):
		y = array(y)
		V = y[self.V_index]
		h = y[self.h_index]
		m = y[self.m_index]
		s_ = array(y[self.s_indices])

		V_pre_s_ = array(y[self.V_pre_s_indices])
		V_pre_gj_ = array(y[self.V_pre_gj_indices])

		dydt = zeros(size(y))

		dydt[self.V_index] = (1./self.C) * (-self.J - (65.+V)
			- self.gNaF * HH_RS.m0(V)**3 * h * (-50.+V)
			- self.gKDR * m**4 * (100. + V)
			- sum(self.g_s_ * (V - self.V_rev_s_) * s_)
			- sum(self.g_gj_ * (V - V_pre_gj_)))
		dydt[self.h_index] = (HH_FS.hinf(V) - h)/HH_FS.tauh(V)
		dydt[self.m_index] = (HH_FS.minf(V) - m)/HH_FS.taum(V)
		dydt[self.s_indices] = -s_/self.tau_d_ + (1.-s_) * (1.+tanh(V_pre_s_/10.)) / self.tau_r_

		return dydt

	def augment_y0(self,y0):
		y0[[self.V_index, self.h_index, self.m_index] + list(self.s_indices)] = [self.V0, self.h0, self.m0]+ list(self.s0_)

class EI(IntStruct):

	def __init__(self,
		EForcing = (lambda t: 0),
		IForcing = (lambda t: 0),
		Je = -8.,
		Ji = 3.8,
		g_ie = .5,
		g_ei = .15,
		g_ii = 20.,
		g_ee = 0.,
		g_iChR2 = 1.,
		g_eChR2 = 1.,
		tau_ie_r = .5,
		tau_eChR2_r = .1,
		tau_ei_r = .25,
		tau_ii_r = .5,
		tau_iChR2_r = .1,
		tau_ee_r = .25,
		tau_ie_d = 6.,
		tau_eChR2_d = .1,
		tau_ei_d = 1.,
		tau_ii_d = 5.,
		tau_iChR2_d = .1,
		tau_ee_d = 1.):

		self.EForcing = EForcing
		self.IForcing = IForcing

		self.RS = HH_RS(
			J = Je,
			initial_index = 0,
			V_pre_s_indices = [7, 0, 13],
			g_s_ = [g_ie, g_ee, g_eChR2],
			tau_r_ = [tau_ie_r, tau_ee_r, tau_eChR2_r],
			tau_d_ = [tau_ie_d, tau_ee_d, tau_eChR2_d],
			V_rev_s_ = [-80., 0., 0.], 
			s0_ = [0,0,0])
		self.FS = HH_FS(
			J = Ji,
			initial_index = 7,
			V_pre_s_indices = [0, 7, 14],
			g_s_ = [g_ei, g_ii, g_iChR2],
			tau_r_ = [tau_ei_r, tau_ii_r, tau_iChR2_r],
			tau_d_ = [tau_ei_d, tau_ii_d, tau_iChR2_d],
			V_rev_s_ = [0., -75., 0.],
			s0_ = [0,0,0])
		"""[-g_ie, g_ee, g_eChR2]"""
		"""[g_ei, -g_ii, g_iChR2]"""
	def get_size(self):
		return self.RS.get_size() + self.FS.get_size()

	def ddt(self, y, t):
		x = list(y)+[self.EForcing(t), self.IForcing(t)];
		return array(self.RS.ddt(x, t)) + array(self.FS.ddt(x, t))

	def augment_y0(self, y0):
		self.RS.augment_y0(y0)
		self.FS.augment_y0(y0)
		

Exp1 = EI()
y0 = zeros(Exp1.get_size())
Exp1.augment_y0(y0)


tlength = 50
dt = .01
	
t = arange(0., tlength/dt)*dt
y = odeint(Exp1.ddt, y0, t, rtol = .0001, atol = .0001)

plt.subplot(611)
plt.plot(t, y[:, Exp1.RS.V_index], 'r-', t, y[:, Exp1.FS.V_index], 'b-')
plt.ylabel('V')
plt.subplot(612)
plt.plot(t, y[:, Exp1.RS.m_index], 'r-', t, y[:, Exp1.FS.m_index], 'b-')
plt.ylabel('m')
plt.subplot(613)
plt.plot(t, y[:, Exp1.RS.h_index], 'r-', t, y[:, Exp1.FS.h_index], 'b-')
plt.ylabel('h')
plt.subplot(614)
plt.plot(t, y[:, Exp1.RS.s_indices[1]], 'r-', t, y[:, Exp1.FS.s_indices[0]], 'b-')
plt.ylabel('E')
plt.subplot(615)
plt.plot(t, y[:, Exp1.RS.s_indices[0]], 'r-', t, y[:, Exp1.FS.s_indices[1]], 'b-')
plt.ylabel('I')
plt.subplot(616)
plt.plot(t, y[:, Exp1.RS.m_h_index], 'r-')
plt.ylabel('m_h')



plt.show()

