from Psience.AnalyticModels import SymbolicCaller, AnalyticModelBase
sym = AnalyticModelBase.sym
m = SymbolicCaller(AnalyticModelBase.symbolic_m)
r = SymbolicCaller(AnalyticModelBase.symbolic_r)
a = SymbolicCaller(AnalyticModelBase.symbolic_a)
t = SymbolicCaller(AnalyticModelBase.symbolic_t)
y = SymbolicCaller(AnalyticModelBase.symbolic_y)
sin = sym.sin
cos = sym.cos
cot = sym.cot
tan = sym.tan
csc = sym.csc
L = SymbolicCaller(AnalyticModelBase.lam)
data = {k: (g[k], vp[k]) for k in g}