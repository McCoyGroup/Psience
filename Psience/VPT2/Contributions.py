"""Optional bra/ket/operator order diagnostics using native analytic evaluation."""
import collections
from collections.abc import Mapping
import numbers
import numpy as np
from McUtils.Formatters import TableFormatter

__all__=['PerturbationTermOrder','AnalyticTermContributions']

class PerturbationTermOrder(collections.namedtuple('_PerturbationTermOrder',
        ['bra_order','ket_order','operator_order'])):
    """Matrix-code key convention: (bra order, ket order, operator order)."""
    __slots__=()
    def __new__(cls,bra_order,ket_order,operator_order):
        values=(bra_order,ket_order,operator_order)
        if any(isinstance(v,(bool,np.bool_)) or not isinstance(v,numbers.Integral) or v<0 for v in values):
            raise ValueError('bra, ket and operator orders must be nonnegative integers')
        return super().__new__(cls,*(int(v) for v in values))
    @property
    def order(self):return sum(self)
    @property
    def analytic_filter(self):return self.ket_order,self.operator_order,self.bra_order
    def label(self,symbol='A'):
        return '<{}|{}({})|{}>'.format(self.bra_order,symbol,self.operator_order,self.ket_order)

def _copy_pack(pack):return [[np.array(matrix,copy=True) for matrix in component] for component in pack]
def _sum_packs(packs):
    iterator=iter(packs);result=_copy_pack(next(iterator))
    for pack in iterator:
        if len(pack)!=len(result):raise ValueError('inconsistent operator component count')
        for i,(left,right) in enumerate(zip(result,pack)):
            if len(left)!=len(right):raise ValueError('inconsistent state block count')
            for j,(a,b) in enumerate(zip(left,right)):
                if a.shape!=np.shape(b):raise ValueError('inconsistent state matrix shapes')
                result[i][j]=a+b
    return result

class AnalyticTermContributions(Mapping):
    """Physical matrix elements by (bra,ket,operator) order, before final resonance rotation.

    Each value is component -> state block -> array[bra state,ket state].
    Individual terms are not symmetrized; Hermitian partners remain distinct.
    """
    def __init__(self,contributions,*,symbol,units,state_blocks,component_names=None,
                 metadata=None,reference=None,atol=1e-11,rtol=1e-9):
        self.contributions=dict(contributions);self.symbol=symbol;self.units=units
        if not self.contributions:raise ValueError('at least one contribution is required')
        self.state_blocks=[dict(bra_states=np.asarray(b['bra_states'],dtype=int),
                                ket_states=np.asarray(b['ket_states'],dtype=int)) for b in state_blocks]
        first=next(iter(self.contributions.values()));count=len(first)
        self.component_names=tuple(component_names or [symbol if count==1 else '{}_{}'.format(symbol,i) for i in range(count)])
        if len(self.component_names)!=count:raise ValueError('component_names does not match operator components')
        self.metadata=dict(metadata or {},term_key_order=['bra_order','ket_order','operator_order'],
            native_filter_order=['ket_order','operator_order','bra_order'],matrix_axes=['bra_state','ket_state'],
            native_operator_correction_to_matrix_element=-1,final_resonance_rotation_applied=False)
        self.reference=reference;self.reconstruction_error=None
        for key,pack in self.contributions.items():
            if len(pack)!=count:raise ValueError('inconsistent component count')
            for component in pack:
                if len(component)!=len(self.state_blocks):raise ValueError('inconsistent block count')
                for matrix,states in zip(component,self.state_blocks):
                    if np.shape(matrix)!=(len(states['bra_states']),len(states['ket_states'])):
                        raise ValueError('contribution matrix does not match its explicit state axes')
        if reference is not None:self.validate(atol=atol,rtol=rtol)
    def __getitem__(self,key):return self.contributions[key]
    def __iter__(self):return iter(self.contributions)
    def __len__(self):return len(self.contributions)
    @property
    def orders(self):return tuple(sorted({key.order for key in self}))
    @property
    def by_order(self):return {n:{key:value for key,value in self.items() if key.order==n} for n in self.orders}
    def term(self,key,component=0,block=0):return self[key][component][block]
    def sum(self,order,component=None,block=None):
        if order not in self.orders:raise KeyError(order)
        pack=_sum_packs(value for key,value in self.items() if key.order==order)
        if component is not None:pack=pack[component]
        if block is not None:
            if component is None:raise ValueError('specify component when selecting a block')
            pack=pack[block]
        return pack
    def validate(self,*,atol=1e-11,rtol=1e-9):
        if self.reference is None:raise ValueError('no independent aggregate reference was requested')
        largest=0.
        for n in self.orders:
            actual=self.sum(n);expected=self.reference[n]
            if len(actual)!=len(expected):raise ValueError('reference component count does not match contributions')
            for a_component,e_component in zip(actual,expected):
                if len(a_component)!=len(e_component):raise ValueError('reference block count does not match contributions')
                for a,e in zip(a_component,e_component):
                    if np.shape(a)!=np.shape(e):raise ValueError('reference matrix shape does not match contributions')
                    largest=max(largest,float(np.max(np.abs(a-e),initial=0.)))
                    np.testing.assert_allclose(a,e,atol=atol,rtol=rtol,
                        err_msg='bra/ket/operator terms do not reconstruct native aggregate at order {}'.format(n))
        self.reconstruction_error=largest
        return largest
    def rows(self,order=None,include_zeros=True):
        records=[]
        for key,pack in self.items():
            if order is not None and key.order!=order:continue
            for component,blocks in enumerate(pack):
                for block,(matrix,states) in enumerate(zip(blocks,self.state_blocks)):
                    for bra,ket in np.ndindex(matrix.shape):
                        value=matrix[bra,ket].item()
                        if not include_zeros and value==0:continue
                        records.append(dict(order=key.order,term=tuple(key),label=key.label(self.symbol),
                            component=self.component_names[component],block=block,
                            bra_state=tuple(int(n) for n in states['bra_states'][bra]),
                            ket_state=tuple(int(n) for n in states['ket_states'][ket]),value=value))
        return records
    def format_table(self,order=None,*,conversion=1.,units=None,include_zeros=True,number_format='{:.12g}'):
        """Native TableFormatter output; conversion affects display only."""
        rows=self.rows(order,include_zeros)
        body=[[r['order'],r['label'],r['component'],r['block'],str(r['bra_state']),str(r['ket_state']),r['value']*conversion] for r in rows]
        formatter=TableFormatter(['{}','{}','{}','{}','{}','{}',number_format],
            headers=[['Order','Term','Component','Block','Bra state','Ket state','Value [{}]'.format(self.units if units is None else units)]],column_join=' | ')
        return formatter.format(body)
    def to_dict(self):
        def array(value):
            value=np.asarray(value)
            return dict(real=value.real.tolist(),imag=value.imag.tolist()) if np.iscomplexobj(value) else value.tolist()
        def pack(value):return [[array(matrix) for matrix in component] for component in value]
        return dict(symbol=self.symbol,units=self.units,component_names=list(self.component_names),metadata=self.metadata,
            state_blocks=[{name:value.tolist() for name,value in states.items()} for states in self.state_blocks],
            contributions=[dict(term=list(key),order=key.order,label=key.label(self.symbol),values=pack(value)) for key,value in self.items()],
            reference=None if self.reference is None else {str(n):pack(value) for n,value in self.reference.items()},
            reconstruction_max_abs_error=self.reconstruction_error)
    @classmethod
    def collect(cls,evaluate,order,*,orders=None,term_orders=None,validate=True,**options):
        if isinstance(order,(bool,np.bool_)) or not isinstance(order,numbers.Integral) or order<0:raise ValueError('order must be a nonnegative integer')
        order=int(order)
        selected_orders=list(range(order+1)) if orders is None else list(orders)
        if not selected_orders or any(not isinstance(n,numbers.Integral) or isinstance(n,(bool,np.bool_)) or n<0 or n>order for n in selected_orders):
            raise ValueError('orders must be nonempty integers between zero and order')
        if term_orders is None:
            keys=[PerturbationTermOrder(bra,ket,n-bra-ket) for n in sorted(set(selected_orders)) for bra in range(n+1) for ket in range(n-bra+1)]
        else:
            keys=list(dict.fromkeys(PerturbationTermOrder(*key) for key in term_orders))
            if not keys or any(key.order>order or key.order not in selected_orders for key in keys):raise ValueError('selected term orders lie outside the requested orders')
        contributions={key:evaluate([key.analytic_filter],key.order)[key.order] for key in keys}
        reference=None
        if validate:
            reference=evaluate(None if term_orders is None else [key.analytic_filter for key in keys],order)
            reference={n:reference[n] for n in sorted({key.order for key in keys})}
        return cls(contributions,reference=reference,**options)
