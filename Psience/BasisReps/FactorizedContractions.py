"""Exact coefficient products for ordered harmonic words of arbitrary degree.

Only coefficient indices are factorized. No projector or resolvent is
inserted between factors of the quantum operator.
"""
import numpy as np
from numba import njit

__all__=['FactorizedCoefficients','FactorizedHarmonicContraction']


class FactorizedCoefficients:
    """C[i...,c...] = sum_r weights[r,c...] prod_f factors[f][r,i[axes[f]]].

    Each factor has a leading rank axis (length one broadcasts). Mode axes
    may overlap between factors, but must collectively cover every tensor
    mode axis. Component axes live in weights. Arrays are real; mode axes
    all have the same dimension. ``asarray`` explicitly materializes C.
    """
    _factorized_coefficients=True

    def __init__(self,factors,axes,weights=None):
        factors=tuple(np.asanyarray(f) for f in factors)
        axes=tuple(tuple(a) for a in axes)
        if not factors or len(factors)!=len(axes):
            raise ValueError('one axis tuple is required per coefficient factor')
        if any(np.iscomplexobj(f) for f in factors):
            raise TypeError('coefficient factors must be real')
        if any(f.ndim!=len(a)+1 or f.shape[0]<1 for f,a in zip(factors,axes)):
            raise ValueError('factors require a nonempty rank axis followed by their mode axes')
        flat=[i for a in axes for i in a]
        if any(not isinstance(i,(int,np.integer)) or i<0 for i in flat):
            raise ValueError('factor axes must be nonnegative integer positions')
        if any(len(set(a))!=len(a) for a in axes):
            raise ValueError('a factor cannot repeat a mode-axis label internally')
        degree=max(flat,default=-1)+1
        if set(flat)!=set(range(degree)):
            raise ValueError('factors must cover consecutive mode-axis labels')
        dimensions=[size for f,a in zip(factors,axes) for size in f.shape[1:]]
        if dimensions and (min(dimensions)<1 or len(set(dimensions))!=1):
            raise ValueError('all coefficient mode axes must have the same positive dimension')
        rank=max(f.shape[0] for f in factors)
        if weights is None:weights=np.ones(rank)
        weights=np.asanyarray(weights)
        if np.iscomplexobj(weights):raise TypeError('coefficient weights must be real')
        if weights.ndim<1 or weights.shape[0]<1:
            raise ValueError('weights require a nonempty leading rank axis')
        rank=max(rank,weights.shape[0])
        if weights.shape[0] not in [1,rank] or any(f.shape[0] not in [1,rank] for f in factors):
            raise ValueError('factor and weight ranks must agree or broadcast from one')
        self.factor_tensors=tuple(np.ascontiguousarray(f,dtype=float) for f in factors)
        self.factor_axes=axes
        self.weights=np.array(np.broadcast_to(weights,(rank,)+weights.shape[1:]),dtype=float,order='C',copy=True)
        self.rank,self.degree=rank,degree
        self.mode_shape=(dimensions[0],)*degree if dimensions else ()
        self.component_shape=self.weights.shape[1:]
        self.shape=self.mode_shape+self.component_shape
        self.ndim=len(self.shape)

    def asarray(self):
        result=np.empty(self.shape)
        for index in np.ndindex(self.mode_shape):
            value=np.ones(self.rank)
            for factor,axes in zip(self.factor_tensors,self.factor_axes):
                value*=factor[(slice(None),)+tuple(index[i] for i in axes)]
            result[index]=np.tensordot(value,self.weights,axes=1)
        return result


@njit(cache=True)
def _factor_contract(bras,kets,values,offsets,strides,ranks,axes,counts,weights,
                     words,transitions,tables,degree,phase,threshold):
    size,nmodes=bras.shape
    output=np.zeros((size,weights.shape[1]))
    minimum=np.empty(degree,dtype=np.int64)
    slots=np.empty(degree,dtype=np.int64)
    ordered=np.empty(degree,dtype=np.int64)
    choices=np.empty(degree//2,dtype=np.int64)
    modes=np.empty(degree,dtype=np.int64)
    subwords=np.empty(degree,dtype=np.int64)
    coefficient=np.empty(weights.shape[1])
    for pair in range(size):
        used=0;possible=True
        for mode in range(nmodes):
            delta=abs(int(bras[pair,mode])-int(kets[pair,mode]))
            if used+delta>degree:possible=False;break
            for j in range(delta):minimum[used]=mode;used+=1
        if not possible or (degree-used)%2:continue
        extra=(degree-used)//2
        for j in range(extra):choices[j]=0
        finished=False
        while not finished:
            for j in range(used):slots[j]=minimum[j]
            for j in range(extra):slots[used+2*j]=slots[used+2*j+1]=choices[j]
            slots.sort()
            for j in range(degree):ordered[j]=slots[j]
            while True:
                for c in range(weights.shape[1]):coefficient[c]=0.
                for r in range(weights.shape[0]):
                    product=1.
                    for f in range(len(offsets)):
                        index=0
                        for a in range(counts[f]):index=index*nmodes+ordered[axes[f,a]]
                        product*=values[offsets[f]+(0 if ranks[f]==1 else r)*strides[f]+index]
                    for c in range(weights.shape[1]):coefficient[c]+=weights[r,c]*product
                count=0
                for j in range(degree):
                    found=-1
                    for u in range(count):
                        if modes[u]==ordered[j]:found=u;break
                    if found<0:
                        found=count;modes[count]=ordered[j];subwords[count]=0;count+=1
                    subwords[found]=transitions[subwords[found],words[j]]
                value=phase
                for u in range(count):
                    mode=modes[u]
                    value*=tables[subwords[u],bras[pair,mode],kets[pair,mode]]
                if abs(value)>=threshold:
                    for c in range(weights.shape[1]):
                        output[pair,c]+=value*coefficient[c]
                # Next distinct multiset permutation; no factorial-size table.
                pivot=degree-2
                while pivot>=0 and ordered[pivot]>=ordered[pivot+1]:pivot-=1
                if pivot<0:break
                swap=degree-1
                while ordered[swap]<=ordered[pivot]:swap-=1
                ordered[pivot],ordered[swap]=ordered[swap],ordered[pivot]
                left,right=pivot+1,degree-1
                while left<right:
                    ordered[left],ordered[right]=ordered[right],ordered[left]
                    left+=1;right-=1
            if extra==0:finished=True
            else:
                j=extra-1
                while j>=0 and choices[j]==nmodes-1:j-=1
                if j<0:finished=True
                else:
                    next_value=choices[j]+1
                    for k in range(j,extra):choices[k]=next_value
    return output


class FactorizedHarmonicContraction:
    """Selection-driven, ordered contraction with no coefficient densification."""
    @classmethod
    def from_dense_operator(cls,operator):
        """Use the generic enumeration for a high-degree supplied dense tensor.

        This is an exact interface adapter, not a low-rank approximation.
        Components are selected by rank labels; the input mode tensor is
        still stored densely. The small-degree dense specialization remains
        available separately for faster canonical coefficient grouping.
        """
        import copy
        degree=len(operator.funcs.terms)
        c=np.asanyarray(operator.coeffs)
        if np.iscomplexobj(c):raise TypeError('direct contraction requires real coefficients')
        axes=operator.axes
        if axes is None:axes=(tuple(range(c.ndim)),)*2
        left,right=tuple(axes[0]),tuple(axes[1])
        if len(left)!=degree or sorted(left)!=list(range(degree)):
            raise ValueError('direct contractions must contract every operator axis')
        ordered=[right[left.index(i)]%c.ndim for i in range(degree)]
        remaining=[i for i in range(c.ndim) if i not in ordered]
        c=np.transpose(c,ordered+remaining)
        if c.shape[:degree]!=(operator.mode_n,)*degree:
            raise ValueError('coefficient mode axes do not match the oscillator basis')
        components=c.shape[degree:]
        component_count=int(np.prod(components)) if components else 1
        tensors=np.moveaxis(c.reshape((operator.mode_n,)*degree+(component_count,)),-1,0)
        weights=np.eye(component_count).reshape((component_count,)+components)
        adapted=copy.copy(operator)
        adapted.coeffs=FactorizedCoefficients([tensors],[tuple(range(degree))],weights)
        adapted.axes=None
        return cls(adapted)

    def __init__(self,operator):
        from .HarmonicOscillator import HarmonicProductOperatorTermEvaluator
        if not isinstance(operator.funcs,HarmonicProductOperatorTermEvaluator):
            raise TypeError('factor contractions require native harmonic evaluators')
        if operator.skipped_indices is not None or operator.coefficient_threshold is not None:
            raise ValueError('factorized coefficients require exact unfiltered support')
        self.word=tuple('x' if x=='q' else x for x in operator.funcs.terms)
        c=operator.coeffs
        self.degree=len(self.word)
        self.nmodes=operator.mode_n
        if self.degree!=c.degree or c.mode_shape!=(self.nmodes,)*self.degree:
            raise ValueError('factor mode axes do not match the harmonic operator')
        self.component_shape=c.component_shape
        ordered_axes=tuple(range(self.degree))
        if operator.axes is not None:
            left,right=tuple(operator.axes[0]),tuple(operator.axes[1])
            if sorted(left)!=list(range(self.degree)):
                raise ValueError('factor contractions must contract every operator axis')
            ordered_axes=tuple(right[left.index(i)]%c.ndim for i in range(self.degree))
            if sorted(ordered_axes)!=list(range(self.degree)):
                raise ValueError('factor contractions must contract exactly the coefficient mode axes')
        self.weights=c.weights.reshape(c.rank,-1)
        factors=c.factor_tensors
        self.values=np.concatenate([f.ravel() for f in factors])
        sizes=np.array([f.size for f in factors],dtype=np.int64)
        self.offsets=np.concatenate([np.zeros(1,dtype=np.int64),np.cumsum(sizes)[:-1]])
        self.ranks=np.array([f.shape[0] for f in factors],dtype=np.int64)
        self.strides=sizes//self.ranks
        self.counts=np.array([len(a) for a in c.factor_axes],dtype=np.int64)
        self.axes=np.full((len(factors),max(self.counts,default=0)),-1,dtype=np.int64)
        for i,axes in enumerate(c.factor_axes):
            self.axes[i,:len(axes)]=[ordered_axes.index(a) for a in axes]
        self.threshold=operator.zero_threshold
        self.phase=operator.funcs.i_phase
        # Intern distinct native subwords rather than allocating 2**degree
        # bit masks, many of which represent the same ordered one-mode word.
        self.subwords=[()];lookup={():0}
        for op in self.word:
            for previous in tuple(self.subwords):
                word=previous+(op,)
                if word not in lookup:
                    lookup[word]=len(self.subwords);self.subwords.append(word)
        self.words=np.array([0 if op=='x' else 1 for op in self.word],dtype=np.int64)
        self.transitions=np.full((len(self.subwords),2),-1,dtype=np.int64)
        for i,word in enumerate(self.subwords):
            for j,op in enumerate(['x','p']):self.transitions[i,j]=lookup.get(word+(op,),-1)
        self.extent=0;self.tables=None

    def evaluate(self,pairs):
        from .HarmonicOscillator import HarmonicOscillatorMatrixGenerator
        bras,kets=pairs.state_pairs
        bras,kets=np.ascontiguousarray(bras.T),np.ascontiguousarray(kets.T)
        if self.degree==0:
            scalar=np.ones(self.weights.shape[0])
            for i,rank in enumerate(self.ranks):
                scalar*=self.values[self.offsets[i]+np.arange(len(scalar))%rank]
            value=scalar@self.weights
            return (np.all(bras==kets,axis=1)[:,None]*value).reshape((len(bras),)+self.component_shape)
        if not len(bras):return np.empty((0,)+self.component_shape)
        extent=int(max(np.max(bras),np.max(kets)))+1
        if extent>self.extent:
            self.tables=np.zeros((len(self.subwords),extent,extent))
            i,j=np.indices((extent,extent))
            for index,word in enumerate(self.subwords[1:],1):
                generator=HarmonicOscillatorMatrixGenerator.load_cached(word)
                self.tables[index]=generator.evaluate_state_terms((i.ravel(),j.ravel())).reshape(extent,extent)
            self.extent=extent
        result=_factor_contract(bras,kets,self.values,self.offsets,self.strides,self.ranks,
            self.axes,self.counts,self.weights,self.words,self.transitions,self.tables,
            self.degree,self.phase,self.threshold)
        return result.reshape((len(bras),)+self.component_shape)
