"""Effective vertices preserve raw V/G inputs, contacts and saved execution."""
from itertools import permutations, product
from fractions import Fraction
from pathlib import Path
import json
import subprocess
import sys
from unittest import mock

import numpy as np
import pytest
from Psience.VPT2 import Analytic as a, Evaluators as e, AnalyticVPTRunner
from PortableAnalyticCompilationTests import simple_run_inputs


def inputs(n=3, order=4, seed=471):
    rng=np.random.default_rng(seed)
    frequencies=.017*np.sqrt(np.arange(2,n+2))
    coefficients=[[np.diag(frequencies)/2,np.diag(frequencies)/2]]
    for degree in range(1,order+1):
        rank=degree+2
        f=rng.normal(size=(n,)*rank)*1.e-5
        # Assign one value per sorted multi-index without factorial work.
        f=np.array([f[tuple(sorted(i))] for i in product(range(n),repeat=rank)]).reshape((n,)*rank)
        g=rng.normal(size=(n,)*rank)*1.e-5
        g=np.array([g[(min(i[0],i[-1]),)+tuple(sorted(i[1:-1]))+(max(i[0],i[-1]),)]
                    for i in product(range(n),repeat=rank)]).reshape((n,)*rank)
        coefficients.append([f,g,0,0,0])
    states=[(np.zeros(n,dtype=int),np.arange(n)[None]),
            (np.arange(n)+2,np.array(list(permutations(range(n)))))]
    return states,[coefficients],frequencies


def primitive(order,change,combined):
    solver=a.AnalyticPerturbationTheorySolver.from_order(max(2,order),logger=False,
        combine_kinetic_coefficients=combined)
    term=solver.hamiltonian_expansion[order]
    return a.PerturbationTheoryExpressionEvaluator(term,term.get_poly_terms(change),change,logger=False)


def resonant_inputs():
    states,coefficients,frequencies=inputs(n=2,order=3)
    states=states[:1]
    for row in coefficients[0]:row.append(0)
    for order in range(1,4):coefficients[0][order][5]=np.ones((2,)*order)*.03**order
    groups=(((1,0),(0,2)),)
    return (states,coefficients,frequencies),dict(projectors=(groups,))


def transition(branch,combined):
    solver=a.AnalyticPerturbationTheorySolver.from_order(2,logger=False,combine_kinetic_coefficients=combined)
    solver.filtered_block_recurrence=True
    source=solver._action_expansion_source('hamiltonian',2,resonance_correction=True)
    term=(source.operator_correction(2,operator_type='transition_moment') if branch=='main' else
        getattr(source.operator_degenerate_correction(2,operator_type='transition_moment',degenerate_changes=()),branch))
    return a.PerturbationTheoryExpressionEvaluator(term,term.get_poly_terms((1,)),(1,),logger=False)


@pytest.mark.parametrize('branch',['main','left','right'])
@pytest.mark.parametrize('mode',['direct','actions'])
def test_full_filtered_transition_and_normalization_branches(branch,mode):
    arguments,bindings=resonant_inputs()
    raw=transition(branch,False).evaluate(*arguments,projection_bindings=bindings,evaluation_mode=mode,zero_cutoff=0)
    folded=transition(branch,True).evaluate(*arguments,projection_bindings=bindings,evaluation_mode=mode,zero_cutoff=0)
    for left,right in zip(raw,folded):np.testing.assert_allclose(right,left,rtol=3.e-11,atol=2.e-13)


def test_fourth_order_recurrence_with_all_internal_vertices():
    arguments=inputs(n=1,order=4)
    values=[]
    for combined in (False,True):
        solver=a.AnalyticPerturbationTheorySolver.from_order(4,logger=False,combine_kinetic_coefficients=combined)
        wrapper=solver.energy_correction(4)(())
        compiled=e.AnalyticExpressionStore().compile(wrapper,mode_dimension=1)
        values.append(compiled.evaluate(*arguments,zero_cutoff=0))
    for left,right in zip(*values):np.testing.assert_allclose(right,left,rtol=3.e-11,atol=2.e-15)


@pytest.mark.parametrize('rank',range(3,11))
def test_exact_polynomial_size_channel_basis(rank):
    signs,inverse=a._HamiltonianCoefficient.basis(rank)
    pairs=list(__import__('itertools').combinations(range(rank),2))
    assert len(signs)==1+len(pairs)
    matrix=[(Fraction(1),)+tuple(Fraction(-s[i]*s[j],len(pairs)) for i,j in pairs) for s in signs]
    for i,j in product(range(len(matrix)),repeat=2):
        assert sum(matrix[i][k]*inverse[k][j] for k in range(len(matrix))) == (i==j)


@pytest.mark.parametrize('order,changes',[(1,[(1,),(-1,),(3,),(-3,),(1,1,1),(1,-2)]),
    (2,[(),(2,),(-2,),(1,1),(1,-1),(4,)]),(3,[(1,),(-1,),(3,),(5,)]),(4,[(),(2,),(-2,),(4,),(6,)])])
@pytest.mark.parametrize('source', ['both','potential','kinetic'])
def test_primitive_through_sextic_repeated_indices_and_hot_bands(order,changes,source):
    arguments=inputs(order=max(2,order))
    if source!='both':
        for row in arguments[1][0][1:]:row[0 if source=='kinetic' else 1]=0
    for change in changes:
        if len(change)>3:continue
        raw=primitive(order,change,False).evaluate(*arguments,zero_cutoff=0,evaluation_backend='indexed')
        folded=primitive(order,change,True).evaluate(*arguments,zero_cutoff=0,evaluation_backend='compare')
        for left,right in zip(raw,folded):np.testing.assert_allclose(right,left,rtol=2.e-12,atol=2.e-14)


@pytest.mark.parametrize('order',range(1,5))
def test_contacts_through_sextic_against_native_oscillator_matrices(order):
    from Psience.BasisReps import HarmonicOscillatorBasis
    rank=order+2
    _,coefficients,frequencies=inputs(n=1,order=4)
    q=HarmonicOscillatorBasis.qmatrix_ho(16).toarray()
    d=HarmonicOscillatorBasis.pmatrix_ho(16).toarray()
    v,g=[tensor.reshape(-1)[0] for tensor in coefficients[0][order][:2]]
    matrix=v*np.linalg.matrix_power(q,rank)-g*(d@np.linalg.matrix_power(q,rank-2)@d)
    for delta in range(rank%2,rank+1,2):
        change=() if delta==0 else (delta,)
        states=[(np.array([n]),np.array([[0]])) for n in range(4)]
        actual=primitive(order,change,True).evaluate(states,coefficients,frequencies,zero_cutoff=0)
        np.testing.assert_allclose(np.asarray(actual).reshape(-1),
            [matrix[n+delta,n] for n in range(4)],rtol=2.e-12,atol=2.e-14)


@pytest.mark.parametrize('ir_backend,plan_backend',list(product(('portable','pickle'),repeat=2)))
def test_standalone_transfer_never_rebuilds_recipe(tmp_path,ir_backend,plan_backend):
    store=e.AnalyticExpressionStore()
    arguments=inputs(order=2)
    expected=primitive(1,(1,),False).evaluate(*arguments)
    ir=tmp_path/'vertex.ir.tar';plan=tmp_path/'vertex.plan.tar'
    with store.export_ir(primitive(1,(1,),True),ir,storage_backend=ir_backend):pass
    # Neither independent compilation nor recipient evaluation needs a solver,
    # an oscillator polynomial generator or the channel-basis constructor.
    with mock.patch.object(a._HamiltonianCoefficient,'basis',side_effect=AssertionError('basis reconstruction')), \
         mock.patch.object(a.AnalyticPerturbationTheorySolver,'__init__',side_effect=AssertionError('solver')):
        saved=store.compile_ir(ir,plan,storage_backend=plan_backend)
        actual=saved.evaluate(*arguments)
    for left,right in zip(expected,actual):np.testing.assert_allclose(right,left,rtol=2.e-12,atol=2.e-14)
    ir.unlink()
    recipient_expected=primitive(1,(1,),False).evaluate(*inputs(n=4,order=2))
    code='''
import sys,json
from unittest import mock
sys.path.insert(0,sys.argv[1])
from EffectiveHamiltonianTests import inputs
from Psience.VPT2 import Analytic as a,Evaluators as e
def fail(*args,**kwargs):raise AssertionError('structural reconstruction')
with mock.patch.object(a._HamiltonianCoefficient,'basis',fail), \\
     mock.patch.object(a.AnalyticPerturbationTheorySolver,'__init__',fail), \\
     mock.patch.object(e.IndexedBlockEvaluator,'compile_plan',fail), \\
     mock.patch.object(e.IndexedCoefficientGatherPlan,'_get_request_pattern_map',fail):
    saved=e.CompiledAnalyticExpression.load(sys.argv[2])
    print(json.dumps([v.tolist() for v in saved.evaluate(*inputs(n=4,order=2))]))
'''
    child=subprocess.run([sys.executable,'-c',code,str(Path(__file__).parent),str(plan)],capture_output=True,text=True)
    assert child.returncode==0,child.stderr
    for left,right in zip(recipient_expected,json.loads(child.stdout)):np.testing.assert_allclose(right,left,rtol=2.e-12,atol=2.e-14)


@pytest.mark.parametrize('backend',['portable','pickle'])
@pytest.mark.parametrize('mode',['direct','actions'])
def test_run_simple_internal_compiled_reuse_and_raw_control(tmp_path,backend,mode):
    system,states,options=simple_run_inputs(2)
    options['internals']=[[0,-1,-1,-1],[1,0,-1,-1],[2,0,1,-1]]
    states.append(np.array([0,2]))
    options.update(evaluation_mode=mode,
        degeneracy_specs={'groups':[[states[1],states[2]]],'group_filter':'unfiltered'})
    options['kinetic_terms'][1]=np.ones((2,)*3)*1.e-5
    options['kinetic_terms'][2]=np.ones((2,)*4)*1.e-6
    raw=AnalyticVPTRunner.run_simple(system,states,combine_kinetic_coefficients=False,**options)
    target=tmp_path/'compiled.tar'
    actual=AnalyticVPTRunner.run_simple(system,states,compiled_file=target,storage_backend=backend,**options)
    has_recipe=False
    for name in e.CompiledAnalyticExpression.expression_names(target):
        with e._CompiledAnalyticArchive.load(target,expression=name) as reader:
            for part in reader.iter_parts():
                has_recipe |= any(isinstance(c[1],a._HamiltonianCoefficient) for key in part.terms for c in key)
    assert has_recipe
    for field in ('energy_corrections','transition_moment_corrections','degenerate_hamiltonian_corrections'):
        np.testing.assert_allclose(getattr(actual,field),getattr(raw,field),rtol=2.e-10,atol=2.e-13)
    with mock.patch.object(a.PerturbationTheoryTerm,'get_poly_terms',side_effect=AssertionError('rebuild')), \
         mock.patch.object(a._HamiltonianCoefficient,'basis',side_effect=AssertionError('basis')):
        again=AnalyticVPTRunner.run_simple(system,states,compiled_file=target,compilation_options={'reuse_only':True},**options)
    np.testing.assert_allclose(again.transition_moment_corrections,actual.transition_moment_corrections,atol=2.e-13)
