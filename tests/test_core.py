import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np
import pytest
from core import (PROFILES,Geometry,Conditions,material_library,outer_radius,compile_load,
                  _slab_power,temperature_model,optimize_one)
from exporters import stl_bytes,step_bytes


def test_rf_air_material_energy_conservation():
    t=np.array([.001,.003,.007])[None,None,:]
    for pol in ("TE","TM"):
        p=_slab_power(1.0,0.0,t,10e9,np.deg2rad(0)*0+30,pol)
        assert np.allclose(p,1.,atol=1e-11)


def test_rf_loss_is_passive():
    for pol in ("TE","TM"):
        p=_slab_power(3.2,.01,np.linspace(.001,.014,50),10e9,45,pol)
        assert np.all(p<=1.) and np.all(p>0.)


@pytest.mark.parametrize('name',PROFILES)
def test_profile_endpoints_and_monotonicity(name):
    geo=Geometry(profile=name)
    x=np.linspace(0,geo.length_m,301)
    r=outer_radius(geo,x)
    assert np.isclose(r[0],0,atol=1e-8)
    assert np.isclose(r[-1],geo.base_diameter_m/2,atol=1e-8)
    assert np.min(np.diff(r))>-1e-8


def test_blunted_tip_endpoints():
    g=Geometry(spherical_tip_radius_m=.015)
    x=np.linspace(0,g.length_m,301)
    y=outer_radius(g,x)
    assert abs(y[0])<1e-10
    assert abs(y[-1]-g.base_diameter_m/2)<1e-8


def test_no_unsafe_expressions():
    expr=compile_load('1000*u + 20*sin(pi*u)')
    assert np.allclose(expr(np.array([0.,.4]),.4),[0.,1000.])
    for bad in ('__import__("os").system("x")','x.__class__','[1,2,3]','open("file")'):
        with pytest.raises(ValueError): compile_load(bad)


def test_temperature_inner_between_skin_and_air():
    c=Conditions()
    g=Geometry()
    x=np.linspace(0,g.length_m,101)
    t=np.full_like(x,.006)
    Ts,Ti=temperature_model(t,x,c,material_library()[1])
    assert np.all(Ti<Ts)
    assert np.all(Ti>c.cabin_C)


def test_optimized_geometry_export():
    cond=Conditions(samples=90,knots=6)
    res=optimize_one(Geometry(),cond,material_library()[1])
    assert res['feasible']
    import trimesh
    mesh=trimesh.load_mesh(file_obj=__import__('io').BytesIO(stl_bytes(res)),file_type='stl')
    assert mesh.is_watertight
    assert mesh.volume>0
    assert mesh.bounds[1,2]-mesh.bounds[0,2]>.35*1000


def test_step_contains_brep():
    pytest.importorskip('cadquery')
    cond=Conditions(samples=60,knots=5)
    res=optimize_one(Geometry(),cond,material_library()[1])
    buf=step_bytes(res)
    assert b'ISO-10303-21' in buf[:100]
    assert len(buf)>10_000
