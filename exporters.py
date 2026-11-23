"""STL / STEP export of concept radome shell; dimensions stored in mm for CAD interoperability."""
from __future__ import annotations
from io import BytesIO
import numpy as np


def _cross_section_mm(result):
    x=np.asarray(result["x_m"],float)
    r=np.asarray(result["outer_radius_m"],float)
    t=np.asarray(result["thickness_m"],float)
    slope=np.asarray(result["slope"],float)
    if min(t)<=0 or max(t)>=max(r):
        raise ValueError("Cannot export: invalid thickness relative to outer radius")
    # Offset by prescribed thickness along inner normal, then remap to common axial stations.
    norm=np.sqrt(1+slope*slope)
    xi=x+t*slope/norm
    ri=r-t/norm
    # Large curvature relative to wall thickness can produce an invalid offset shell.
    if np.min(np.diff(xi))<=0:
        raise ValueError("Normal-offset curve folds over itself; decrease thickness or change nose shape")
    if ri[-1]<=0:
        raise ValueError("No base aperture remains at selected thickness")
    # Interior cavity begins at the axial intercept of the normal-offset inner curve.
    idx=np.flatnonzero(ri>1e-10)
    if len(idx)<2: raise ValueError("Shell has no valid hollow interior")
    k=int(idx[0]); first_i=max(k-1,0)
    if ri[first_i]<0:
        x0=xi[first_i]+(0-ri[first_i])*(xi[k]-xi[first_i])/(ri[k]-ri[first_i])
    else: x0=xi[first_i]
    # Remove offset curve beyond the base and end cavity precisely at the base plane.
    target=x[(x>x0)&(x<x[-1])]
    ir=np.interp(target,xi,ri)
    base_r=float(np.interp(x[-1],xi,ri))
    if not np.isfinite(base_r) or base_r<=0: raise ValueError("Invalid inside base diameter")
    outer=list(zip(r*1000,x*1000))
    inside=[(base_r*1000,x[-1]*1000)]
    inside.extend((float(rr*1000),float(xx*1000)) for rr,xx in zip(ir[::-1],target[::-1]))
    inside.append((0.0,float(x0*1000)))
    section=np.asarray(outer+inside,float)
    if np.any(section[:,0]<-1e-8): raise ValueError("Negative interior radius, invalid CAD")
    # Outline must be simple to support solid revolve.
    try:
        from shapely.geometry import Polygon
        poly=Polygon(section)
        if not poly.is_valid or poly.area<=1e-6:
            raise ValueError("Self-intersecting wall section; use smaller thickness or another profile")
    except ImportError: pass
    return section


def stl_bytes(result, circumferential_segments=96):
    import trimesh
    section=_cross_section_mm(result)
    nang=int(circumferential_segments)
    if nang<12: raise ValueError("Need at least 12 angular segments")
    angle=2*np.pi*np.arange(nang)/nang
    vertices=[]
    for rad,z in section:
        vertices.extend((rad*np.cos(a),rad*np.sin(a),z) for a in angle)
    vertices=np.asarray(vertices,dtype=float)
    tris=[]
    m=len(section)
    for i in range(m):
        j=(i+1)%m
        r0,r1=section[i,0],section[j,0]
        for a in range(nang):
            b=(a+1)%nang
            ia=i*nang+a; ib=i*nang+b
            ja=j*nang+a; jb=j*nang+b
            if r0>1e-9 and r1>1e-9:
                tris.extend(((ia,ib,jb),(ia,jb,ja)))
            elif r0<=1e-9 and r1>1e-9:
                tris.append((ia,jb,ja))
            elif r0>1e-9 and r1<=1e-9:
                tris.append((ia,ib,ja))
    mesh=trimesh.Trimesh(vertices=vertices,faces=np.asarray(tris,int),process=True)
    mesh.remove_unreferenced_vertices()
    if mesh.volume<0: mesh.invert()
    if not mesh.is_watertight:
        raise ValueError("Mesh is not watertight; try more sections or smaller wall thickness")
    return mesh.export(file_type="stl")


def step_bytes(result):
    """CadQuery optional locally; returns STEP binary for full solid rotational shell."""
    try:
        import cadquery as cq
    except ImportError as exc:
        raise RuntimeError("STEP export requires cadquery (see README); STL remains available") from exc
    import tempfile
    from pathlib import Path
    section=_cross_section_mm(result)
    # XZ workplane: 2D x = global radial X, 2D y = global longitudinal Z.
    points=[(float(r),float(z)) for r,z in section]
    body=cq.Workplane("XZ").polyline(points).close().revolve(360,(0,0),(0,1))
    with tempfile.TemporaryDirectory() as tmp:
        path=Path(tmp)/"radome.step"
        cq.exporters.export(body,str(path),exportType="STEP")
        output=path.read_bytes()
    if len(output)<1000: raise ValueError("STEP file appears incomplete")
    return output
