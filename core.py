"""Reduced-order civil/industrial RF radome concept screening.
SI internally. Not a structural certification, antenna solver or manufacturing release.
"""
from __future__ import annotations
import ast
import json
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Callable
import numpy as np
from scipy.integrate import trapezoid
from scipy.interpolate import PchipInterpolator
from scipy.optimize import minimize

C0 = 299_792_458.0
PROFILES = ("Conical", "Tangent ogive", "Secant ogive", "Elliptical", "Parabolic", "Power series", "Von Karman", "Haack LV", "Bi-conic")


@dataclass
class Geometry:
    length_m: float = 0.4
    base_diameter_m: float = 0.2
    profile: str = "Tangent ogive"
    parabolic_k: float = 0.75
    power_n: float = 0.75
    secant_center_ratio: float = 1.10
    bicone_break_fraction: float = 0.45
    bicone_radius_fraction: float = 0.65
    spherical_tip_radius_m: float = 0.0
    spherical_tip_join_deg: float = 60.0


@dataclass
class Conditions:
    frequency_GHz: float = 10.0
    bandwidth_fraction: float = 0.10
    angle_max_deg: float = 45.0
    angle_count: int = 5
    rf_angle_mode: str = "User-defined local incidence envelope"
    skin_tip_C: float = 150.0
    skin_base_C: float = 50.0
    skin_shape_exp: float = 1.5
    cabin_C: float = 25.0
    h_inner_W_m2K: float = 12.0
    differential_pressure_Pa: float = 20_000.0
    axial_load_expr: str = "0"
    shear_load_expr: str = "0"
    moment_load_expr: str = "0"
    restraint_fraction: float = 0.03
    gradient_stress_fraction: float = 0.03
    safety_factor: float = 2.0
    thickness_min_m: float = 0.0015
    thickness_max_m: float = 0.014
    knots: int = 8
    rf_weight: float = 1.0
    mass_weight: float = 0.12
    smoothness_weight: float = 0.12
    max_loss_dB: float = 1.5
    samples: int = 120


def material_library(path=None):
    return json.loads(Path(path or Path(__file__).with_name("materials.json")).read_text(encoding="utf-8"))


def _base_radius(geo: Geometry, x):
    x = np.asarray(x, dtype=float)
    L, R = geo.length_m, geo.base_diameter_m / 2
    if L <= 0 or R <= 0:
        raise ValueError("Length and diameter must be positive")
    u = np.clip(x / L, 0, 1)
    kind = geo.profile
    if kind == "Conical":
        y = R * u
    elif kind == "Tangent ogive":
        rho = (L*L + R*R)/(2*R)
        y = np.sqrt(np.maximum(0.0, rho*rho - (L-x)**2)) + R-rho
    elif kind == "Secant ogive":
        xc = geo.secant_center_ratio*L
        yc = (L*L+R*R-2*L*xc)/(2*R)
        if yc >= 0:
            raise ValueError("Secant center factor incompatible with this L/D: choose larger center ratio")
        rho = np.hypot(xc, yc)
        y = yc + np.sqrt(np.maximum(0.0, rho*rho-(x-xc)**2))
    elif kind == "Elliptical":
        y = R * np.sqrt(np.maximum(0.0, 2*u-u*u))
    elif kind == "Parabolic":
        K = geo.parabolic_k
        if not 0 <= K <= 1: raise ValueError("Parabolic K must be 0..1")
        y = R*(2*u-K*u*u)/(2-K)
    elif kind == "Power series":
        n=geo.power_n
        if not 0.3 <= n <= 1: raise ValueError("Power exponent must be 0.3..1")
        y = R*u**n
    elif kind in ("Von Karman", "Haack LV"):
        th = np.arccos(np.clip(1-2*u,-1,1))
        C = 0 if kind == "Von Karman" else 1/3
        y = R*np.sqrt(np.maximum(0,(th-0.5*np.sin(2*th)+C*np.sin(th)**3)/np.pi))
    elif kind == "Bi-conic":
        b, frac = geo.bicone_break_fraction, geo.bicone_radius_fraction
        if not 0.05<b<0.95 or not 0.05<frac<0.95: raise ValueError("Bi-conic fractions out of range")
        y = R*np.where(u<=b, frac*u/b, frac+(1-frac)*(u-b)/(1-b))
    else:
        raise ValueError(f"Unknown profile: {kind}")
    y[0] = 0.0 if np.isclose(x[0],0) else y[0]
    return y


def outer_radius(geo: Geometry, x):
    """Optionally splice an exact spherical cap with C1 Hermite transition to main profile.
    Geometry is a concept blending approximation; inspect monotonicity.
    """
    x = np.asarray(x, dtype=float)
    r = _base_radius(geo,x)
    S = geo.spherical_tip_radius_m
    if S <= 0: return r
    theta = np.deg2rad(geo.spherical_tip_join_deg)
    if not 15 < geo.spherical_tip_join_deg < 85:
        raise ValueError("Spherical cap join angle must lie strictly between 15 and 85 degrees")
    xj = S*(1-np.cos(theta)); yj=S*np.sin(theta)
    L, R = geo.length_m, geo.base_diameter_m/2
    if xj >= 0.3*L or yj >= 0.85*R:
        raise ValueError("Spherical cap too large for nose geometry; reduce tip radius")
    rj = float(_base_radius(geo,np.array([xj]))[0])
    dx=max(1e-7,L*1e-5)
    drj=(float(_base_radius(geo,np.array([xj+dx]))[0])-float(_base_radius(geo,np.array([xj-dx]))[0]))/(2*dx)
    # value and derivative correction, zeroing smoothly at base
    v = np.clip((x-xj)/(L-xj),0,1)
    correction = (2*v**3-3*v**2+1)*(yj-rj) + (v**3-2*v**2+v)*(L-xj)*(1/np.tan(theta)-drj)
    r = np.where(x<=xj, np.sqrt(np.maximum(0,2*S*x-x*x)), r+correction)
    return r


_ALLOWED_FUNCS = {"sin":np.sin,"cos":np.cos,"tan":np.tan,"exp":np.exp,"sqrt":np.sqrt,
                  "log":np.log,"abs":np.abs,"minimum":np.minimum,"maximum":np.maximum,"where":np.where}
_ALLOWED_VARS={"x","u","L","pi"}
_ALLOWED_AST=(ast.Expression,ast.BinOp,ast.UnaryOp,ast.Constant,ast.Name,ast.Call,
              ast.Add,ast.Sub,ast.Mult,ast.Div,ast.Pow,ast.USub,ast.UAdd,ast.Mod,
              ast.Compare,ast.Lt,ast.LtE,ast.Gt,ast.GtE,ast.Eq,ast.NotEq,ast.BoolOp,ast.And,ast.Or,ast.Load)


def compile_load(expr: str) -> Callable:
    """Safe arithmetic parser; expressions in x [m], u=x/L, L [m], output SI."""
    root=ast.parse(expr,mode="eval")
    nodes=list(ast.walk(root))
    if len(nodes)>100: raise ValueError("Load expression is too long")
    for node in nodes:
        if not isinstance(node,_ALLOWED_AST): raise ValueError(f"Disallowed expression feature: {type(node).__name__}")
        if isinstance(node,ast.Call) and (not isinstance(node.func,ast.Name) or node.func.id not in _ALLOWED_FUNCS):
            raise ValueError("Only documented mathematical functions are allowed")
        if isinstance(node,ast.Name) and node.id not in _ALLOWED_VARS and node.id not in _ALLOWED_FUNCS:
            raise ValueError(f"Unknown symbol: {node.id}")
        if isinstance(node,ast.Constant) and not isinstance(node.value,(int,float,bool)):
            raise ValueError("Only numeric constants permitted")
    compiled=compile(root,"<load_expression>","eval")
    def fn(x,L):
        x=np.asarray(x,float)
        val=eval(compiled,{"__builtins__":{}},{"x":x,"u":x/L,"L":L,"pi":np.pi,**_ALLOWED_FUNCS})
        a=np.broadcast_to(np.asarray(val,dtype=float),x.shape)
        if not np.all(np.isfinite(a)): raise ValueError("Load expression produced nonfinite values")
        return a
    return fn


def surface_grid(geo:Geometry, cond:Conditions):
    x=geo.length_m*np.linspace(0,1,cond.samples)
    r=outer_radius(geo,x)
    if np.any(~np.isfinite(r)) or np.any(r < -1e-8) or np.any(np.diff(r)<-1e-6):
        raise ValueError("Non-monotonic or invalid shell profile; try another profile, smaller cap or a smaller blend angle")
    r=np.maximum(0,r)
    slope=np.gradient(r,x,edge_order=2)
    arc=np.sqrt(1+slope*slope)
    area_weights=2*np.pi*np.maximum(r,0)*arc
    return x,r,slope,area_weights


def temperature_model(t,x,cond:Conditions,mat):
    """Steady 1D planar through-thickness heat flow against cabin convection; no ablation/radiation."""
    u=x/x[-1]
    outer=cond.skin_base_C+(cond.skin_tip_C-cond.skin_base_C)*(1-u)**cond.skin_shape_exp
    k=float(mat["thermal_k_W_mK"])
    h=max(cond.h_inner_W_m2K,1e-6)
    inner=(outer*(k/t)+cond.cabin_C*h)/(k/t+h)
    return outer,inner


def local_allowable(mat,temp_C,cond:Conditions):
    frac=np.clip((temp_C-20)/max(mat["max_service_C"]-20,1),0,1)
    factor=1-(1-mat["strength_fraction_at_limit"])*frac
    return mat["strength_MPa"]*1e6*factor/cond.safety_factor


def mechanical_metrics(x,r,t,cond:Conditions,mat):
    outer_T,inner_T=temperature_model(t,x,cond,mat)
    mid=(outer_T+inner_T)/2
    hot=np.maximum(outer_T,inner_T)
    if np.any(hot > mat["max_service_C"]):
        # feasibility handled separately
        pass
    N=compile_load(cond.axial_load_expr)(x,x[-1])
    V=compile_load(cond.shear_load_expr)(x,x[-1])
    M=compile_load(cond.moment_load_expr)(x,x[-1])
    # Tube membrane proxies, NOT a shell solution, meaningless at pointed tip.
    effective_radius=np.maximum(r,0.015*cond.thickness_max_m)
    sx_mech=N/(2*np.pi*effective_radius*t) + np.abs(M)/(np.pi*effective_radius**2*t)
    hoop=np.abs(cond.differential_pressure_Pa)*effective_radius/t
    tau=2*np.abs(V)/(np.pi*effective_radius*t)
    # Very approximate restrained + through-wall-gradient thermal stress envelopes.
    E=mat["modulus_GPa"]*1e9
    alpha=mat["cte_1_per_K"]
    restraint=cond.restraint_fraction*np.maximum(mid-cond.cabin_C,0)
    gradient=cond.gradient_stress_fraction*np.abs(outer_T-inner_T)
    sigma_T=E*alpha*(restraint+gradient)/max(1-mat["poisson"],.05)
    sx=sx_mech+sigma_T
    vm=np.sqrt(np.maximum(0,sx*sx+hoop*hoop-sx*hoop+3*tau*tau))
    if mat.get("is_brittle",False):
        # Positive maximum principal stress proxy for brittle ceramics (not a Weibull reliability model).
        equivalent=0.5*(sx+hoop)+np.sqrt(((sx-hoop)*0.5)**2+tau*tau)
    else:
        equivalent=vm
    allow=local_allowable(mat,hot,cond)
    utilization=equivalent/np.maximum(allow,1)
    return {"temperature_outer_C":outer_T,"temperature_inner_C":inner_T,
            "failure_metric_Pa":equivalent, "von_mises_proxy_Pa":vm,"allowable_Pa":allow,"utilization":utilization,
            "axial_N":N,"shear_N":V,"moment_Nm":M}


def _slab_power(er,loss,t,freq_Hz,inc_angle_deg,pol):
    """Air-dielectric-air complex transfer matrix, power transmission, TE/TM.
    Broadcast t/frequency/angle to common shape.
    """
    sin2=np.sin(np.deg2rad(inc_angle_deg))**2
    e=er*(1-1j*loss)
    kz=np.sqrt(e-sin2+0j)
    q0=np.sqrt(1-sin2+0j)
    q=kz if pol=="TE" else e/kz
    qair=q0 if pol=="TE" else 1/q0
    phase=2*np.pi*freq_Hz/C0*kz*t
    A=np.cos(phase); B=1j*np.sin(phase)/q
    C=1j*q*np.sin(phase); D=A
    denom=qair*A+qair*qair*B+C+qair*D
    amp=2*qair/denom
    power=np.abs(amp)**2
    return np.clip(np.real(power),1e-12,1)


def em_metrics(t,x,slope,area_w,cond:Conditions,mat):
    # Samples stations, frequencies, and incidence angles; local planar homogeneous wall only.
    nfreq=3
    freq=cond.frequency_GHz*1e9*np.linspace(1-cond.bandwidth_fraction/2,1+cond.bandwidth_fraction/2,nfreq)
    angles=np.linspace(0,cond.angle_max_deg,int(cond.angle_count))
    if cond.rf_angle_mode=="Axial-ray normal-angle proxy":
        intrinsic=np.rad2deg(np.arctan2(1,np.maximum(np.abs(slope),1e-6)))
        angles=np.clip(intrinsic[None,:]+np.linspace(0,cond.angle_max_deg,int(cond.angle_count))[:,None],0,85)
    else:
        angles=angles[:,None]+np.zeros((1,len(x)))
    freq=freq[:,None,None]
    ang=angles[None,:,:]
    thick=t[None,None,:]
    # Mid-wall temperature derivative model. Coefficients illustrative and require test data.
    t_outer,t_inner=temperature_model(t,x,cond,mat)
    midT=(t_outer+t_inner)/2
    eps=mat["epsilon_r"]*(1+mat["eps_temp_coeff_per_K"]*(midT-20))
    tand=mat["loss_tangent"]*np.maximum(0.05,1+mat["loss_temp_coeff_per_K"]*(midT-20))
    eps=eps[None,None,:]; tand=tand[None,None,:]
    pte=_slab_power(eps,tand,thick,freq,ang,"TE")
    ptm=_slab_power(eps,tand,thick,freq,ang,"TM")
    # Beam/ray distribution unknown. Surface area weight is a screening assumption, NOT antenna gain.
    weights=area_w/np.maximum(np.sum(area_w),1e-20)
    avgte=np.sum(pte*weights,axis=-1)
    avgtm=np.sum(ptm*weights,axis=-1)
    ilte=-10*np.log10(np.maximum(avgte,1e-12))
    iltm=-10*np.log10(np.maximum(avgtm,1e-12))
    return {"mean_loss_dB":float(0.5*(ilte.mean()+iltm.mean())),
            "worst_loss_dB":float(max(ilte.max(),iltm.max())),
            "mean_power_transmission":float(0.5*(avgte.mean()+avgtm.mean())),
            "TE_loss_dB":ilte,"TM_loss_dB":iltm,
            "freq_GHz":freq[:,0,0]/1e9,"angles_deg":np.linspace(0,cond.angle_max_deg,int(cond.angle_count))}


def mass_integral(x,r,t,slope,mat):
    # Thin shell midsurface mass approximation, t = normal thickness
    area=2*np.pi*r*np.sqrt(1+slope*slope)
    return float(trapezoid(t*area,x)*mat["density_kg_m3"])


def _evaluate_thickness(t,geo,cond,mat,struct_mask=None, pre=None):
    x,r,slope,area=pre if pre is not None else surface_grid(geo,cond)
    if struct_mask is None: struct_mask=x>=0.06*geo.length_m
    mech=mechanical_metrics(x,r,t,cond,mat)
    util=np.max(mech["utilization"][struct_mask])
    em=em_metrics(t,x,slope,area,cond,mat)
    mass=mass_integral(x,r,t,slope,mat)
    return em,mech,mass,util


def optimize_one(geo: Geometry,cond:Conditions,mat):
    pre=surface_grid(geo,cond)
    x,r,slope,area=pre
    if cond.skin_tip_C>mat["max_service_C"] or cond.skin_base_C>mat["max_service_C"]:
        return {"material":mat["name"],"status":"Rejected: maximum material service temperature exceeded","feasible":False}
    if cond.thickness_min_m<=0 or cond.thickness_max_m<=cond.thickness_min_m:
        raise ValueError("Invalid thickness bounds")
    if not 3<=cond.knots<=16: raise ValueError("Control-point count must be 3..16")
    mask=(x>=.06*geo.length_m) & (r > 2*cond.thickness_max_m)
    if not np.any(mask):
        return {"material":mat["name"],"status":"Rejected: no region where thin-shell stress proxy is applicable","feasible":False}
    knot_x=np.linspace(0,geo.length_m,cond.knots)
    def expand(knots):
        return np.asarray(PchipInterpolator(knot_x,knots)(x))
    def structural_slack(knots):
        t=expand(knots)
        mech=mechanical_metrics(x,r,t,cond,mat)
        # Reserve negative utilization for pass; interpolate conservative check on all sampled stations.
        return 1-mech["utilization"][mask]
    lo,hi=cond.thickness_min_m,cond.thickness_max_m
    if np.min(structural_slack(np.full(cond.knots,hi))) < -1e-5:
        return {"material":mat["name"],"status":"Rejected: maximum thickness fails structural screening","feasible":False}
    # Find smallest feasible uniform seed on a grid to reduce sensitivity to first guess.
    seed_uniform=hi
    for candidate in np.linspace(lo,hi,24):
        if np.min(structural_slack(np.full(cond.knots,candidate)))>=0:
            seed_uniform=float(candidate)
            break
    def objective(knots):
        t=expand(knots)
        em=em_metrics(t,x,slope,area,cond,mat)
        m=mass_integral(x,r,t,slope,mat)
        relative_mass=m/(2200.0*2*np.pi*geo.length_m*geo.base_diameter_m/2*hi)
        rough=np.mean(np.diff(knots,2)**2)/(hi-lo)**2
        return cond.rf_weight*em["mean_loss_dB"]+cond.mass_weight*relative_mass+cond.smoothness_weight*rough
    def rf_slack(knots):
        t=expand(knots)
        return cond.max_loss_dB-em_metrics(t,x,slope,area,cond,mat)["worst_loss_dB"]
    # Multi-start includes feasible constant seeds and halfwave-ish seed where permissible.
    halfwave=C0/(2*cond.frequency_GHz*1e9*np.sqrt(mat["epsilon_r"]))
    starts=[seed_uniform,hi,min(hi,max(seed_uniform,halfwave))]
    best=None
    # First attempt to meet both structural and RF constraints; if impossible,
    # return the best structural design labeled with RF screening failure.
    for require_rf in (True,False):
        for seed in starts:
            initial=np.full(cond.knots,seed)
            constraints=[{"type":"ineq","fun":structural_slack}]
            if require_rf: constraints.append({"type":"ineq","fun":rf_slack})
            result=minimize(objective,initial,method="SLSQP",bounds=[(lo,hi)]*cond.knots,
                            constraints=constraints,options={"maxiter":110,"ftol":4e-7})
            trial=result.x
            if np.min(structural_slack(trial))>=-1e-3 and (not require_rf or rf_slack(trial)>=-1e-3):
                score=objective(trial)
                if best is None or score<best[0]: best=(score,trial,bool(result.success),result.message)
        if best is not None: break
    if best is None:
        return {"material":mat["name"],"status":"No feasible optimized thickness profile found","feasible":False}
    t=expand(best[1]); em,mech,mass,util=_evaluate_thickness(t,geo,cond,mat,mask,pre)
    # Explicit RF cap is a screening threshold, not imposed as optimizer constraint.
    feasible=util<=1.005 and em["worst_loss_dB"]<=cond.max_loss_dB
    reasons=[]
    if util>1.005: reasons.append("structural screen not met")
    if em["worst_loss_dB"]>cond.max_loss_dB: reasons.append(f"RF loss > {cond.max_loss_dB:g} dB")
    status="Screening pass" if feasible else "Screening limit: "+", ".join(reasons)
    return {"material":mat["name"],"status":status,"feasible":feasible,
            "score":float(best[0]),"optimizer_converged":best[2],"optimizer_message":str(best[3]),
            "x_m":x,"outer_radius_m":r,"thickness_m":t,"slope":slope,
            "mass_kg":mass,"max_utilization":float(util),"em":em,"mechanics":mech,
            "geo":asdict(geo),"conditions":asdict(cond),"material_properties":mat}


def optimize_candidates(geo,cond,materials):
    results=[]
    for mat in materials:
        try: results.append(optimize_one(geo,cond,mat))
        except (ValueError,OverflowError,FloatingPointError) as exc:
            results.append({"material":mat["name"],"status":f"Error: {exc}","feasible":False})
    results.sort(key=lambda v:(not v.get("feasible",False),v.get("score",float('inf'))))
    return results


def alternate_profiles(geo,cond,materials):
    """Small shape shortlist for infeasible requested geometry; does not modify chosen design."""
    outcomes=[]
    for name in PROFILES:
        if name==geo.profile: continue
        test=Geometry(**{**asdict(geo),"profile":name})
        quick_cond=Conditions(**{**asdict(cond),"samples":60,"knots":5})
        res=optimize_candidates(test,quick_cond,materials)
        best=next((a for a in res if a.get("feasible")),None)
        outcomes.append({"profile":name,"status":"Candidate screening pass" if best else "No screening pass",
                         "material":best["material"] if best else None,"score":best["score"] if best else None})
    return sorted(outcomes,key=lambda r:(r["score"] is None,r["score"] if r["score"] is not None else 1e9))
