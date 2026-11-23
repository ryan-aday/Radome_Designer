"""Streamlit interactive radome concept designer.
Run: streamlit run app.py
"""
from __future__ import annotations
import json
from dataclasses import asdict
from io import StringIO
from pathlib import Path
import numpy as np
import pandas as pd
import plotly.graph_objects as go
import streamlit as st
from core import (Geometry,Conditions,PROFILES,material_library,optimize_candidates,
                  alternate_profiles,outer_radius,compile_load)
from exporters import stl_bytes,step_bytes

st.set_page_config(page_title="Radome Studio · concept design",page_icon="◉",layout="wide")
st.markdown("""<style>
.block-container{padding-top:1.25rem;max-width:1450px}
h1{letter-spacing:-.035em} .stMetric {background:rgba(128,128,128,.07);padding:12px;border-radius:8px}
</style>""", unsafe_allow_html=True)

st.title("Radome Studio")
st.caption("Monolithic RF shell · preliminary RF/thermal/mechanical screening · parametric CAD")
st.warning("**Concept-screening tool, not a validated design solver.** Shell membrane, thermal stress and plane-slab RF models omit local buckling, transient heating, impact, erosion, full-wave antenna effects, boresight error and detailed attachments. Material values are illustrative—replace with qualified data before use.")

mats=material_library()
with st.sidebar.expander("Material database · JSON",expanded=False):
    st.caption("Values are illustrative and require grade-specific qualification. Edit materials.json or upload a revised array of records.")
    upload=st.file_uploader("Custom materials.json",type=["json"])
    if upload is not None:
        try:
            edited=json.loads(upload.getvalue().decode("utf-8"))
            expected=set(mats[0].keys())
            if not isinstance(edited,list) or not edited or any(not isinstance(m,dict) or not expected.issubset(m) for m in edited):
                raise ValueError("Expected a nonempty array of materials with the same keys as materials.json")
            for mat in edited:
                for field in ["density_kg_m3","epsilon_r","loss_tangent","max_service_C","strength_MPa","modulus_GPa","poisson","cte_1_per_K","thermal_k_W_mK","strength_fraction_at_limit","eps_temp_coeff_per_K","loss_temp_coeff_per_K"]:
                    if not np.isfinite(float(mat[field])): raise ValueError(f"Nonfinite {field}")
                if mat["epsilon_r"]<1 or mat["density_kg_m3"]<=0 or mat["thermal_k_W_mK"]<=0 or mat["strength_MPa"]<=0 or not 0<mat["poisson"]<.5:
                    raise ValueError("Material has nonphysical parameters")
            mats=edited
            st.success(f"Loaded {len(mats)} material records")
        except Exception as exc:
            st.error(f"Material JSON rejected: {exc}")
    st.download_button("Download database template",data=json.dumps(material_library(),indent=2).encode(),file_name="materials_template.json",mime="application/json")
mat_names=[m["name"] for m in mats]

with st.sidebar:
    st.header("Design inputs")
    with st.expander("01 · Nose geometry",expanded=True):
        st.selectbox("Profile",PROFILES,index=PROFILES.index("Tangent ogive"),key="profile_choice")
        Lmm=st.number_input("Cone length (mm)",min_value=20.0,max_value=5000.0,value=400.0,step=10.0)
        Dmm=st.number_input("Base outside diameter (mm)",min_value=10.0,max_value=3000.0,value=200.0,step=10.0)
        tip=st.number_input("Spherical tip radius (mm; 0 = unblunted)",min_value=0.,max_value=500.,value=0.,step=1.)
        join=st.slider("Spherical cap join angle (deg)",16,84,60)
        if st.session_state.profile_choice=="Parabolic": K=st.slider("Parabolic K",0.,1.,.75,.05)
        else: K=.75
        if st.session_state.profile_choice=="Power series": n=st.slider("Power exponent",.3,1.,.75,.05)
        else: n=.75
        if st.session_state.profile_choice=="Secant ogive": secant=st.slider("Circle-center / length",1.01,1.50,1.10,.01)
        else: secant=1.10
        if st.session_state.profile_choice=="Bi-conic":
            b=st.slider("Break location (fraction L)",.1,.9,.45,.05)
            br=st.slider("Radius at break (fraction R)",.1,.9,.65,.05)
        else: b,br=.45,.65
    with st.expander("02 · RF specification",expanded=True):
        freq=st.number_input("Center frequency (GHz)",min_value=.1,max_value=200.,value=10.,step=1.)
        bandwidth=st.slider("Fractional bandwidth (%)",0,80,10,5)
        angle=st.slider("Maximum incidence / scan offset (deg)",0,80,45,5)
        angle_count=st.slider("Incidence samples",2,9,5)
        mode=st.selectbox("Incidence model",["User-defined local incidence envelope","Axial-ray normal-angle proxy"])
        maxloss=st.number_input("Maximum worst-case patch loss (dB)",min_value=.05,max_value=15.,value=1.5,step=.1)
    with st.expander("03 · Temperature / environment",expanded=True):
        Ttip=st.number_input("Tip skin temperature (°C)",min_value=-30.,max_value=1800.,value=150.,step=10.)
        Tbase=st.number_input("Base skin temperature (°C)",min_value=-30.,max_value=1800.,value=50.,step=10.)
        Tcab=st.number_input("Interior air temperature (°C)",min_value=-30.,max_value=300.,value=25.,step=5.)
        h=st.number_input("Interior convection h (W/m²K)",min_value=.1,max_value=1000.,value=12.,step=2.)
        Tshape=st.slider("Skin temperature decay exponent",.2,5.,1.5,.1)
        pressure=st.number_input("Differential pressure magnitude (kPa)",min_value=0.,max_value=1000.,value=20.,step=5.)
        restraint=st.slider("Thermal restraint stress proxy factor",0.,.3,.03,.01)
        grad_stress=st.slider("Thermal-gradient stress proxy factor",0.,.3,.03,.01)
        sf=st.number_input("Structural safety factor",min_value=1.,max_value=8.,value=2.,step=.25)
    with st.expander("04 · Section resultants",expanded=False):
        st.caption("Position x from tip in meters; u=x/L (0 tip, 1 base). Output axial/shear in N, bending in N·m. e.g. `1000*u` or `50*(1-u)**2`. These are resultant loads at each station, not distributed load intensity.")
        Nexpr=st.text_input("Axial force N(x) [N]","0")
        Vexpr=st.text_input("Shear force V(x) [N]","0")
        Mexpr=st.text_input("Bending moment M(x) [N·m]","0")
    with st.expander("05 · Optimization",expanded=False):
        tmin=st.number_input("Minimum wall (mm)",min_value=.2,max_value=30.,value=1.5,step=.25)
        tmax=st.number_input("Maximum wall (mm)",min_value=.3,max_value=60.,value=14.,step=.5)
        knots=st.slider("Thickness control points",4,12,8)
        samples=st.slider("Longitudinal samples",60,240,120,20)
        w_rf=st.slider("RF loss weight",.05,5.,1.,.05)
        w_mass=st.slider("Mass weight",0.,3.,.12,.02)
        w_smooth=st.slider("Smoothness weight",0.,3.,.12,.02)
        selected=st.multiselect("Materials to consider",mat_names,default=mat_names)
    launch=st.button("Optimize selected design",type="primary",use_container_width=True)

geo=Geometry(length_m=Lmm/1000,base_diameter_m=Dmm/1000,profile=st.session_state.profile_choice,
             spherical_tip_radius_m=tip/1000,spherical_tip_join_deg=join,parabolic_k=K,power_n=n,
             secant_center_ratio=secant,bicone_break_fraction=b,bicone_radius_fraction=br)
cond=Conditions(frequency_GHz=freq,bandwidth_fraction=bandwidth/100,angle_max_deg=angle,angle_count=angle_count,
                rf_angle_mode=mode,skin_tip_C=Ttip,skin_base_C=Tbase,skin_shape_exp=Tshape,
                cabin_C=Tcab,h_inner_W_m2K=h,differential_pressure_Pa=pressure*1000,
                axial_load_expr=Nexpr,shear_load_expr=Vexpr,moment_load_expr=Mexpr,
                restraint_fraction=restraint,gradient_stress_fraction=grad_stress,safety_factor=sf,
                thickness_min_m=tmin/1000,thickness_max_m=tmax/1000,knots=knots,
                rf_weight=w_rf,mass_weight=w_mass,smoothness_weight=w_smooth,max_loss_dB=maxloss,
                samples=samples)

# Preview geometry even before optimizing.
left,right=st.columns([1.6,1])
with left:
    st.subheader("Shape preview")
    try:
        xp=np.linspace(0,geo.length_m,250); rp=outer_radius(geo,xp)
        fig=go.Figure()
        fig.add_trace(go.Scatter(x=xp*1000,y=rp*1000,mode="lines",name="Outer profile",line=dict(width=3)))
        fig.add_trace(go.Scatter(x=xp*1000,y=-rp*1000,mode="lines",showlegend=False,line=dict(width=3)))
        fig.update_layout(xaxis_title="Position from tip (mm)",yaxis_title="Radius (mm)",height=310,margin=dict(l=30,r=20,t=20,b=40),yaxis_scaleanchor="x",yaxis_scaleratio=1)
        st.plotly_chart(fig,use_container_width=True)
    except ValueError as e: st.error(str(e))
with right:
    st.subheader("Reference thickness")
    wl=299.792458/freq # mm wavelength in air
    st.metric("Free-space wavelength",f"{wl:.2f} mm")
    st.metric("Nominal dielectric half-wave",f"{(wl/(2*np.sqrt(3.9))):.2f} mm",help="Illustrative reference for εr=3.9 only; actual optimization includes materials and angles")
    st.caption("Near-zero thickness may minimize reflection but can fail mechanics. The global response is oscillatory in wall thickness; three optimization starts are used.")

if launch:
    if not selected:
        st.error("Choose at least one material")
    else:
        try:
            for exp in (Nexpr,Vexpr,Mexpr): compile_load(exp)(np.linspace(0,geo.length_m,8),geo.length_m)
            if tmax<=tmin: raise ValueError("Maximum thickness must exceed minimum")
            with st.spinner("Evaluating material and thickness candidates …"):
                output=optimize_candidates(geo,cond,[m for m in mats if m["name"] in selected])
            st.session_state.optimized={"results":output,"geo":geo,"cond":cond,"selected":selected}
        except Exception as exc:
            st.error(f"Input or optimizer failure: {exc}")

result_bundle=st.session_state.get("optimized")
if result_bundle:
    if result_bundle["geo"]!=geo or result_bundle["cond"]!=cond or result_bundle["selected"]!=selected:
        st.info("Inputs have changed. Results below belong to the last completed optimization; press Optimize to refresh.")
    results=result_bundle["results"]
    st.divider(); st.subheader("Material comparison")
    summary=[]
    for entry in results:
        summary.append({"Material":entry["material"],"Status":entry["status"],
                        "Mass (kg)":round(entry["mass_kg"],3) if "mass_kg" in entry else None,
                        "Mean loss (dB)":round(entry["em"]["mean_loss_dB"],3) if "em" in entry else None,
                        "Worst sampled loss (dB)":round(entry["em"]["worst_loss_dB"],3) if "em" in entry else None,
                        "Max stress utilization":round(entry["max_utilization"],3) if "max_utilization" in entry else None})
    st.dataframe(pd.DataFrame(summary),use_container_width=True,hide_index=True)
    candidates=[entry for entry in results if "thickness_m" in entry]
    good=[entry for entry in candidates if entry["feasible"]]
    if not good:
        st.error("No candidate met *all* specified screening criteria. The best result can still be inspected, but should not be treated as feasible.")
        if st.button("Try alternative nose profiles"):
            with st.spinner("Screening alternate profiles …"):
                alts=alternate_profiles(result_bundle["geo"],result_bundle["cond"],[m for m in mats if m["name"] in selected])
            st.session_state["alt_results"]=alts
        if "alt_results" in st.session_state:
            st.dataframe(pd.DataFrame(st.session_state.alt_results),hide_index=True,use_container_width=True)
            options=[a["profile"] for a in st.session_state.alt_results if a["score"] is not None]
            if options:
                choice=st.selectbox("Suggested alternative",options)
                if st.button("Set selected alternative as active profile"):
                    st.session_state.profile_choice=choice
                    st.rerun()
    if candidates:
        choices=[c["material"] for c in candidates]
        chosen=st.selectbox("Inspect optimized material",choices)
        res=next(r for r in candidates if r["material"]==chosen)
        metrics=st.columns(4)
        metrics[0].metric("Concept mass",f"{res['mass_kg']:.3f} kg")
        metrics[1].metric("Mean insertion loss",f"{res['em']['mean_loss_dB']:.3f} dB")
        metrics[2].metric("Worst sampled loss",f"{res['em']['worst_loss_dB']:.3f} dB")
        metrics[3].metric("Max utilization",f"{res['max_utilization']:.2f}×")
        if not res["optimizer_converged"]: st.warning("Local optimizer reported non-convergence. A feasible stationary result may still be shown; review and rerun.")
        x=np.asarray(res["x_m"])*1000
        t=np.asarray(res["thickness_m"])*1000
        mech=res["mechanics"]
        c1,c2=st.columns(2)
        with c1:
            f=go.Figure();f.add_scatter(x=x,y=t,mode="lines",name="Normal wall thickness")
            f.update_layout(title="Optimized wall thickness",xaxis_title="Axial position (mm)",yaxis_title="Thickness (mm)",height=330)
            st.plotly_chart(f,use_container_width=True)
        with c2:
            f=go.Figure();f.add_scatter(x=x,y=mech["utilization"],name="Stress utilization")
            f.add_hline(y=1,line_dash="dash")
            f.add_vrect(x0=0,x1=.06*x[-1],opacity=.12,line_width=0,annotation_text="tip excluded")
            f.update_layout(title="Structural screening (tip excluded)",xaxis_title="Axial position (mm)",yaxis_title="Utilization",height=330)
            st.plotly_chart(f,use_container_width=True)
        c3,c4=st.columns(2)
        with c3:
            f=go.Figure();f.add_scatter(x=x,y=mech["temperature_outer_C"],name="Hot skin")
            f.add_scatter(x=x,y=mech["temperature_inner_C"],name="Predicted inner surface")
            f.update_layout(title="Steady-state wall temperatures",xaxis_title="Axial position (mm)",yaxis_title="°C",height=330)
            st.plotly_chart(f,use_container_width=True)
        with c4:
            f=go.Figure()
            for i,frequency in enumerate(res["em"]["freq_GHz"]):
                f.add_scatter(x=res["em"]["angles_deg"],y=res["em"]["TE_loss_dB"][i],name=f"TE {frequency:.2f} GHz")
                f.add_scatter(x=res["em"]["angles_deg"],y=res["em"]["TM_loss_dB"][i],name=f"TM {frequency:.2f} GHz",line=dict(dash="dot"))
            f.update_layout(title="Area-averaged patch insertion loss",xaxis_title="Incidence/scan offset (deg)",yaxis_title="One-way loss (dB)",height=330)
            st.plotly_chart(f,use_container_width=True)
        with st.expander("3D shell visualization",expanded=False):
            idx=np.linspace(0,len(x)-1,min(len(x),65)).astype(int)
            phi=np.linspace(0,2*np.pi,48)
            rr=np.asarray(res["outer_radius_m"])[idx]*1000
            zz=x[idx]
            tt=t[idx]
            ss=np.asarray(res["slope"])[idx]
            ri=np.maximum(0,rr-tt/np.sqrt(1+ss*ss))
            zi=zz+tt*ss/np.sqrt(1+ss*ss)
            scene=go.Figure()
            for radii,axial,opacity,col in [(rr,zz,.95,"Viridis"),(ri,zi,.42,"Greys")]:
                scene.add_surface(x=np.outer(np.cos(phi),radii),y=np.outer(np.sin(phi),radii),
                                  z=np.tile(axial,(len(phi),1)),showscale=False,opacity=opacity,colorscale=col)
            scene.update_layout(height=510,scene=dict(xaxis_title="X (mm)",yaxis_title="Y (mm)",zaxis_title="Axial mm",aspectmode="data"),margin=dict(l=0,r=0,t=0,b=0))
            st.plotly_chart(scene,use_container_width=True)
            st.caption("Rendering is visual only. CAD export uses an interpolated normal-offset inner wall and exact base plane closure.")
        st.markdown("#### Section data and exports")
        st.caption("STL/STEP use millimeters and a **normal-offset** interior wall. Solid geometry is valid only if the resulting 2D profile can be revolved without self-intersection. No flange, mount, coating, adhesive, or tolerance is included.")
        df=pd.DataFrame({"x_mm":x,"outer_radius_mm":np.asarray(res["outer_radius_m"])*1000,
                         "normal_thickness_mm":t,"skin_C":mech["temperature_outer_C"],
                         "inner_C":mech["temperature_inner_C"],"failure_metric_MPa":mech["failure_metric_Pa"]/1e6,
                         "stress_allowable_MPa":mech["allowable_Pa"]/1e6,"utilization":mech["utilization"]})
        with st.expander("Detailed section table"):
            st.dataframe(df,use_container_width=True,hide_index=True)
        st.download_button("Download section CSV",data=df.to_csv(index=False).encode(),file_name="radome_sections.csv",mime="text/csv")
        config={"geometry":res["geo"],"conditions":res["conditions"],"material":res["material_properties"],
                "result":{"screening_pass":res["feasible"],"mass_kg":res["mass_kg"],"mean_loss_dB":res["em"]["mean_loss_dB"],
                          "worst_loss_dB":res["em"]["worst_loss_dB"],"max_utilization":res["max_utilization"]}}
        st.download_button("Download design JSON",data=json.dumps(config,indent=2),file_name="radome_design.json",mime="application/json")
        with st.expander("Generate CAD geometry",expanded=True):
            col1,col2=st.columns(2)
            with col1:
                try:
                    st.download_button("Export STL (watertight)",data=stl_bytes(res),file_name="radome_optimized.stl",mime="model/stl",use_container_width=True)
                except Exception as exc: st.warning(f"STL unavailable: {exc}")
            with col2:
                try:
                    st.download_button("Export STEP (solid)",data=step_bytes(res),file_name="radome_optimized.step",mime="application/step",use_container_width=True)
                except Exception as exc: st.warning(f"STEP unavailable: {exc}")
        if not res["feasible"]:
            st.error("This geometry is **not** a screening pass. CAD export is for inspection only.")

with st.expander("Model theory, limitations and citations",expanded=False):
    st.markdown(r"""
**RF wall:** Homogeneous nonmagnetic dielectric, complex relative permittivity $\epsilon_c=\epsilon_r(1-j\tan\delta)$, propagation constant $k_{z}=k_0\sqrt{\epsilon_c-\sin^2\theta}$, and TE/TM two-port transmission-line matrix. Air is assumed both sides, no layers/coatings. Reported power transmission $P_T=|S_{21}|^2$ and insertion loss $IL=-10\log_{10}P_T$. Area-weighted planar patches are **not** antenna patterns or boresight errors.

**Heat:** Prescribed outer skin temperature $T_s(x)$. One-dimensional steady conduction in the wall coupled to interior air by convection: $q''=(T_s-T_\infty)/(t/k+1/h_i)$. No transient, radiation, ablative recession, or spatial heat redistribution.

**Stress:** Axial and bending stress approximated by thin circular-section resultants $N/(2\pi rt)+|M|/(\pi r^2 t)$, membrane hoop pressure stress $\Delta p r/t$, and shear proxy $2|V|/(\pi rt)$; combined with *user-adjustable* thermal stress envelopes. Ductile materials use a von Mises-like metric; brittle ceramics use a maximum-principal-tension proxy. This is **not** curved-shell elasticity; pointed tip region and non-thin sections are excluded from the pass/fail screen. Brittle reliability, fracture initiation, and Weibull statistics are not modeled.

**Optimization:** PCHIP thickness interpolation between axial control points; SLSQP local multi-start under sampled structural constraints. Objective = RF mean loss + normalized mass + thickness curvature penalty. The worst sampled angle/frequency **area-averaged** loss is constrained; if no joint solution exists, the solver can display a structural-only candidate with explicit RF failure. RF response can have multiple thickness optima.

**References and PDFs:** See **References & papers** immediately below, with direct source links, bundled PDFs and the role each work plays in the current model.
""")

# The bibliography is data-driven so the main interface and bundled sources stay aligned.
st.divider()
st.subheader("References & papers")
st.caption("Four original web links (R01–R04), two original user-uploaded PDFs (R05–R06), and one supplemental ceramics review (R07). Bundled PDF downloads work offline; links require internet.")

_catalog_path = Path(__file__).resolve().parent / "references" / "catalog.json"
if _catalog_path.is_file():
    _refs = json.loads(_catalog_path.read_text(encoding="utf-8"))["references"]
    for _ref in _refs:
        with st.expander(f"{_ref['id']} · {_ref['title']} ({_ref['year']})", expanded=False):
            st.markdown(f"**Citation:** {_ref['authors']}. *{_ref['publication']}*.")
            st.markdown(f"**Source:** [{_ref['url']}]({_ref['url']})")
            if _ref.get("pdf_url") and _ref["pdf_url"] != _ref["url"]:
                st.markdown(f"**Original PDF:** [{_ref['pdf_url']}]({_ref['pdf_url']})")
            st.markdown(f"**Use:** {_ref['relevance']}")
            st.caption(f"App implementation: {_ref['implementation']}")
            if _ref.get("local_pdf"):
                _paper_path=Path(__file__).resolve().parent / _ref["local_pdf"]
                if _paper_path.is_file():
                    st.download_button(
                        f"Download bundled PDF — {_ref['id']}",
                        data=_paper_path.read_bytes(),
                        file_name=_paper_path.name,
                        mime="application/pdf",
                        key=f"ref_pdf_{_ref['id']}",
                    )
                else:
                    st.warning(f"Bundled paper is missing: {_ref['local_pdf']}")
            elif _ref.get("pdf_url"):
                st.caption("Link-only publisher PDF; not redistributed as part of the ZIP. An optional personal-use fetch script is included.")
    st.caption("Scientific references describe model foundations and future work; citing them does not imply that experimental material allowables, antenna boresight errors, wet radomes, ablation, or multilayer optimizations have been validated here.")
else:
    st.error("Reference catalog not found. Check references/catalog.json in your installation.")
