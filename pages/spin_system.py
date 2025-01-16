import numpy as np
import streamlit as st
from mrsimulator import Site, Coupling
from mrsimulator.spin_system.tensors import SymmetricTensor
import os
import pandas as pd


# Nuclei listing and adding chemical shift/quadrupolar interaction etc.
def add_site(list_of_nuclei, key):
    isotope = st.text_input ( f"Nucleus : ", value="1H", key = f'isotope_{key}')
    isotropic_chemical_shift = st.number_input(
        "Isotropic chemical shift (ppm)",
        value=0.0,
        format="%.2f", key = f'chemical_shift_{key}'
    )

    st.markdown("##### CSA Parameters (optional)")
    colcs1, colcs2, colcs3, colcs4, colcs5 = st.columns(5)
    with colcs1:
        zeta = st.number_input("ζ (ppm)", value=0.0, format="%.2f", key = f'cszeta_{key}')
    with colcs2:
        eta_shielding = st.number_input("η ( 0 ≤ η ≤ 1)", min_value=0.0, max_value=1.0, value=0.0, format="%.1f", key = f'cseta_{key}')
    with colcs3:
        alpha_shielding = st.number_input("α (degrees)", min_value=0., max_value=360., value=0., format="%.1f", key = f'csalpha_{key}')
    with colcs4:
        beta_shielding = st.number_input("β (degrees)", min_value=0., max_value=360., value=0., format="%.1f", key = f'csbeta_{key}')
    with colcs5:
        gamma_shielding = st.number_input("γ (degrees)", min_value=0., max_value=360., value=0., format="%.1f", key = f'csgamma_{key}')

    shielding_symmetric = SymmetricTensor(
        zeta=zeta,
        eta=eta_shielding,
        alpha=alpha_shielding*np.pi/180.0,
        beta=beta_shielding*np.pi/180.0,
        gamma=gamma_shielding*np.pi/180.0
    )

    quadrupolar = None
    mask = list_of_nuclei["Name"].isin ([isotope])
    spin_num = list_of_nuclei[ mask]["Spin"]
    if spin_num.gt ( 0.5 ).any ():
        # Section: Quadrupolar tensor
        st.markdown("##### Quadrupolar Tensor Parameters")
        col1, col2, col3, col4, col5 = st.columns(5)
        with col1:
            Cq = st.number_input("Cq (Hz)", value=0, format="%.3e", key = f'cq_{key}')
        with col2:
            eta_quadrupolar = st.number_input("η (0 ≤ η ≤ 1)", min_value=0.0, max_value=1.0, value=0., key = f'ceta_{key}')
        with col3:
            alpha_quadrupolar = st.number_input("α (degree)", min_value=0.0, max_value=360.0, value = 0.0, format="%.2f", key = f'calpha_{key}')
        with col4:
            beta_quadrupolar = st.number_input("β (degree)", min_value=0.0, max_value=360.0, value = 0.0, format="%.2f", key = f'cbeta_{key}')
        with col5:
            gamma_quadrupolar = st.number_input("γ (degree)", min_value=0.0, max_value=360.0, value = 0.0, format="%.2f", key = f'cgamma_{key}')

        quadrupolar = SymmetricTensor(
            Cq=Cq,
            eta=eta_quadrupolar,
            alpha=alpha_quadrupolar*np.pi/180.0,
            beta=beta_quadrupolar*np.pi/180.0,
            gamma=gamma_quadrupolar*np.pi/180.0
        )

    site_instance = Site(
        isotope=isotope,
        isotropic_chemical_shift=isotropic_chemical_shift,
        shielding_symmetric=shielding_symmetric,
    )


    if quadrupolar is not None:
        site_instance = Site(
            isotope=isotope,
            isotropic_chemical_shift=isotropic_chemical_shift,
            shielding_symmetric=shielding_symmetric,
            quadrupolar=quadrupolar
        )

    return site_instance


# Dipolar and J-Coupling generator
def add_couplings_DJ(num_sites, key):
    st.subheader(f"Define Advanced Coupling for Interaction Pair {key + 1}")

    col1, col2 = st.columns(2)
    with col1:
        spin1 = st.selectbox(
            "Select Site 1 (Index)", options=list(range(num_sites)), key=f"spin1_{key}", index = 0,
        )
    with col2:
        spin2 = st.selectbox(
            "Select Site 2 (Index)", options=list(range(num_sites)), key=f"spin2_{key}", index=1,
        )
    st.markdown("##### Define J-Coupling Parameters (Optional)")
    isotropic_j = st.number_input(
        "Isotropic J-Coupling (Hz)", value=10.0, format="%.3f", key=f"isotropic_j_{key}"
    )

    # Additional parameters for J-symmetric tensor
    st.markdown("##### Define J-Symmetric Tensor Parameters (Optional)")
    col_js1, col_js2, col_js3, col_js4, col_js5 = st.columns(5)
    with col_js1:
        zeta_j = st.number_input("ζ (Hz)", value=0.0, format="%.3f", key=f"zeta_j_{key}")
    with col_js2:
        eta_j = st.number_input("η", min_value=0.0, max_value=1.0, value=0.0, format="%.2f", key=f"eta_j_{key}")
    with col_js3:
        alpha_j = st.number_input("α (radians)", min_value=0.0, max_value=6.283, value=0.0, format="%.3f",
                                  key=f"alpha_j_{key}")
    with col_js4:
        beta_j = st.number_input("β (radians)", min_value=0.0, max_value=6.283, value=0.0, format="%.3f",
                                 key=f"beta_j_{key}")
    with col_js5:
        gamma_j = st.number_input("γ (radians)", min_value=0.0, max_value=6.283, value=0.0, format="%.3f",
                                  key=f"gamma_j_{key}")

    j_symmetric_tensor = SymmetricTensor(
        zeta=zeta_j, eta=eta_j, alpha=alpha_j, beta=beta_j, gamma=gamma_j
    )

    # Additional parameters for Dipolar coupling
    st.markdown("##### Define Dipolar Tensor Parameters (Optional)")
    col_ds1, col_ds2, col_ds3, col_ds4 = st.columns(4)
    with col_ds1:
        D_dipolar = st.number_input("D (Hz)", value=0.0, format="%.3f", key=f"D_dipolar_{key}")
    with col_ds2:
        alpha_d = st.number_input("α (radians)", min_value=0.0, max_value=6.283, value=0.0, format="%.3f",
                                  key=f"alpha_d_{key}")
    with col_ds3:
        beta_d = st.number_input("β (radians)", min_value=0.0, max_value=6.283, value=0.0, format="%.3f",
                                 key=f"beta_d_{key}")
    with col_ds4:
        gamma_d = st.number_input("γ (radians)", min_value=0.0, max_value=6.283, value=0.0, format="%.3f",
                                  key=f"gamma_d_{key}")

    dipolar_tensor = SymmetricTensor(
        D=D_dipolar, alpha=alpha_d, beta=beta_d, gamma=gamma_d
    )

    # Create the Coupling instance
    coupling_instance = Coupling(
        site_index=[spin1, spin2],
        isotropic_j=isotropic_j,
        j_symmetric=j_symmetric_tensor,
        dipolar=dipolar_tensor,
    )

    return coupling_instance




script_dir = os.path.dirname(__file__)
csv_file = os.path.join(script_dir, '../resources/NMR_freq_table.csv')

table_of_nuclei = pd.read_csv(csv_file)
st.header("**Add nuclei**")
sites = []


add_nuclei = "yes"
num_site = 0
while add_nuclei == "yes":
    site_added = add_site(table_of_nuclei, num_site)
    sites.append(site_added)
    num_site += 1
    st.write(f"Number of sites added: {num_site}")
    add_nuclei = st.selectbox("Add more sites?", ["yes", "no"], index = None, key = f'add_nuclei_{num_site}')
    st.divider()

if "sites" in st.session_state :
    pass
else :
    st.session_state.sites = sites

st.subheader("Summary of Spin System Parameters")

# Check if spin sites have been added
if sites:
    for idx, site in enumerate(sites, start=1):
        st.write(f"**Site {idx}:**")
        st.write(f"- **Nucleus**: {site.isotope} \t -- **σ (ppm)**: {site.isotropic_chemical_shift} ")

        # Display Shielding Symmetric Tensor Parameters
        if site.shielding_symmetric is not None:
            st.write("**Shielding Tensor Parameters:**")
            st.write(f"   **ζ (ppm)**: {site.shielding_symmetric.zeta} -- **η**: {site.shielding_symmetric.eta} -- **α (deg)**: {np.degrees(site.shielding_symmetric.alpha)} -- **β (deg)**: {np.degrees(site.shielding_symmetric.beta)} -- **γ (deg)**: {np.degrees(site.shielding_symmetric.gamma)}")

        # Display Quadrupolar Tensor Parameters if available
        if site.quadrupolar is not None:
            st.write("**Quadrupolar Tensor Parameters:**")
            st.markdown(f"  **Cq (Hz):** {site.quadrupolar.Cq} -- **η**: {site.quadrupolar.eta}  -- **α (deg)**: {np.degrees(site.quadrupolar.alpha)} -- **β (deg)**: {np.degrees(site.quadrupolar.beta)} -- **γ (deg)**: {np.degrees(site.quadrupolar.gamma)}")


        st.divider()  # Horizontal line for better readability

else:
    st.write("No spin sites have been added yet.")









# Adding J- and/or Dipolar Coupling Instances


dipolar_and_J_couplings = []

add_couplings_DJ_pairs = "yes"
num_advanced_coupling = 0
if len ( sites ) >= 2 :
    st.subheader("Add J- or Dipolar Coupling Between Spin Sites")
    while add_couplings_DJ_pairs == "yes":

        advanced_coupling_added = add_couplings_DJ(len(sites), num_advanced_coupling)
        dipolar_and_J_couplings.append(advanced_coupling_added)
        num_advanced_coupling += 1

        st.write(f"Number of advanced couplings defined: {num_advanced_coupling}")
        add_couplings_DJ_pairs = st.selectbox(
            "Add more couplings?", ["yes", "no"], index=None, key=f"add_couplings_DJ_{num_advanced_coupling}"
        )
        st.divider()

# Showing Summary of Advanced Couplings
st.subheader("Summary of Couplings")
if dipolar_and_J_couplings:
    for idx, coupling in enumerate(dipolar_and_J_couplings, start=1):
        st.write(f"**Advanced Coupling {idx}:**")
        st.write(f"- **Site 1 Index**: {coupling.site_index[0]}")
        st.write(f"- **Site 2 Index**: {coupling.site_index[1]}")
        st.write(f"- **Isotropic J-Coupling (Hz)**: {coupling.isotropic_j}")
        st.write("**J-Symmetric Tensor Parameters:**")
        st.write(f"  - ζ: {coupling.j_symmetric.zeta}")
        st.write(f"  - η: {coupling.j_symmetric.eta}")
        st.write(f"  - α: {coupling.j_symmetric.alpha}")
        st.write(f"  - β: {coupling.j_symmetric.beta}")
        st.write(f"  - γ: {coupling.j_symmetric.gamma}")
        st.write("**Dipolar Tensor Parameters:**")
        st.write(f"  - D: {coupling.dipolar.D}")
        st.write(f"  - α: {coupling.dipolar.alpha}")
        st.write(f"  - β: {coupling.dipolar.beta}")
        st.write(f"  - γ: {coupling.dipolar.gamma}")
        st.divider()
    st.session_state.dipolar_and_J_couplings = dipolar_and_J_couplings
else:
    st.write("No advanced couplings have been defined yet.")

