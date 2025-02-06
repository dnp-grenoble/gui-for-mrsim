import numpy as np
import streamlit as st
from mrsimulator import Site, Coupling, SpinSystem
from mrsimulator.spin_system.tensors import SymmetricTensor
import os
import pandas as pd
from mrsimulator.method.lib import (BlochDecaySpectrum, BlochDecayCTSpectrum, ThreeQ_VAS,
                                    FiveQ_VAS, SevenQ_VAS, ST1_VAS, ST2_VAS, SSB2D)
from mrsimulator.method import SpectralDimension
from mrsimulator import Simulator
from mrsimulator import signal_processor as sp
import time
import plotly.graph_objects as go



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
            cq_value = st.number_input("cq_value (MHz)", value=0.0, format="%.2f", key = f'cq_{key}')
        with col2:
            eta_quadrupolar = st.number_input("η (0 ≤ η ≤ 1)", min_value=0.0, max_value=1.0, value=0., key = f'ceta_{key}')
        with col3:
            alpha_quadrupolar = st.number_input("α (degree)", min_value=0.0, max_value=360.0, value = 0.0, format="%.2f", key = f'calpha_{key}')
        with col4:
            beta_quadrupolar = st.number_input("β (degree)", min_value=0.0, max_value=360.0, value = 0.0, format="%.2f", key = f'cbeta_{key}')
        with col5:
            gamma_quadrupolar = st.number_input("γ (degree)", min_value=0.0, max_value=360.0, value = 0.0, format="%.2f", key = f'cgamma_{key}')

        quadrupolar = SymmetricTensor(
            Cq=cq_value*1e6,
            eta=eta_quadrupolar,
            alpha=alpha_quadrupolar*np.pi/180.0,
            beta=beta_quadrupolar*np.pi/180.0,
            gamma=gamma_quadrupolar*np.pi/180.0
        )

    if quadrupolar is None:
        site_instance = Site(
            name=f"site-{isotope}-{key}",
            isotope=isotope,
            isotropic_chemical_shift=isotropic_chemical_shift,
            shielding_symmetric=shielding_symmetric,
        )
    else:
        site_instance = Site(
            name=f"site-{isotope}-{key}" ,
            isotope=isotope,
            isotropic_chemical_shift=isotropic_chemical_shift,
            shielding_symmetric=shielding_symmetric,
            quadrupolar=quadrupolar
        )

    return site_instance


# Dipolar and J-Coupling generator
def add_couplings_DJ(num_sites_fn , key):
    st.subheader(f"Define Advanced Coupling for Interaction Pair {key + 1}")

    col1, col2 = st.columns(2)
    with col1:
        spin1 = st.selectbox(
            "Select Site 1 (Index)", options=list( range( num_sites_fn ) ), key=f"spin1_{key}", index = 0,
        )
    with col2:
        spin2 = st.selectbox(
            "Select Site 2 (Index)", options=list( range( num_sites_fn ) ), key=f"spin2_{key}", index=1,
        )
    st.markdown("##### Define J-Coupling Parameters (Optional)")
    isotropic_j = st.number_input(
        "Isotropic J-Coupling (Hz)", value=0.0, format="%.3f", key=f"isotropic_j_{key}"
    )


    # Additional parameters for Dipolar coupling
    st.markdown("##### Define Dipolar Tensor Parameters (Optional)")
    col_ds1, col_ds2, col_ds3, col_ds4 = st.columns(4)
    with col_ds1:
        D_dipolar = st.number_input("D (Hz)", value=0.0, format="%.3f", key=f"D_dipolar_{key}")

    with col_ds2:
        alpha_d = st.number_input("α (degrees)", min_value=0.0, max_value=360.0, value=0.0, format="%.3f",
                                  key=f"alpha_d_{key}")
    with col_ds3:
        beta_d = st.number_input("β (degrees)", min_value=0.0, max_value=360.0, value=0.0, format="%.3f",
                                 key=f"beta_d_{key}")
    with col_ds4:
        gamma_d = st.number_input("γ (degrees)", min_value=0.0, max_value=360.0, value=0.0, format="%.3f",
                                  key=f"gamma_d_{key}")

    dipolar_tensor = SymmetricTensor(
        D=D_dipolar, alpha=alpha_d*np.pi/180.0, beta=beta_d*np.pi/180.0, gamma=gamma_d*np.pi/180.0
    )

    # Create the Coupling instance
    coupling_instance = Coupling(
        site_index=[spin1, spin2],
        isotropic_j=isotropic_j,
        dipolar=dipolar_tensor,
    )

    return coupling_instance

def get_sites(table_of_nuclei_fn):
    sites_fn = [ ]
    add_nuclei = "yes"
    num_site = 0
    while add_nuclei == "yes":
        site_added = add_site( table_of_nuclei_fn , num_site )
        sites_fn.append( site_added )
        num_site += 1
        st.write(f"Number of sites added: {num_site}")
        add_nuclei = st.selectbox("Add more sites?", ["yes", "no"], index = None, key = f'add_nuclei_{num_site}')
        st.divider()
    return sites_fn, num_site


def add_advanced_couplings(sites_coup_fn) :
    """
    Function to add J- or Dipolar Coupling between spin sites.

    Args:
        sites_coup_fn (list): List of spin sites; the length defines the number of available sites.

    Returns:
        list: A list of dipolar and J couplings added by the user.
    """
    dipolar_and_J_couplings = [ ]

    # Check if the number of sites is 2 or greater
    if len ( sites_coup_fn ) >= 2 :
        st.subheader ( "Add J- or Dipolar Coupling Between Spin Sites" )

        add_couplings_DJ_pairs = "yes"
        num_advanced_coupling = 0

        # Loop to keep adding couplings as long as the user chooses "yes"
        while add_couplings_DJ_pairs == "yes" :
            # Function for adding advanced coupling logic (assumes add_couplings_DJ exists)
            advanced_coupling_added = add_couplings_DJ ( len ( sites_coup_fn ) , num_advanced_coupling )
            dipolar_and_J_couplings.append ( advanced_coupling_added )
            num_advanced_coupling += 1

            # Inform user about the number of couplings defined so far
            st.write ( f"Number of advanced couplings defined: {num_advanced_coupling}" )

            # Ask user if they would like to add more couplings (select box)
            add_couplings_DJ_pairs = st.selectbox (
                "Add more couplings?" , [ "yes" , "no" ] , index=None , key=f"add_couplings_DJ_{num_advanced_coupling}"
            )
            st.divider ()

    return dipolar_and_J_couplings



def get_spin_system_parameters_from_gui_input_fn () :
    script_dir = os.path.dirname(__file__)
    csv_file = os.path.join(script_dir, '../resources/NMR_freq_table.csv')

    table_of_nuclei = pd.read_csv(csv_file)
    st.header("**Add nuclei**")

    sites, num_sites = get_sites(table_of_nuclei)
    if num_sites >= 2:
        dipolar_and_J_couplings = add_advanced_couplings(sites)


    st.subheader("Summary of Spin System Parameters")
    for idx, site in enumerate(sites, start=1):
        st.write(f"**Site {idx}:**")
        st.write(f"- **Nucleus**: {site.name} \t -- **σ (ppm)**: {site.isotropic_chemical_shift} ")

        # Display Shielding Symmetric Tensor Parameters
        if site.shielding_symmetric is not None:
            st.write("**Shielding Tensor Parameters:**")
            st.write(f"   **ζ (ppm)**: {site.shielding_symmetric.zeta} -- **η**: {site.shielding_symmetric.eta} -- **α (deg)**: {np.degrees(site.shielding_symmetric.alpha)} -- **β (deg)**: {np.degrees(site.shielding_symmetric.beta)} -- **γ (deg)**: {np.degrees(site.shielding_symmetric.gamma)}")

        # Display Quadrupolar Tensor Parameters if available
        if site.quadrupolar is not None:
            st.write("**Quadrupolar Tensor Parameters:**")
            st.markdown(f"  **Cq (Hz):** {site.quadrupolar.Cq} -- **η**: {site.quadrupolar.eta}  -- **α (deg)**: {np.degrees(site.quadrupolar.alpha)} -- **β (deg)**: {np.degrees(site.quadrupolar.beta)} -- **γ (deg)**: {np.degrees(site.quadrupolar.gamma)}")

        st.divider()  # Horizontal line for better readability


    # Showing Summary of Advanced Couplings
    st.subheader("Summary of Couplings")
    if num_sites >= 2:
        for idx, coupling in enumerate(dipolar_and_J_couplings, start=1):
            st.write(f"**Advanced Coupling {idx}:**")
            st.write(f"- **Site 1 Index**: {sites[coupling.site_index[0]].name} \t -- **Site 2 Index**: {sites[coupling.site_index[1]].name}")
            st.write("--")
            st.write(f"- **Isotropic J-Coupling (Hz)**: {coupling.isotropic_j}")
            st.write(f"- **Dipolar Coupling (Hz)**: {coupling.dipolar.D}")
            st.divider()
    else:
        st.write("No advanced couplings have been defined yet.")

    if num_sites >= 2:
        return SpinSystem(sites=sites, couplings=dipolar_and_J_couplings)
    else:
        return SpinSystem(sites=sites)

def common_1d_experimental_parameters(type_of_method) :
    # Introduce constants for readability
    ROTOR_ANGLE_CONVERSION = np.pi / 180
    ROTOR_FREQ_CONVERSION = 1000.0

    # Renamed and extracted user inputs for clarity
    experiment_channels = st.text_input ( "Channels (comma-separated)" , "1H" )
    rotor_frequency_khz = st.number_input ( "Rotor Frequency (kHz)" , min_value=0.0 , max_value=200.0 , value=10.0 )
    rotor_angle_degrees = st.number_input ( "Rotor Angle (degrees)" , min_value=0.0 , max_value=90.0 , value=54.735 )
    magnetic_flux_density_t = st.number_input ( "Magnetic Flux Density (T)" , min_value=0.0 , max_value=40.0 ,
                                                value=9.4 )
    spectral_dimension_count = st.number_input ( "Spectral Dimension Count" , min_value=1 , value=1024 )
    spectral_width_hz = st.number_input ( "Spectral Width (Hz)" , min_value=100.0 , value=25000.0 )
    spectral_dimension_label = st.text_input ( "Spectral Dimension Label" , "1D Spectrum" )
    reference_offset_hz = st.number_input ( "Reference Offset (Hz)" , value=0.0 )

    # Extracted reusable method generation logic
    def create_spectrum(method_class) :
        return method_class (
            channels=[experiment_channels] ,
            rotor_frequency=rotor_frequency_khz * ROTOR_FREQ_CONVERSION ,
            rotor_angle=rotor_angle_degrees * ROTOR_ANGLE_CONVERSION ,
            magnetic_flux_density=magnetic_flux_density_t ,
            spectral_dimensions=[
                SpectralDimension (
                    count=int ( spectral_dimension_count ) ,
                    spectral_width=spectral_width_hz ,
                    reference_offset=reference_offset_hz ,
                    label=spectral_dimension_label,
                )
            ] ,
        )

    # Modified method generator logic with simplified mapping
    method_classes = {
        "BlochDecaySpectrum" : BlochDecaySpectrum ,
        "BlochDecayCTSpectrum" : BlochDecayCTSpectrum ,
    }

    if type_of_method in method_classes :
        return create_spectrum ( method_classes[ type_of_method ] )
    else :
        raise ValueError ( f"Method {type_of_method} not found." )


def common_2d_experimental_parameters(type_of_method) :

    ROTOR_ANGLE_CONVERSION = np.pi / 180
    ROTOR_FREQ_CONVERSION = 1000.0

    # Renamed and extracted user inputs for clarity
    experiment_channels = st.text_input ( "Channels (comma-separated)" , "1H" )
    rotor_frequency_khz = st.number_input ( "Rotor Frequency (kHz)" , min_value=0.0 , max_value=200.0 , value=10.0 )
    rotor_angle_degrees = st.number_input ( "Rotor Angle (degrees)" , min_value=0.0 , max_value=90.0 , value=54.735 )
    magnetic_flux_density_t = st.number_input ( "Magnetic Flux Density (T)" , min_value=0.0 , max_value=40.0 ,
                                                value=9.4 )
    spectral_dimension1_count = st.number_input ( "Spectral Dimension Count" , min_value=1 , value=1024 , key="td1")
    spectral_width1_hz = st.number_input ( "Spectral Width (Hz)" , min_value=100.0 , value=25000.0 , key="sw1")
    reference_offset1_hz = st.number_input ( "Reference Offset (Hz)" , value=0.0, key="r1" )

    spectral_dimension2_count = st.number_input ( "Spectral Dimension Count" , min_value=1 , value=1024, key="td2" )
    spectral_width2_hz = st.number_input ( "Spectral Width (Hz)" , min_value=100.0 , value=25000.0, key="sw2" )
    reference_offset2_hz = st.number_input ( "Reference Offset (Hz)" , value=0.0, key="r2" )

    def create_spectrum_VAS(method_class) :
        return method_class (
            channels=[experiment_channels] ,
            magnetic_flux_density=magnetic_flux_density_t ,
            spectral_dimensions=[
                SpectralDimension (
                    count=int ( spectral_dimension1_count ) ,
                    spectral_width=spectral_width1_hz ,
                    reference_offset=reference_offset1_hz ,
                ),
                SpectralDimension (
                    count=int ( spectral_dimension2_count ) ,
                    spectral_width=spectral_width2_hz ,
                    reference_offset=reference_offset2_hz ,
                ),
            ]
        )

    def create_spectrum_SSB2D(method_class) :
        return method_class (
            channels=experiment_channels ,
            rotor_frequency=rotor_frequency_khz * ROTOR_FREQ_CONVERSION ,
            rotor_angle=rotor_angle_degrees * ROTOR_ANGLE_CONVERSION ,
            spectral_dimensions=[
                SpectralDimension (
                    count=int ( spectral_dimension1_count ) ,
                    spectral_width=spectral_width1_hz ,
                    reference_offset=reference_offset1_hz ,
                ),
                SpectralDimension (
                    count=int ( spectral_dimension2_count ) ,
                    spectral_width=spectral_width2_hz ,
                    reference_offset=reference_offset2_hz ,
                ),

            ]
        )

    method_classes = {
        "ThreeQ_VAS" : ThreeQ_VAS ,
        "FiveQ_VAS" : FiveQ_VAS ,
        "SevenQ_VAS" : SevenQ_VAS ,
        "ST1_VAS" : ST1_VAS ,
        "ST2_VAS" : ST2_VAS ,
        "SSB2D" : SSB2D ,
    }

    if type_of_method in method_classes :
        return create_spectrum_VAS ( method_classes[ type_of_method ] )
    elif type_of_method == "SS2BD" :
        return create_spectrum_SSB2D ( SSB2D )
    else :
        raise ValueError ( f"Method {type_of_method} not found." )



def get_method_parameters_from_gui_input_fn():
    st.title("Methods and Parameters")
    options_list = ["BlochDecaySpectrum", "BlochDecayCTSpectrum", "ThreeQ_VAS", "FiveQ_VAS",
                    "SevenQ_VAS", "ST1_VAS", "ST2_VAS", "SS2BD"]

    kind_of_spectrum = st.selectbox("What kind of spectrum do you want to simulate?", options_list, index=0)

    if kind_of_spectrum == "BlochDecaySpectrum" or kind_of_spectrum == "BlochDecayCTSpectrum":
        st.write("## Experimental Parameters")
        st.write("### 1D Spectrum")
        method = common_1d_experimental_parameters(kind_of_spectrum)
    else:
        st.write("## Experimental Parameters")
        st.write("### 2D Spectrum")
        method = common_2d_experimental_parameters(kind_of_spectrum)

    return method

def click_button():
    st.session_state.clicked = True


st.title("Mrsimulator GUI")

spin_system_tab, method_tab, simulation_tab, process_and_plot = st.tabs(["Spin System", "Method", "Simulation", "Plot Data"])

with spin_system_tab:
    spin_system = get_spin_system_parameters_from_gui_input_fn ()
with method_tab:
    method = get_method_parameters_from_gui_input_fn()

with simulation_tab:
    sim = Simulator(spin_systems=[spin_system], methods=[method])

    sim.config.integration_volume = st.selectbox("Average over: ", ["octant", "hemisphere", "sphere"], index=0)
    sim.config.integration_density = st.number_input("Integration Density: ", min_value=1, value=70, format="%d")
    sim.config.number_of_sidebands = st.number_input("Number of sidebands: ", min_value=1, value=64, format="%d")
    sim.config.number_of_gamma_angles = st.number_input("Number of gamma angles: ", min_value=1, value=32, format="%d")
    sim.config.decompose_spectrum = st.selectbox("Decompose Spectrum: ", ["none", "spin_system"], index=0)
    sim.config.isotropic_interpolation=st.selectbox("Isotropic Interpolation: ", ["linear", "gaussian"], index=0)

    if 'clicked' not in st.session_state :
        st.session_state.clicked = False

    st.button ( 'Run Simulation' , on_click=click_button )

    if st.session_state.clicked :
        with st.spinner("Running Simulation..."):
            sim.run()
            time.sleep ( 1 )
            st.success("Simulation Completed!")

with process_and_plot:
    plot_data_choice = st.selectbox ( "Process and Plot Data:" , [ "No" , "Yes" ] , index=None , key="plot_data" )
    if plot_data_choice == "Yes" :
        one_or_two_dim = st.selectbox("One or two dimensional data?", ["1D", "2D"], index=None, key="plot_data_dim")
        if one_or_two_dim == "1D":
            processed_dataset = None
            min_line_broadening_hz = st.number_input("Min Line Broadening in Hz:", value=0.0)
            max_line_broadening_hz = st.number_input("Max Line Broadening in Hz:", value=100.0)
            line_broadening_hz = st.slider( "Line Broadening in Hz:" , min_value=min_line_broadening_hz, max_value=max_line_broadening_hz , value=10.0, format="%f", step=2.0 )
            processor = sp.SignalProcessor(
                operations=[
                    sp.IFFT(),
                    sp.apodization.Gaussian( FWHM=f"{line_broadening_hz} Hz" ),
                    sp.FFT(),
                ]
            )
            processed_dataset = processor.apply_operations(dataset=sim.methods[0].simulation)

            fig = go.Figure()
            hz_or_ppm = st.selectbox("Axis in Hz or ppm?", ["Hz", "ppm"], index=None, key="x_scale_choice")
            if hz_or_ppm == "Hz":
                x_oned = sim.methods[0].spectral_dimensions[0].coordinates_Hz()
            else:
                x_oned = sim.methods[0].spectral_dimensions[0].coordinates_ppm()
            # Add trace for real part of the dataset
            y_data = np.array(processed_dataset.real.dependent_variables[0].components[0])

            fig.add_trace(go.Scatter(
                x=x_oned,  # Inverting x-axis
                y=y_data,
                mode='lines',
                line=dict(color='rgb(22, 128, 178)', width=2)
            ))

            # Update layout
            fig.update_layout(
                width=600,
                height=400,
                xaxis_title=f"Chemical Shift / {hz_or_ppm}",
                yaxis_title="Arbitrary Units",
                template="plotly_white"
            )

            fig.update_xaxes(autorange="reversed")

            # Display in Streamlit
            st.plotly_chart(fig)

        elif one_or_two_dim == "2D":
            processed_dataset = None
            min_line_broadening_hz_dim1 = st.number_input("Min Line Broadening in Hz:", value=0.0, key='lb11')
            max_line_broadening_hz_dim1 = st.number_input("Max Line Broadening in Hz:", value=100.0, key='lb12')
            min_line_broadening_hz_dim2 = st.number_input ( "Min Line Broadening in Hz:" , value=0.0, key='lb13' )
            max_line_broadening_hz_dim2 = st.number_input ( "Max Line Broadening in Hz:" , value=100.0 , key='lb14' )

            line_broadening_hz_dim1 = st.slider( "Line Broadening in Hz:" , min_value=min_line_broadening_hz_dim1, max_value=max_line_broadening_hz_dim1 , value=10.0, format="%f", step=2.0 , key='lb15' )
            line_broadening_hz_dim2 = st.slider( "Line Broadening in Hz:" , min_value=min_line_broadening_hz_dim2, max_value=max_line_broadening_hz_dim2 , value=10.0, format="%f", step=2.0 , key='lb16' )

            processor = sp.SignalProcessor (
                operations=[
                    # Gaussian convolution along both dimensions.
                    sp.IFFT ( dim_index=(0 , 1) ) ,
                    sp.apodization.Gaussian ( FWHM=f"{line_broadening_hz_dim1} Hz" , dim_index=0 ) ,
                    sp.apodization.Gaussian ( FWHM=f"{line_broadening_hz_dim2} Hz" , dim_index=1 ) ,
                    sp.FFT ( dim_index=(0 , 1) ) ,
                ]
            )
            processed_dataset = processor.apply_operations ( dataset=sim.methods[ 0 ].simulation )
            processed_dataset /= processed_dataset.max () * 100.0

            fig = go.Figure()
            hz_or_ppm_d1 = st.selectbox("Axis in Hz or ppm in F1 dim?", ["Hz", "ppm"], index=None, key="x_scale_choice")
            hz_or_ppm_d2 = st.selectbox("Axis in Hz or ppm in F2 dim?", ["Hz", "ppm"], index=None, key="y_scale_choice")


            if hz_or_ppm_d1 == "Hz":
                x_scale = sim.methods[0].spectral_dimensions[0].coordinates_Hz()
            else:
                x_scale = sim.methods[0].spectral_dimensions[0].coordinates_ppm()

            if hz_or_ppm_d2 == "Hz":
                y_scale = sim.methods[0].spectral_dimensions[1].coordinates_Hz()
            else:
                y_scale = sim.methods[0].spectral_dimensions[1].coordinates_ppm()

            z_data = np.array(processed_dataset.real.dependent_variables)

            # data = np.array(processed_dataset.real.dependent_variables[0].components[0])
            st.write(z_data)
            fig = go.Figure(data=go.Contour(
                x=x_scale,
                y=y_scale,
                z=z_data,
                colorscale="viridis",
                contours=dict(
                    start=0.,
                    end=100.0,
                    size=50,
                ),
            ))

            fig.update_layout(
                width=425,
                height=300
            )

            st.plotly_chart(fig)

