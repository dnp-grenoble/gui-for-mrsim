import streamlit as st

def app():
    st.title("MRSimulator GUI")
    st.markdown("This is a simple GUI for the MRSimulator package.")
    st.markdown("This is still in development. Please report any issues to the GitHub repository.")
    st.markdown("---")
    st.markdown("## About")
    st.markdown("MRSimulator is a Python package for simulating magnetic resonance (MR) spectra.")
    st.markdown("MRSimulator is developed by the [MRSimulator Team](https://mrsimulator.readthedocs.io/en/stable/).")

    st.markdown('''
    *Helpful Features*


    - You can enter nuclei information in a table format. 
    - You can add couplings between sites. 
    - It delineates the options needed for each method.
    - Helps to generate 1D and 2D spectra.''')

    st.markdown("---")

    st.markdown("## Features to be added")

    st.markdown('''
    - Add special methods like simulation of INADEQUATE spectra.
    - Add custom methods.
    
    For complex simulation and learn how it functions please read:  
    Srivastava, Deepansh J., and Philip J. Grandinetti. ‘Simulating Multipulse NMR Spectra of Polycrystalline Solids in the Frequency Domain’. The Journal of Chemical Physics 160, no. 23 (21 June 2024): 234110. https://doi.org/10.1063/5.0209887.
    
    If you want to run the package off line, please follow the steps below:
    ## Installation and Usage

        ### Prerequisites
        
        1. Install **[Streamlit](https://docs.streamlit.io/get-started/installation)**.
        2. Clone the repository to your local machine:
           ```bash
           git clone https://github.com/dnp-grenoble/simpson_gui.git
           ```
        3. Install dependencies:
           ```bash
           pip install -r requirements.txt
           ```
        
        ### Running the Application
        
        1. Navigate to the project directory.
        2. Launch the application using **Streamlit**:
           ```bash
           streamlit run Homepage.py
    
    ''')

if __name__ == "__main__":
    app()
