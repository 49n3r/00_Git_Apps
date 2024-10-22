import streamlit as st
import random as ran
import numpy as np
import pandas as pd
import time


### Set up Page ###########
st.set_page_config(page_title="Amorphous Graphite Initiator", page_icon="https://chinonsougwumadu.com/wp-content/uploads/2024/05/microsoftteams-image-17.jpg")
st.markdown("# Amorphous Graphite Constructor")

logo_url = "https://chinonsougwumadu.com/wp-content/uploads/2024/05/microsoftteams-image-17.jpg"
st.sidebar.image(logo_url)
#st.sidebar.link_button("About Us", "https://daviddrabold.com/")
st.sidebar.markdown("## Paramater initialization")
st.markdown(
    """
    This app generates and distributes carbon atoms randomly in 3D space while maintaining periodic boundary conditions. 
    It creates initial configurations for molecular dynamics simulations of amorphous graphite. 
    You can download the POSCAR file directly from the app!

    > Associated Publications:
    - [Ab initio simulation of amorphous graphite](http://people.ohio.edu/drabold/pubs/246.pdf), Phys. Rev. Lett., 128 236402 (2022).
    - [Atomistic nature of amorphous graphite](http://people.ohio.edu/drabold/pubs/250.pdf), Euro. J. Glass Sci. and Tech.: B, 64 16 (2023).
    
    """
)



### To disable the generate key until the pore distribution button is clicked on
def disable(b):
    st.session_state["disabled"] = b


## Initialize Data
st.sidebar.number_input("Number of Atoms:", min_value=60, value=1000, step = 50, key="num_atoms",on_change=disable, args=(False,), help="The number of C atoms required")
st.sidebar.number_input("Density [g/cm$^3$]:", min_value=2.0, max_value = 4.0, value = 2.44, step = 0.02, key="density",on_change=disable, args=(False,), help ="The desired density of the model")
st.sidebar.slider("Carbon Bonds Initial Cutoff [\u212B]", min_value=1.0,max_value=1.4, step = 0.1, value =1.2, key='cutoff',on_change=disable, args=(False,), help="C-C cutoff, 1.2 \u212B is a good choice")


# You can access the value at any point with:
# st.session_state.num_atoms
# st.session_state.density
# st.session_state.cutoff

st.session_state.num_atoms_arr = np.array([st.session_state.num_atoms], dtype = np.int64)

### Naming convention for saved files ##########################################################
stringDensity = str(np.round(st.session_state.density,2)).replace(".","p")+"gcc_"
stringNumAtoms = str(st.session_state.num_atoms)+"atoms_"


st.session_state.atoms_vasp = "POSCAR_"+stringNumAtoms+stringDensity
#################################################################################################

##################### FUNCTIONS ##############################################################################################

#def boxSize(*,density=st.session_state.density):
def boxSize(density):
    '''
    This function predicts the box size for the model, in units of angstrom
    '''

    atom_mass_amu = np.array([12.0107])
    atom_mass_gram = [atom*1.66054e-24 for atom in atom_mass_amu]  #Convert amu to grams (1 amu = 1.66054e-24 g)

    total_mass = sum([atom_mass_gram[i]*st.session_state.num_atoms_arr[i] for i in range(np.size(atom_mass_gram))])
    volume = total_mass / density

    box_cm = volume**(1/3) #(cm)
    ##convert box length  in cm to armstrong
    box_arm = box_cm/1e-8
    return box_arm, volume

# Some Parameter Initialization
st.session_state.box, volume = boxSize(st.session_state.density)
st.session_state.final_output = ""

################################### MAIN ################################################################

col1, col2 = st.columns(2)      
        
with col1:
    st.markdown("##### Generate Amorphous Graphite Model")

    if st.button('Generate Model', key="generateButton",on_click=disable, args=(False,)):
        with col1:
            
            st.write(f"The box length for {st.session_state.num_atoms} C atoms is {st.session_state.box:.2f} \u212B")


        ######################## Amorphous C constructor Algorithm starts here ####################################################

        pos = np.zeros([st.session_state.num_atoms,3],float)   ## A list that takes in the position of the atoms
        
        box = st.session_state.box
        num_atoms = st.session_state.num_atoms
        cutoff = st.session_state.cutoff

        ct = 0
        rnd = lambda i: i-round(i/box)*box  ## This makes rnd a function that takes i and perform i - round(i/box)*box
        vec_rnd = np.vectorize(rnd)         ## takes  a function and gives result in a callable vectorised function
 
     
        with st.spinner("Creating carbon atoms. Please wait..."):
            my_bar = st.progress(0, text="Progress Status.")

            while ct < num_atoms:
                atoms =[box*ran.random(),box*ran.random(),box*ran.random()]
                test = 0
    
                ### Uncomment for debugging
                #print (atoms-pos[test][:])
                ###########################
    
                while test < ct:
                    atom_distance = atoms-pos[test][:]
                    atom_distance = vec_rnd(atom_distance)
                   
                    ### Position the atoms if atoms are not close to center    
                    if sum(map(lambda i: i*i, atom_distance)) < cutoff**2:
                        atoms =np.array([box*ran.random(),box*ran.random(),box*ran.random()])
                        
                        ### Uncomment for debugging##################
                        #print ("OOPS!!! TOO CLOSE to another atom")
                        #############################################
            
                        test = 0
                    else:
                        test += 1
            
                pos[ct][:] = atoms
    
    
                ct +=1
                my_bar.progress(ct/(num_atoms), text="Progress Status.")         
                print (f"Placing Atom number {ct} of {num_atoms}", end='\r')
         ################################################################################################

        my_bar.empty()
        with st.spinner("Preparing file for download. Please wait..."):
            time.sleep(5)

    
        atom_position = pos/st.session_state.box

        #Convert the NumPy array to a list of lists
        list_of_lists = atom_position.tolist()

        #Convert each sublist to a string with elements separated by spaces
        lines = [" ".join(map(str, sublist)) for sublist in list_of_lists]

        #Join the resulting lines with newline characters
        st.session_state.final_output = "\n".join(lines)

with col2:
    st.header("")
    #my_anime = st.markdown("![Alt Text](https://chinonsougwumadu.com/wp-content/uploads/2024/06/ag_anime-2.gif?w=1024)")
    my_anime =st.markdown(
    f'<img src="https://chinonsougwumadu.com/wp-content/uploads/2024/06/ag_anime-2.gif?w=1024" width="500" alt="Amorphous Graphite Animation">',
    unsafe_allow_html=True,)

    if st.session_state.generateButton:
        st.session_state.txt1 = st.text_area("POSCAR File", f"{stringNumAtoms} {stringDensity}\n\
{st.session_state.box:10.6f}\n\
{1.0:2.6f} {0.0:2.6f} {0.0:2.6f}\n\
{0.0:2.6f} {1.0:2.6f} {0.0:2.6f}\n\
{0.0:2.6f} {0.0:2.6f} {1.0:2.6f}\n\
C \n\
{st.session_state.num_atoms} \n\
Direct\n{st.session_state.final_output}")
                                               
        my_anime.empty()
        st.success('Done!')
           
        if st.download_button(label="Download POSCAR",key="downloadPOSCAR",data=st.session_state.txt1, file_name=st.session_state.atoms_vasp,disabled=st.session_state.get("disabled", True)):
            st.write("Download Complete.")
