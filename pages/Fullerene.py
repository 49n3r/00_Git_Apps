import streamlit as st
import random as ran
import numpy as np
import pandas as pd
import time


### Set up Page ###########
st.set_page_config(page_title="Multi-shell Fullerene", page_icon="https://chinonsougwumadu.com/wp-content/uploads/2024/05/microsoftteams-image-17.jpg")
st.markdown("# Multi-shell Fullerene Constructor")

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
    - [Simulation of multi-shell fullerenes using Machine-Learning Gaussian Approximation Potential](https://doi.org/10.1016/j.cartre.2022.100239), Carbon Trends, 10 100239 (2023).
    
    """
)



### To disable the generate key until the pore distribution button is clicked on
def disable(b):
    st.session_state["disabled"] = b


## Initialize Data
st.sidebar.number_input("Number of Atoms:", min_value=60, value=1000, step = 50, key="num_atoms",on_change=disable, args=(False,), help="The number of C atoms required")
st.sidebar.number_input("Density [g/cm$^3$]:", min_value=1.4, max_value = 2.8, value = 2.26, step = 0.02, key="density",on_change=disable, args=(False,), help ="The density of the fullerene model")
st.sidebar.slider("Carbon Bonds Initial Cutoff [\u212B]", min_value=1.0,max_value=1.4, step = 0.1, value =1.2, key='cutoff',on_change=disable, args=(False,), help="C-C cutoff, 1.2 \u212B is a good choice")
st.sidebar.slider("3D Vaccum [\u212B]", min_value=3.0,max_value=8.0, step = 0.2, value =6.0, key='vacuum',on_change=disable, args=(False,), help="Vacuum added to all 3 dimensions.\n 6 \u212B is a good choice")

# You can access the value at any point with:

num_atoms = st.session_state.num_atoms
density = st.session_state.density
cutoff = st.session_state.cutoff
vacuum = st.session_state.vacuum



st.session_state.num_atoms_arr = np.array([num_atoms], dtype = np.int64)

### Naming convention for saved files ##########################################################
stringDensity = str(np.round(density,2)).replace(".","p")+"gcc_"
stringNumAtoms = str(num_atoms)+"atoms_"


st.session_state.atoms_vasp = "POSCAR_"+stringNumAtoms+stringDensity
#################################################################################################

##################### FUNCTIONS ##############################################################################################


def getRadius(density):
    '''
    This function predicts the box size for the model, in units of angstrom
    '''

    atom_mass_amu = np.array([12.0107])
    atom_mass_gram = [atom*1.66054e-24 for atom in atom_mass_amu]  #Convert amu to grams (1 amu = 1.66054e-24 g)

    total_mass = sum([atom_mass_gram[i]*st.session_state.num_atoms_arr[i] for i in range(np.size(atom_mass_gram))])
    volume = total_mass / density

    ######### For spherical radius #########
    radius_cm = np.cbrt((3*volume) / (4*np.pi))
    
    ##convert box lenght  in cm to armstrong
    
    radius_angs = radius_cm/1e-8
    vol_sphere_angs = volume/(1e-8)**3
    return radius_angs,vol_sphere_angs

# Some Parameter Initialization
st.session_state.radius, _ = getRadius(density)
box = st.session_state.radius
radius = st.session_state.radius
st.session_state.big_box = 2*radius + vacuum
big_box = st.session_state.big_box
st.session_state.final_output = ""

def isInSphere(_x,_y,_z,box):
    '''
    Takes a specified co-ordinate and determine if the points fall within the spherical constrains
    Returns one that works if the answer is FALSE
    '''
    center = big_box/2
    ## if axis not within circle reject it
    while (_x - center)**2 + (_y - center)**2 + (_z- center)**2 > radius**2:
        _x = big_box*ran.random()
        _y = big_box*ran.random()
        _z = big_box*ran.random()        
    return _x,_y,_z 


################################### MAIN ################################################################

col1, col2 = st.columns(2)      
        
with col1:
    st.markdown("##### Generate Fullerene Model")
    st.write(f"The radius for {num_atoms} C atoms is {radius:.2f} \u212B")

    if st.button('Generate Model', key="generateButton",on_click=disable, args=(False,)):

        st.write(f"A {vacuum} \u212B vaccum will be added in all 3 dimensions.")           
            


        ######################## Amorphous C constructor Algorithm starts here ####################################################

        pos = np.zeros([num_atoms,3],float)   ## A list that takes in the position of the atoms
        
        ct = 0
        rnd = lambda i: i - round(i/(2 * radius)) * (2 * radius)  
        vec_rnd = np.vectorize(rnd)         ## takes  a function and gives result in a callable vectorised function
 
     
        with st.spinner("Creating carbon atoms. Please wait..."):
            my_bar = st.progress(0, text="Progress Status.")

            while ct < num_atoms:
                atoms =[box*ran.random(),box*ran.random(),box*ran.random()]
                atoms[0],atoms[1], atoms[2] = isInSphere(atoms[0],atoms[1],atoms[2],box)
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
                        atoms[0],atoms[1], atoms[2] = isInSphere(atoms[0],atoms[1],atoms[2],box)
                        
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

    
        atom_position = pos/big_box

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
    f'<img src="https://chinonsougwumadu.com/wp-content/uploads/2024/06/afullerenes_anime.gif?w=1024" width="500" alt="Multi-shell Fullerene Animation">',
    unsafe_allow_html=True,)

    if st.session_state.generateButton:
        st.session_state.txt1 = st.text_area("POSCAR File", f"{stringNumAtoms} {stringDensity}\n\
{1.0:2.6f}\n\
{big_box:2.6f} {0.0:2.6f} {0.0:2.6f}\n\
{0.0:2.6f} {big_box:2.6f} {0.0:2.6f}\n\
{0.0:2.6f} {0.0:2.6f} {big_box:2.6f}\n\
C \n\
{num_atoms} \n\
Direct\n{st.session_state.final_output}")
                                               
        my_anime.empty()
        st.success('Done!')
           
        if st.download_button(label="Download POSCAR",key="downloadPOSCAR",data=st.session_state.txt1, file_name=st.session_state.atoms_vasp,disabled=st.session_state.get("disabled", True)):
            st.write("Download Complete.")