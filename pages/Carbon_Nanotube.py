import streamlit as st
import random as ran
import numpy as np
import pandas as pd
import time


st.set_page_config(page_title="Carbon Nanotube Initializer", page_icon="https://chinonsougwumadu.com/wp-content/uploads/2024/05/microsoftteams-image-17.jpg")
st.markdown("# Multi-walled CNT Constructor")

logo_url = "https://chinonsougwumadu.com/wp-content/uploads/2024/05/microsoftteams-image-17.jpg"
st.sidebar.image(logo_url)
st.sidebar.markdown("## Paramater initialization")
st.markdown(
    """
    This app generates and distributes carbon atoms randomly in 3D space while maintaining periodic boundary conditions. 
    It creates initial configurations for molecular dynamics simulations of amorphous graphite. 
    You can download the POSCAR file directly from the app!

    > Associated Publications:
    - [Formation of amorphous carbon multi-walled nanotubes from random initial configurations](https://onlinelibrary.wiley.com/toc/15213951/2023/260/3), Phys. Status Solidi B, 260: 2200527. 
   
    """
)



### To disable the generate key until the pore distribution button is clicked on
def disable(b):
    st.session_state["disabled"] = b


## Initialize Data
st.sidebar.number_input("Number of Atoms:", min_value=60, value=1000, step = 50, key="num_atoms",on_change=disable, args=(False,), help="The number of C atoms required")
st.sidebar.slider("Aspect Ratio [height/diameter]:", min_value=1.0, max_value=5.5, value = 1.3, step = 0.05, key="aspect_ratio",on_change=disable, args=(False,), help ="The height/diameter of the CNT")

num_atoms = st.session_state.num_atoms
st.session_state.num_atoms_arr = np.array([num_atoms], dtype = np.int64)
aspect_ratio = st.session_state.aspect_ratio

## Required to define the radius, Density can be between 0.15 - 2 g/cc
def getRadiusRange():
    '''
    Required to decide what Radius range should be
    '''
    atom_mass_amu = np.array([12.0107])
    atom_mass_gram = [atom*1.66054e-24 for atom in atom_mass_amu]  #Convert amu to grams (1 amu = 1.66054e-24 g)
    total_mass = sum([atom_mass_gram[i]*st.session_state.num_atoms_arr[i] for i in range(np.size(atom_mass_gram))])
    r_max = (total_mass/(0.15*aspect_ratio*2*np.pi))**(1/3) * 1e8 # convert cm to angs
    r_val = (total_mass/(1.7*aspect_ratio*2*np.pi))**(1/3) * 1e8 # convert cm to angs
    r_min = (total_mass/(2.1*aspect_ratio*2*np.pi))**(1/3)    * 1e8 # convert cm to angs
    return r_min,r_val, r_max

r_min, r_val, r_max = getRadiusRange()

st.sidebar.slider("Radius [cm]:", min_value=r_min, value = r_val, max_value =r_max, step = 0.05, key="radius",on_change=disable, args=(False,), help ="The radius of the CNT")
st.sidebar.slider("Carbon Bonds Initial Cutoff [\u212B]", min_value=1.0,max_value=1.4, step = 0.1, value =1.2, key='cutoff',on_change=disable, args=(False,), help="C-C cutoff, 1.2 \u212B is a good choice")
st.sidebar.slider("XY Plane Vacuum [\u212B]", min_value=3.0,max_value=8.0, step = 0.2, value =6.0, key='vacuum',on_change=disable, args=(False,), help="Vacuum added to XY plane.\n3\u212B is a good choice")

# You can access the value at any point with:

radius = st.session_state.radius
cutoff = st.session_state.cutoff
vacuum = st.session_state.vacuum



### Naming convention for saved files ##########################################################
stringRadius = str(np.round(radius,1)).replace(".","p")+"A_"
stringNumAtoms = str(num_atoms)+"atoms_"
stringAspectRatio = str(aspect_ratio).replace(".","p")+"aspect_ratio"


st.session_state.atoms_vasp = "POSCAR_"+stringNumAtoms+stringRadius+stringAspectRatio
#################################################################################################

##################### FUNCTIONS ##############################################################################################


def cntHeight(radius, aspect_ratio):
    '''
    This function predicts the height of the cylinder 
    '''
    h_angs = 2*radius * aspect_ratio
    r_cm = radius * 1e-8
    h_cm = h_angs * 1e-8
    volume = np.pi * r_cm**2 * h_cm #in cm
    atom_mass_amu = np.array([12.0107])
    atom_mass_gram = [atom*1.66054e-24 for atom in atom_mass_amu]  #Convert amu to grams (1 amu = 1.66054e-24 g)
    total_mass = sum([atom_mass_gram[i]*st.session_state.num_atoms_arr[i] for i in range(np.size(atom_mass_gram))])
    density = total_mass/volume # in g/cc

    return h_angs, density


def isInCylinder(_x,_y,box):
    '''
    Takes a specified co-ordinate and determine if the points fall within the spherical constrains
    Returns one that works if the answer is FALSE
    '''
    center = box/2
    ## if axis not within circle reject it
    while (_x - center)**2 + (_y - center)**2  > radius**2:
        _x = box*ran.random()
        _y = box*ran.random()
       #_z = box*ran.random()        
    return _x,_y 

# Some Parameter Initialization
st.session_state.height, st.session_state.density = cntHeight(radius,aspect_ratio)
st.session_state.final_output = ""
box = st.session_state.height
height = st.session_state.height
density = st.session_state.density

################################### MAIN ################################################################

col1, col2 = st.columns(2)      
        
with col1:
    st.markdown("##### Generate Multi-Walled CNT")
    st.write(f"The denisty of the CNT (no vacuum) is {density:.2f} g/cm$^3$.")
            


    if st.button('Generate Model', key="generateButton",on_click=disable, args=(False,)):          
           
        st.write(f"A {vacuum} \u212B vaccum will be added in xy plane")
     
        ######################## Amorphous C constructor Algorithm starts here ####################################################

        pos = np.zeros([num_atoms,3],float)   ## A list that takes in the position of the atoms

        ct = 0
        rnd_xy = lambda i: i - round(i/(2 * radius)) * (2 * radius)
        rnd_z = lambda i: i-round(i/height)*height  ## This makes rnd a function that takes i and perform i - round(i/box)*box
        vec_rnd_xy = np.vectorize(rnd_xy)
        vec_rnd_z = np.vectorize(rnd_z)         ## takes  a function and gives result in a callable vectorised function
 
     
        with st.spinner("Creating carbon atoms. Please wait..."):
            my_bar = st.progress(0, text="Progress Status.")

            while ct < num_atoms:
                atoms =[box*ran.random(),box*ran.random(),box*ran.random()]
                atoms[0],atoms[1] = isInCylinder(atoms[0],atoms[1],box)
                test = 0
    
                ### Uncomment for debugging
                #print (atoms-pos[test][:])
                ###########################
    
                while test < ct:
                    atom_distance = atoms-pos[test][:]
                    atom_distance[0] = vec_rnd_xy(atom_distance[0]) # doing for x
                    atom_distance[1] = vec_rnd_xy(atom_distance[1]) # doing for x
                    atom_distance[2] = vec_rnd_z(atom_distance[2])    # doing for z
                   
                    ### Position the atoms if atoms are not close to center    
                    if sum(map(lambda i: i*i, atom_distance)) < cutoff**2:
                        atoms =np.array([box*ran.random(),box*ran.random(),box*ran.random()])
                        atoms[0],atoms[1] = isInCylinder(atoms[0],atoms[1],box)
                        
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

    
        atom_position = pos/height

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
    f'<img src="https://chinonsougwumadu.com/wp-content/uploads/2024/06/acnt_anime-1.gif?w=1024" width="500" alt="Multi-walled CNT Animation">',
    unsafe_allow_html=True,)

    if st.session_state.generateButton:
        st.session_state.txt1 = st.text_area("POSCAR File", f"{stringNumAtoms}atoms Radius: {stringRadius}\u212B Aspect-ratio: {stringAspectRatio}\n\
{1.0:10.6f}\n\
{(radius + vacuum):2.6f} {0.0:2.6f} {0.0:2.6f}\n\
{0.0:2.6f} {(radius + vacuum):2.6f} {0.0:2.6f}\n\
{0.0:2.6f} {0.0:2.6f} {box:2.6f}\n\
C \n\
{st.session_state.num_atoms} \n\
Direct\n{st.session_state.final_output}")
                                               
        my_anime.empty()
        st.success('Done!')
           
        if st.download_button(label="Download POSCAR",key="downloadPOSCAR",data=st.session_state.txt1, file_name=st.session_state.atoms_vasp,disabled=st.session_state.get("disabled", True)):
            st.write("Download Complete.")
